function outDir = test_ideal_point_bench_ca()
% Fixed, low-speed linear-stage fixture. No cone, classifier or FPGA model.
% Exactly one isotropic unit scatterer at the carriage centre. Simulate
% the active RF receive window; PRI is 1 ms (empty RF tail not simulated).
root = fileparts(fileparts(mfilename('fullpath')));
oldPath = path; restore = onCleanup(@() path(oldPath)); %#ok<NASGU>
addpath(fullfile(root,'functions')); ensure_md_zt_paths();
sigmaR = 0.0003; sigmaV = 0.001; seed = 17;
info = isar_ca_run_info(fullfile(root,'results','ideal_point_bench'), ...
    mfilename,"centre","measured",sigmaR,sigmaV);
outDir = char(info.outputDirectory);
fprintf('Output: %s\n',outDir);
try
    cfg = struct('fs',1e9,'fs_scope',200e6,'Tr',100e-6,'Ta',1e-3,'Na',32, ...
        'fL',11.8e9,'fH',17.8e9,'fc',14.8e9,'propagationSpeed',3e8,'kr',6e13);
    cfg.lambda = cfg.propagationSpeed/cfg.fc;
    cfg.T = 128e-6;
    p = isarpath_settings(cfg);
    p.phaseCompensationMode = "MODEL_CA";
    assert(p.rangeAlignmentMode=="INTEGER" && p.rangeMigrationMode=="FFT2048_SINGLE" && ...
        p.dispersionCompensationEnabled && p.windowType=="NO");
    epoch = (0:.1:.6)'; n = numel(epoch);
    % Truth generates RF and independent measurement observations, never
    % passed to the processing prior. A single 5 MHz bench LO is intentional.
    r0=20; v0=.22; a0=.6; bias=5e6;
    rec = isar_cone_measurement_record(epoch,r0+v0*epoch+0.5*a0*epoch.^2, ...
        v0+a0*epoch,struct('Mode',"measured",'RangeSigmaM',sigmaR, ...
        'VelocitySigmaMps',sigmaV,'Seed',seed));
    aFit = isar_ca_record_acceleration(rec);
    tf = (0:round(cfg.T*cfg.fs)-1)/cfg.fs;
    eo = struct('PulseWidthSec',cfg.Tr,'PropagationSpeed',cfg.propagationSpeed, ...
        'FastTimeChunkSize',32768);
    modes = ["NONE","MODEL_CA"]; results = cell(n,2); audits = cell(n,1);
    f = (-128:127)'/(256*cfg.Ta);
    theoryHz = .885893/(cfg.Na*cfg.Ta);
    save(fullfile(outDir,'config.mat'),'cfg','p','rec','info','seed','r0','v0','a0','bias','aFit');
    for k=1:n
        t = epoch(k)+(0:cfg.Na-1)'*cfg.Ta;
        prior = isar_cone_processing_prior(rec,k,aFit,bias,cfg);
        audits{k} = isar_ca_lo_audit(0,k,bias,prior);
        br = struct('fastTimeSec',tf,'referenceStartTimeSec',prior.loStartSec, ...
            'f0Hz',cfg.fL,'f1Hz',cfg.fH,'pulseWidthSec',cfg.Tr, ...
            'outputSampleRateHz',cfg.fs_scope);
        bo = struct('Mode','ANALYTIC_IF','RemoveDc',true,'ReturnPerTarget',false);
        H = complex(zeros(cfg.Na,p.hrrfBins));
        for j=1:cfg.Na
            echo = radar.echo_generate_stateFcn_lfm(@point_state,1,tf,t(j),cfg.fL,cfg.fH,eo);
            beat = radar.lfm_local_reference_beat({echo},br,bo);
            segment = isarpath_extract_segment(beat.total,cfg,prior,t(j),p);
            H(j,:) = isarpath_range_profiles(segment,cfg,prior,t(j),p);
        end
        [aligned,shifts] = isarpath_align_profiles(H,p);
        assert(all(isfinite(H(:))) && any(H(:)~=0) && all(shifts==round(shifts)) && ...
            all(abs(shifts)<p.rcmcSearchRadiusBins),'Invalid HRRF or RCMC search boundary.');
        for mi=1:2
            pp=p; pp.phaseCompensationMode=modes(mi); ref=[];
            if mi==2, ref=isarpath_prior_phase(cfg,prior,t,p,shifts,"MODEL_CA"); end
            [post,img,factor] = isarpath_phase_image(aligned,pp,ref);
            rd = fftshift(fft(post,256,1),1)/256;
            amp = abs(rd); [pk,ix]=max(amp(:)); [~,col]=ind2sub(size(rd),ix);
            assert(isequal(size(rd),[256 1024]) && pk>0 && all(isfinite(rd(:))));
            assert(max(abs(amp(:)/pk-img(:)))<1e-9 && max(abs(abs(factor)-1))<1e-9);
            % Denser frequency samples only for numerical half-power crossings.
            N=4096; df=1/(N*cfg.Ta); mid=N/2+1;
            line=abs(fftshift(fft(post(:,col),N))); [peak,row]=max(line);
            offset=(-N/2:N/2-1)'*df; peakHz=offset(row);
            line=circshift(line/peak,mid-row);
            assert(any(line(1:mid-1)<1/sqrt(2)) && any(line(mid+1:end)<1/sqrt(2)));
            width=signal.resolution_find_3db_width_linear(line,offset,mid);
            results{k,mi}=struct('frame',k,'epochSec',epoch(k),'mode',modes(mi), ...
                'peakCol1',col,'peakHz',peakHz,'fwhmHz',width,'idealHz',theoryHz, ...
                'minShift',min(shifts),'maxShift',max(shifts));
        end
        save(fullfile(outDir,sprintf('frame_%02d.mat',k)), ...
            'H','aligned','shifts','factor','rd','f','prior');
        fprintf('CPI %d/%d: CA peak %.3f Hz, FWHM %.3f Hz\n',k,n,peakHz,width);
    end
    metrics=struct2table(vertcat(results{:})); ca=metrics(metrics.mode=="MODEL_CA",:);
    writetable(metrics,fullfile(outDir,'frame_metrics.csv'));
    writetable(struct2table(vertcat(audits{:})),fullfile(outDir,'lo_per_frame.csv'));
    fig=figure('Visible','off','Color','w'); closeFig=onCleanup(@() close(fig)); %#ok<NASGU>
    cols=max(1,col-20):min(1024,col+20);
    imagesc(cols,f,20*log10(max(amp(:,cols)/pk,1e-8)));axis xy;clim([-40 0]);colorbar;
    xlabel('HRRF column (one-based)');ylabel('Doppler (Hz)');title('Ideal point, final CPI, MODEL CA');
    print(fig,fullfile(outDir,'point_rd.png'),'-dpng','-r140');
    assert(abs(aFit-a0)<0.05,'Bench range-fit error exceeds 0.05 m/s^2.');
    assert(all(abs(ca.peakHz)<2/(256*cfg.Ta)),'Residual peak exceeds two display bins.');
    assert(all(abs(ca.fwhmHz/theoryHz-1)<.10),'CA FWHM outside 10%% of the rectangular aperture.');
    save(fullfile(outDir,'PASS.mat'),'metrics','theoryHz','aFit');
    fprintf('PASS: 7 CPI, 224 pulse echoes; truth %.5f, fitted %.5f m/s^2\n',a0,aFit);
catch err
    failure=getReport(err,'extended','hyperlinks','off');
    save(fullfile(outDir,'FAIL.mat'),'failure'); rethrow(err);
end

    function [r,w]=point_state(tAbs,index)
        assert(isequal(index,1));
        u=tAbs(:).'; r=r0+v0*u+0.5*a0*u.^2; w=ones(size(r));
    end
end
