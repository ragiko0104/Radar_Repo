function result=run_isar_tip_slowtime_tracks(options)
%RUN_ISAR_TIP_SLOWTIME_TRACKS Range-isolated tip slices from full-cone RF.
% Mirror the field experiment's range-gated slow-time analysis, while
% retaining the current native ISARpath frontend and CA prior correction.
% Diagnostic point channels quantify leakage; they never select the gate or
% supply compensation. Only tip spectra, tracks, audits and figures survive.
if nargin<1,options=struct();end
[baseCfg,sensor,base,cases]=isar_grid_setup();
if isfield(options,'Cases')
    cases=options.Cases;
else
    cases=cases(cases.spinFrequencyHz==1,:);
end
root=fileparts(fileparts(mfilename('fullpath')));
out=char(opt(options,'OutputDirectory',fullfile(root,'results','isar_tip_slowtime_20261020')));
if ~isfolder(out),mkdir(out);end
duration=opt(options,'DurationSec',5);
frames=opt(options,'FrameIndices',1:round(duration/.1));
cfg=struct('fs',1e9,'fs_scope',200e6,'T',1e-3,'Tr',100e-6, ...
    'fL',11.8e9,'fH',17.8e9,'kr',6e13,'propagationSpeed',3e8, ...
    'fc',14.8e9,'lambda',3e8/14.8e9,'Ta',.001,'Na',32, ...
    'startTimeSec',0,'nativeWorkers',opt(options,'NativeWorkers',8), ...
    'generatePriorWaveform',false,'recordPolicy',"TIGHT_EMPTY_TAIL");
p=isarpath_settings(cfg);p.rangeAlignmentMode="FRACTIONAL";
p.phaseCompensationMode="MODEL_CA";
nf=64*cfg.Na;fd=(-nf/2:nf/2-1)'/(cfg.Ta*nf);
keep=abs(fd)<=opt(options,'SavedDopplerHalfSpanHz',90);
rx=(-512:511)*p.rangeBinSpacingM+p.cropCenterOffsetBinsEffective*p.rangeBinSpacingM;
% Search ahead of the PRIOR-PREDICTED processed centre, not coordinate zero:
% the native CV dispersion/crop gauge leaves a velocity-dependent range
% offset even after the CA slow phase correction. Use no tip truth to track.
gateHalfWidthBins=3;gatePolicy="PRIOR_PROCESSING_CENTRE_V2";
diary(fullfile(out,'run_log.txt'));diaryGuard=onCleanup(@()diary('off')); %#ok<NASGU>
timer=tic;tracks=table();audits=table();allTruth=table();sceneResults=cell(height(cases),1);
for ci=1:height(cases)
    q=cases(ci,:);checkpoint=fullfile(out,char(q.caseId+"_tip_result.mat"));
    reuse=false;
    if opt(options,'Resume',true)&&isfile(checkpoint)
        cache=load(checkpoint,'one');
        reuse=isfield(cache.one,'gatePolicy')&&cache.one.gatePolicy==gatePolicy;
    end
    if reuse
        assert(isequal(cache.one.frameIndices,frames),'Checkpoint frame indices differ.');
        assert(isequal(cache.one.dopplerHz,fd(keep)),'Checkpoint Doppler coordinates differ.');
        one=cache.one;fprintf('REUSE %s\n',q.caseId);
    else
        work=baseCfg;work.bodyRotationMode=q.mode;
        work.spinFrequencyHz=q.spinFrequencyHz;work.precessionFrequencyHz=q.precessionFrequencyHz;
        work.precessionHalfAngleDeg=q.precessionHalfAngleDeg;
        work.swingFrequencyHz=q.swingFrequencyHz;
        work.swingAngleMinDeg=-q.swingAmplitudeDeg;
        work.swingAngleMaxDeg=q.swingAmplitudeDeg;
        target=isar_apply_micro_motion(base,work);
        powers=zeros(nnz(keep),numel(frames));tipPowers=powers;
        localTracks=table();localAudits=table();
        denseTime=(0:.001:duration)';
        state=scatterermodel.evaluateTarget(target,sensor,denseTime);
        centreRange=10050+3400*denseTime+50*denseTime.^2;
        centreRate=3400+100*denseTime;
        tipRange=state.rangeM(1,:).'-centreRange;
        tipHz=2/cfg.lambda*(state.rangeRateMps(1,:).'-centreRate);
        separation=min(abs(state.rangeM(2:4,:)-state.rangeM(1,:)),[],'all');
        truthTable=table(repmat(q.caseId,numel(denseTime),1),denseTime, ...
            tipRange,100*(tipRange-mean(tipRange)),tipHz, ...
            'VariableNames',{'caseId','slowTimeSec','tipRelativeCentreRangeM', ...
            'tipRangeDeviationCm','tipPhysicalDopplerHz'});
        for fi=1:numel(frames)
            frame=frames(fi);epoch=(frame-1)*.1;t=epoch+(0:cfg.Na-1)'*cfg.Ta;
            prior=struct('epochSec',epoch,'rangeM',10050+3400*epoch+50*epoch^2, ...
                'velocityMps',3400+100*epoch,'accelerationMps2',100,'beatOffsetHz',20e6);
            prior.loStartSec=2*prior.rangeM/3e8-20e6/cfg.kr;
            cfg.T=ceil((prior.loStartSec+cfg.Tr+16e-6)*cfg.fs/32000)*32000/cfg.fs;
            [H,~,audit]=isarpath_synthesize_profiles(target,sensor,cfg,prior,t,p);
            assert(audit.maximumRawBeatHz<100e6&&audit.coherentSumError<1e-11);
            [~,shifts]=isarpath_align_profiles(H(:,:,1),p);
            H=isarpath_align_profiles(H,p,shifts);
            ref=isarpath_prior_phase(cfg,prior,t,p,shifts,"MODEL_CA");
            H=H.*ref.phaseCorrection;
            predictedCentre=(p.cropCenterOffsetBinsEffective+ ...
                ref.fixedGateOffsetHz/p.frequencyBinSpacingHz)*p.rangeBinSpacingM;
            candidate=find(rx-predictedCentre>=-.70&rx-predictedCentre<=-.36);
            assert(~isempty(candidate));
            [~,peak]=max(sum(abs(H(:,candidate,1)).^2,1));centreBin=candidate(peak);
            gate=centreBin+(-gateHalfWidthBins:gateHalfWidthBins);
            % Oracle channel is an audit AFTER selection, never its input.
            [~,isolatedPeakBin]=max(sum(abs(H(:,:,2)).^2,1));
            mainPeakInGate=any(gate==isolatedPeakBin);
            assert(mainPeakInGate,'Measured tip gate missed the diagnostic tip main peak.');
            % Integrate RD power over the small tip range gate, so a
            % sub-bin shift does not switch the phase of an extracted gate.
            spectrum=fftshift(fft(H(:,gate,:),nf,1),1)/cfg.Na;
            fullPower=sum(abs(spectrum(:,:,1)).^2,2);
            tipPower=sum(abs(spectrum(:,:,2)).^2,2);
            powers(:,fi)=fullPower(keep);tipPowers(:,fi)=tipPower(keep);
            observed=peak_frequency(fullPower,fd);
            isolated=peak_frequency(tipPower,fd);
            receiveMid=mean(t)+2*prior.rangeM/3e8+cfg.Tr/2;
            midState=scatterermodel.evaluateTarget(target,sensor,receiveMid);
            trueHz=2/cfg.lambda*(midState.rangeRateMps(1)-3400-100*receiveMid);
            tail=sum(H(:,gate,3:5),3);tip=H(:,gate,2);
            leakageDb=10*log10(sum(abs(tail).^2,'all')/sum(abs(tip).^2,'all'));
            localTracks=[localTracks;table(q.caseId,frame,epoch,receiveMid,observed,isolated,trueHz, ...
                observed-isolated,observed-trueHz,rx(centreBin),leakageDb,max(abs(shifts)), ...
                predictedCentre,rx(centreBin)-predictedCentre,mainPeakInGate, ...
                'VariableNames',{'caseId','frameIndex','epochSec','receiveMidpointSec', ...
                'extractedTipPeakHz','isolatedTipPeakHz','physicalTipMidpointHz', ...
                'tailInducedPeakDifferenceHz','peakMinusPhysicalHz','tipGateCentreM', ...
                'tailLeakageRelativeTipDb','maximumRangeShiftBins', ...
                'priorProcessedCentreM','tipGateRelativePriorCentreM','tipMainPeakInGate'})]; %#ok<AGROW>
            localAudits=[localAudits;table(q.caseId,frame,audit.maximumRawBeatHz,audit.coherentSumError, ...
                'VariableNames',{'caseId','frameIndex','maximumRawBeatHz','coherentSumError'})]; %#ok<AGROW>
            fprintf('%s | %d/%d | tip %.3f Hz | physical %.3f Hz | tail %.1f dB | %.1f s\n', ...
                q.caseId,fi,numel(frames),observed,trueHz,leakageDb,toc(timer));
        end
        one=struct('scene',q,'frameIndices',frames,'dopplerHz',fd(keep), ...
            'tipGatePower',powers,'isolatedTipGatePower',tipPowers, ...
            'tracks',localTracks,'audits',localAudits,'truth',truthTable, ...
            'minimumTipTailRangeSeparationM',separation, ...
            'gateHalfWidthM',gateHalfWidthBins*p.rangeBinSpacingM, ...
            'rawEchoSaved',false,'hrrpSaved',false,'slowTimeMeanRemoved',false, ...
            'typicalPointCompensationEnabled',false,'rangeGateSelectedFromFullScene',true, ...
            'gatePolicy',gatePolicy);
        save(checkpoint,'one');
    end
    sceneResults{ci}=one;tracks=[tracks;one.tracks];audits=[audits;one.audits]; %#ok<AGROW>
    allTruth=[allTruth;one.truth]; %#ok<AGROW>
    callback=opt(options,'CaseCompletedFcn',[]);
    if ~isempty(callback),callback(one);end
end
writetable(tracks,fullfile(out,'tip_extracted_tracks.csv'));
writetable(audits,fullfile(out,'waveform_audit.csv'));
writetable(allTruth,fullfile(out,'tip_geometry_tracks.csv'));
summary=table();
for ci=1:height(cases)
    one=sceneResults{ci};q=one.scene;X=one.tracks;
    row=q;row.minimumTipTailRangeSeparationM=one.minimumTipTailRangeSeparationM;
    row.maximumTipPhysicalDopplerHz=max(abs(one.truth.tipPhysicalDopplerHz));
    row.tipRangeHalfSpanCm=.5*(max(one.truth.tipRangeDeviationCm)-min(one.truth.tipRangeDeviationCm));
    row.maximumTailLeakageDb=max(X.tailLeakageRelativeTipDb);
    row.maximumTailInducedPeakDifferenceHz=max(abs(X.tailInducedPeakDifferenceHz));
    row.medianPeakMinusPhysicalHz=median(X.peakMinusPhysicalHz);
    row.maximumAbsExtractedPeakHz=max(abs(X.extractedTipPeakHz));
    summary=[summary;row]; %#ok<AGROW>
end
writetable(summary,fullfile(out,'tip_extraction_summary.csv'));
result=struct('configuration',cfg,'processing',p,'cases',cases,'summary',summary, ...
    'tracks',tracks,'outputDirectory',string(out),'rangeGatePolicy',gatePolicy, ...
    'displayNormalizedPerCpi',true,'typicalPointCompensationEnabled',false, ...
    'slowTimeMeanRemoved',false,'durationSec',duration);
save(fullfile(out,'summary.mat'),'result');
if opt(options,'GenerateFigures',true),render_results(sceneResults,out);end
disp(summary);fprintf('DONE %s | %.1f s\n',out,toc(timer));
end

function hz=peak_frequency(power,f)
[~,i]=max(power);offset=0;
if i>1&&i<numel(power)
    y=log(max(power(i+(-1:1)),realmin));den=y(1)-2*y(2)+y(3);
    if den<0,offset=.5*(y(1)-y(3))/den;end
end
hz=f(i)+offset*(f(2)-f(1));
end

function render_results(scenes,out)
names=["tip_slowtime_doppler","tip_physical_doppler","tip_range_track"];
mainTitles=["锥尖距离门慢时间—多普勒图（横向：s；纵向：Hz）", ...
    "锥尖物理多普勒轨迹（横向：s；纵向：Hz）", ...
    "锥尖距离起伏（横向：s；纵向：cm）"];
for kind=1:3
    fig=figure('Visible','off','Color','w','Position',[30 30 2200 1450]);
    layout=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
    for ci=1:numel(scenes)
        one=scenes{ci};q=one.scene;ax=nexttile(layout,ci);hold(ax,'on');
        if kind==1
            db=10*log10(max(one.tipGatePower./max(one.tipGatePower,[],1),1e-4));
            imagesc(ax,one.tracks.receiveMidpointSec,one.dopplerHz,db);
            axis(ax,'xy');colormap(ax,'turbo');clim(ax,[-30 0]);ylim(ax,[-85 85]);
            h1=plot(ax,one.truth.slowTimeSec,one.truth.tipPhysicalDopplerHz,'w--','LineWidth',2);
            h2=plot(ax,one.tracks.receiveMidpointSec,one.tracks.extractedTipPeakHz,'k-','LineWidth',2);
            if ci==1,legend(ax,[h1 h2],{'物理真值','提取峰值'},'Location','southoutside', ...
                'Orientation','horizontal','FontSize',30,'FontName','Microsoft YaHei');end
            bar=colorbar(ax);bar.FontSize=30;bar.Ticks=[-30 -20 -10 0];
        elseif kind==2
            plot(ax,one.truth.slowTimeSec,one.truth.tipPhysicalDopplerHz, ...
                'Color',[.07 .32 .60],'LineWidth',2.2);ylim(ax,[-60 60]);grid(ax,'on');
        else
            plot(ax,one.truth.slowTimeSec,one.truth.tipRangeDeviationCm, ...
                'Color',[.07 .32 .60],'LineWidth',2.2);ylim(ax,[-3.3 3.3]);grid(ax,'on');
        end
        xlim(ax,[0 5]);xticks(ax,0:1:5);
        title(ax,sprintf('f_s = 1 Hz,  f_p = %g Hz,  α = %g°', ...
            q.precessionFrequencyHz,q.precessionHalfAngleDeg), ...
            'FontSize',30,'FontWeight','normal','FontName','Microsoft YaHei');
        ax.FontSize=30;ax.FontName='Microsoft YaHei';ax.LineWidth=1.2;ax.TickDir='out';
        ax.XLabel.String='';ax.YLabel.String='';
    end
    sgtitle(layout,mainTitles(kind),'FontSize',30,'FontName','Microsoft YaHei','FontWeight','normal');
    set(findall(fig,'Type','text'),'FontSize',30);
    savefig(fig,fullfile(out,char(names(kind)+".fig")));
    exportgraphics(fig,fullfile(out,char(names(kind)+".png")),'Resolution',120);
    close(fig);
end
end

function value=opt(s,name,fallback)
if isfield(s,name),value=s.(name);else,value=fallback;end
end
