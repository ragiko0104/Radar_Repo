function results=run_isar_long_micro_motion(options)
%RUN_ISAR_LONG_MICRO_MOTION Stationary / CV-prior cone 5 s RD experiment.
% Default: Fs=200 MHz, CPI=32 ms, hop=100 ms, PRF=[1 4] kHz; spin=1 Hz,
% precession=2 Hz / 10 deg, inertial-Z swing=1 Hz / [-5,+5] deg.
% Modes: spin, coneSpin, swingSpin, coneSwingSpin. No phase autofocus, no
% slow-time mean/MTI, no point-specific correction. Compare noRCMC against
% shared three-pass envelope RCMC. Only figures, CSV trajectories/metrics,
% configuration and small summary tables are saved. NEVER save RF/IF/HRRP.
% Options use the capitalized names in the configuration block below.
if nargin<1,options=struct();end
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'functions')); ensure_md_zt_paths();
out=char(opt(options,'OutputDirectory',fullfile(root,'results','isar_long_micro_motion_20261019')));
if ~isfolder(out),mkdir(out);end
diary(fullfile(out,'run_log.txt'));cleanupDiary=onCleanup(@()diary('off')); %#ok<NASGU>
cfg=struct('sampleRateHz',opt(options,'SampleRateHz',200e6), ...
    'durationSec',opt(options,'DurationSec',5),'cpiSec',opt(options,'CpiSec',.032), ...
    'stepSec',opt(options,'StepSec',.100),'prfHz',opt(options,'PrfHz',[1000 4000]), ...
    'pulseWidthSec',100e-6,'fL',11.8e9,'fH',17.8e9,'propagationSpeed',3e8, ...
    'beatOffsetHz',20e6,'beatGuardHz',90e6,'rangeFftPoints',32768,'rangeBins',1024, ...
    'centerRangeM',10050,'slowPadding',8,'rcmcSearchBins',16, ...
    'spinFrequencyHz',opt(options,'SpinFrequencyHz',1), ...
    'precessionFrequencyHz',opt(options,'PrecessionFrequencyHz',2), ...
    'precessionHalfAngleDeg',opt(options,'PrecessionHalfAngleDeg',10), ...
    'swingFrequencyHz',opt(options,'SwingFrequencyHz',1), ...
    'swingAngleMinDeg',-opt(options,'SwingHalfAmplitudeDeg',5), ...
    'swingAngleMaxDeg',opt(options,'SwingHalfAmplitudeDeg',5), ...
    'batchPulses',opt(options,'BatchPulses',32), ...
    'phaseCompensationEnabled',false,'slowMeanRemovalEnabled',false, ...
    'translationEnabled',false,'rawEchoSaved',false,'hrrfSaved',false);
cfg.initialRadialVelocityMps=opt(options,'InitialRadialVelocityMps',0);
cfg.radialAccelerationMps2=opt(options,'RadialAccelerationMps2',0);
cfg.priorUpdateSec=opt(options,'PriorUpdateSec',.1);
cfg.translationEnabled=cfg.initialRadialVelocityMps~=0||cfg.radialAccelerationMps2~=0;
cfg.accelerationPriorUsed=false;
cfg.dispersionCompensationEnabled=cfg.translationEnabled;
cfg.knownTranslationPhaseRemoved=cfg.translationEnabled;
cfg.exactRangeAndVelocityUpdates=cfg.translationEnabled;
cfg.experimentName=string(opt(options,'ExperimentName',"静止目标5s滑窗ISAR测试"));
validateattributes(cfg.sampleRateHz,{'numeric'},{'scalar','finite','>=',200e6});
validateattributes(cfg.durationSec,{'numeric'},{'scalar','finite','>',cfg.cpiSec});
validateattributes(cfg.batchPulses,{'numeric'},{'scalar','integer','>=',1});
validateattributes(cfg.stepSec,{'numeric'},{'scalar','positive','finite'});
assert(cfg.cpiSec>0&&cfg.cpiSec<=cfg.durationSec);
cfg.prfHz=sort(unique(cfg.prfHz(:).'));
assert(all(cfg.prfHz*cfg.pulseWidthSec<1)&&all(cfg.prfHz>0));
maxPrf=max(cfg.prfHz);
assert(all(abs(maxPrf./cfg.prfHz-round(maxPrf./cfg.prfHz))<1e-10), ...
    'PRFs must divide the maximum PRF so identical pulses can be reused.');
assert(all(abs(cfg.cpiSec*cfg.prfHz-round(cfg.cpiSec*cfg.prfHz))<1e-8));
assert(all(abs(cfg.stepSec*cfg.prfHz-round(cfg.stepSec*cfg.prfHz))<1e-8));
assert(abs(cfg.durationSec*maxPrf-round(cfg.durationSec*maxPrf))<1e-8);
if cfg.translationEnabled
    assert(cfg.cpiSec<=cfg.priorUpdateSec&&abs(cfg.stepSec/cfg.priorUpdateSec-round(cfg.stepSec/cfg.priorUpdateSec))<1e-10, ...
        'Each CPI must start at a prior update and stay inside that update interval.');
end
factors=[cfg.spinFrequencyHz,cfg.precessionFrequencyHz,cfg.swingFrequencyHz];
reuse=opt(options,'ReusePeriodicCycle',true)&&all(abs(factors-round(factors))<1e-12)&&cfg.durationSec>1&&~cfg.translationEnabled;
cfg.periodicCycleReuseEnabled=reuse;
cfg.generatedDurationSec=cfg.durationSec;if reuse,cfg.generatedDurationSec=1;end
cfg.generatedPulses=round(cfg.generatedDurationSec*maxPrf);
cfg.totalRequestedPulsesAtMaxPrf=round(cfg.durationSec*maxPrf);
cfg.fastTimeStateMethod="exact endpoint range/rate cubic Hermite, checked against exact internal states";
cfg.frequencyConvention="Euler spin/precession factors; swing is inertial-Z rigid rotation";
spec=struct('Ntip',1,'radiusM',.3,'phase0Rad',pi/6, ...
    'axialSeparationM',1.5,'alpha0',[1;.3;.3;.3]);
base=scatterermodel.createLegacyDiskTarget(spec,struct('type','accelerated', ...
    'p0_I',[cfg.centerRangeM 0 0],'v0_I',[cfg.initialRadialVelocityMps 0 0], ...
    'a0_I',[cfg.radialAccelerationMps2 0 0],'attitude0YprRad',[-3*pi/4 0 0]));
radar=struct('position_I',[0 0 0],'velocity_I',[0 0 0], ...
    'carrierFrequencyHz',(cfg.fL+cfg.fH)/2,'wavelengthM',cfg.propagationSpeed/((cfg.fL+cfg.fH)/2), ...
    'dopplerSign',-1);
modes=string(opt(options,'Modes',["spin","coneSpin","swingSpin","coneSwingSpin"]));
assert(all(ismember(modes,["spin","coneSpin","swingSpin","coneSwingSpin"])));
generatePlots=opt(options,'GeneratePlots',true);
summary=table();allMetrics=table();verification=struct([]);trajectoryData=struct([]);
for mode=modes(:).'
    timer=tic;motionCfg=cfg;motionCfg.bodyRotationMode=mode;
    target=isar_apply_micro_motion(base,motionCfg);
    trajectory=make_trajectory(target,radar,cfg,mode);
    writetable(trajectory.table,fullfile(out,char(mode+"_trajectories.csv")));
    if reuse
        checks=[0 .137 .731]';
        a=scatterermodel.evaluateTarget(target,radar,checks);
        b=scatterermodel.evaluateTarget(target,radar,checks+1);
        periodicError=max(abs(a.position_I-b.position_I),[],'all');
        assert(periodicError<1e-9,'Cannot reuse a nonperiodic motion cycle.');
    else,periodicError=NaN;end
    fprintf('\n%s: synthesizing %d distinct coherent pulses at %.0f Hz; virtual observation %.3f s\n', ...
        mode,cfg.generatedPulses,maxPrf,cfg.durationSec);
    % Single-precision CROPPED range-profile cache in RAM only. With the
    % default exact 1 s period it is reused to represent all five cycles;
    % lower-PRF pulses are exact subsets of the same high-PRF waveform.
    cache=complex(zeros(cfg.generatedPulses,cfg.rangeBins,4,'single'));
    maxRangeError=0;maxBeat=0;loMin=Inf;loMax=-Inf;maxProcessedBeat=0;
    beforeRangeWidths=[];afterRangeWidths=[];
    for first=1:cfg.batchPulses:cfg.generatedPulses
        rows=first:min(cfg.generatedPulses,first+cfg.batchPulses-1);
        times=(rows'-1)/maxPrf;
        if cfg.translationEnabled
            prior=make_prior(target,radar,times,cfg);
            [cache(rows,:,:),info]=isar_cv_stretch_profiles(target,radar,times,cfg,prior);
            maxProcessedBeat=max(maxProcessedBeat,info.maxProcessedBeatHz);
            beforeRangeWidths=[beforeRangeWidths;info.rangeWidthBeforeDispersionM]; %#ok<AGROW>
            afterRangeWidths=[afterRangeWidths;info.rangeWidthAfterDispersionM]; %#ok<AGROW>
        else
            [cache(rows,:,:),info]=isar_streamed_stretch_profiles(target,radar,times,cfg);
        end
        maxRangeError=max(maxRangeError,info.maxRangeInterpolationErrorM);
        maxBeat=max(maxBeat,info.maxAbsBeatHz);
        loMin=min(loMin,min(info.loDelaySec));loMax=max(loMax,max(info.loDelaySec));
        if first==1||mod(first-1,cfg.batchPulses*25)==0
            fprintf('  generated %d/%d pulses, peak beat %.3f MHz, elapsed %.1f s\n', ...
                rows(end),cfg.generatedPulses,maxBeat/1e6,toc(timer));
        end
    end
    if ~cfg.translationEnabled
        assert(loMax-loMin<1e-15,'Stationary center must yield constant tracking LO delay.');
    end
    verification=[verification;struct('mode',mode,'rangeInterpolationErrorM',maxRangeError, ...
        'phaseErrorBoundRad',4*pi*cfg.fH/cfg.propagationSpeed*maxRangeError, ...
        'maxAbsBeatHz',maxBeat,'loDelaySpanSec',loMax-loMin,'periodicPositionErrorM',periodicError, ...
        'generatedDistinctPulses',cfg.generatedPulses,'runtimeSec',toc(timer), ...
        'maxProcessedBeatHz',maxProcessedBeat, ...
        'medianRangeWidthBeforeDispersionM',median(beforeRangeWidths,1), ...
        'medianRangeWidthAfterDispersionM',median(afterRangeWidths,1))]; %#ok<AGROW>
    caseMetrics=table();
    for prf=cfg.prfHz
        starts=(0:floor((cfg.durationSec-cfg.cpiSec)/cfg.stepSec))*cfg.stepSec;
        nPulse=round(cfg.cpiSec*prf);caseName=sprintf('%s_prf%d',mode,prf);
        gifPath=fullfile(out,[caseName '_rd.gif']);
        if generatePlots
            fig=figure('Visible','off','Color','w','Position',[50 50 1120 510]);
            figureCleanup=onCleanup(@()close_valid(fig));
        end
        stills={};
        for frameIndex=1:numel(starts)
            pulseTime=starts(frameIndex)+(0:nPulse-1)'/prf;
            index=round(pulseTime*maxPrf);
            if reuse,index=mod(index,cfg.generatedPulses);end
            H=cache(index+1,:,:);
            frameCfg=cfg;
            if cfg.translationEnabled
                frameCfg.frameMotionPrior=make_prior(target,radar,starts(frameIndex),cfg);
            end
            frame=isar_long_rd_frame(H,frameCfg,prf,target,radar,pulseTime);
            M=frame.metrics;M.mode=repmat(mode,height(M),1);
            M.prfHz=repmat(prf,height(M),1);M.frameIndex=repmat(frameIndex,height(M),1);
            M.frameStartSec=repmat(starts(frameIndex),height(M),1);
            M.frameCenterSec=repmat(frame.timeSec,height(M),1);
            caseMetrics=[caseMetrics;M]; %#ok<AGROW>
            if generatePlots
                draw_rd_frame(fig,frame,cfg,prf,mode,trajectory);
                rgb=print(fig,'-RGBImage','-r95');
                [indexed,map]=rgb2ind(rgb,256,'nodither');
                if frameIndex==1
                    imwrite(indexed,map,gifPath,'gif','LoopCount',Inf,'DelayTime',cfg.stepSec);
                else
                    imwrite(indexed,map,gifPath,'gif','WriteMode','append','DelayTime',cfg.stepSec);
                end
                if ismember(frameIndex,[1 4 7 10]),stills{end+1}=rgb;end %#ok<AGROW>
            end
        end
        if generatePlots
            if numel(stills)==4
                contact=[stills{1},stills{2};stills{3},stills{4}];
                imwrite(contact,fullfile(out,[caseName '_snapshots.png']));
            else,imwrite(rgb,fullfile(out,[caseName '_snapshots.png']));end
            clear figureCleanup;
        end
        X=caseMetrics(caseMetrics.prfHz==prf,:);
        raw=X(X.method=="noRCMC"&X.pointIndex>1,:);
        aligned=X(X.method=="commonRCMC"&X.pointIndex>1,:);
        maxTruth=max(raw.maxAbsTruthDopplerHz);
        maxTruth=max(maxTruth,max(X.maxAbsTruthDopplerHz));
        entry=struct('mode',mode,'prfHz',prf,'pulseCountPerCpi',nPulse, ...
            'durationSec',cfg.durationSec,'cpiSec',cfg.cpiSec,'stepSec',cfg.stepSec, ...
            'frameCount',numel(starts),'requestedPulseCount',round(prf*cfg.durationSec), ...
            'maxTruthDopplerFcHz',maxTruth,'maxTruthDopplerHighBandHz',maxTruth*cfg.fH/radar.carrierFrequencyHz, ...
            'aliasRiskFc',maxTruth>=prf/2,'aliasRiskHighBand',any(X.frameAliasRiskAtHighBand), ...
            'maxAbsBeatMHz',maxBeat/1e6,'loDelaySpanNs',(loMax-loMin)*1e9, ...
            'maxTailRangePeakSpanRawCm',100*max(raw.residualHrrfPeakSpanM), ...
            'maxTailRangePeakSpanAlignedCm',100*max(aligned.residualHrrfPeakSpanM), ...
            'medianTailWidthRatioRaw',median(raw.widthToIdealRatio,'omitnan'), ...
            'medianTailWidthRatioAligned',median(aligned.widthToIdealRatio,'omitnan'), ...
            'maxTailWidthRatioRaw',max(raw.widthToIdealRatio), ...
            'maxTailWidthRatioAligned',max(aligned.widthToIdealRatio), ...
            'alignedTailFractionRangeSpanOver2p5cm',mean(aligned.residualHrrfPeakSpanM>.025), ...
            'alignedTailFractionOtherPeakAboveMinus3dB',mean(aligned.largestOtherPeakDb>-3), ...
            'rangeInterpolationErrorM',maxRangeError,'phaseAutofocus',false,'slowMeanRemoval',false);
        summary=[summary;struct2table(entry)]; %#ok<AGROW>
        fprintf('  %s: %d frames; alias=%d; tail residual RCM max %.2f cm; median width/ideal %.2f -> %.2f\n', ...
            caseName,numel(starts),entry.aliasRiskHighBand,entry.maxTailRangePeakSpanAlignedCm, ...
            entry.medianTailWidthRatioRaw,entry.medianTailWidthRatioAligned);
    end
    allMetrics=[allMetrics;caseMetrics]; %#ok<AGROW>
    if generatePlots,plot_trajectories(trajectory,caseMetrics,cfg,mode,out);end
    trajectoryData=[trajectoryData;struct('mode',mode,'maxTipRangeSpanM', ...
        max(trajectory.rangeM(:,1))-min(trajectory.rangeM(:,1)), ...
        'maxTipDopplerHz',max(abs(trajectory.fdHz(:,1))))]; %#ok<AGROW>
    clear cache;
    writetable(summary,fullfile(out,'summary.csv'));
    writetable(allMetrics,fullfile(out,'point_frame_metrics.csv'));
end
results=struct('experimentName',cfg.experimentName,'configuration',cfg,'summary',summary,'pointFrameMetrics',allMetrics, ...
    'verification',verification,'motionChecks',trajectoryData,'outputDirectory',string(out), ...
    'storagePolicy',"No raw RF, IF, HRRP, profile cache or RD cube is saved.");
save(fullfile(out,'summary.mat'),'results','-v7');
write_report(results,out);
disp(summary(:,{'mode','prfHz','aliasRiskHighBand','maxTailRangePeakSpanAlignedCm', ...
    'medianTailWidthRatioRaw','medianTailWidthRatioAligned'}));
fprintf('Done. Small results and figures only: %s\n',out);
end

function trajectory=make_trajectory(target,radar,cfg,mode)
t=(0:round(cfg.durationSec*1000))'/1000;
s=scatterermodel.evaluateTarget(target,radar,t);
centre=target.centerStateFcn(target,radar,t);
Rc=vecnorm(centre.position_I,2,2);vc=centre.velocity_I(:,1);
rho=(s.rangeM-Rc.').';fd=s.dopplerHz.'+2*vc/radar.wavelengthM;
position=s.position_I;position(:,1,:)=position(:,1,:)-reshape(Rc,1,1,[]);
prior=make_prior(target,radar,t,cfg);
modelR=prior.rangeM+prior.velocityMps.*(t-prior.epochSec);
residualFd=s.dopplerHz.'+2*prior.velocityMps/radar.wavelengthM;
T=table();
for j=1:4
    xyz=reshape(position(j,:,:),3,[]).';
    row=table(repmat(mode,numel(t),1),repmat(j,numel(t),1),t,xyz(:,1),xyz(:,2),xyz(:,3), ...
        rho(:,j),fd(:,j),mod(fd(:,j)+500,1000)-500,mod(fd(:,j)+2000,4000)-2000, ...
        'VariableNames',{'mode','pointIndex','timeSec','relativeXLosM','relativeYM','relativeZM', ...
        'rangeOffsetM','trueDopplerFcHz','foldedDoppler1kHz','foldedDoppler4kHz'});
    if cfg.translationEnabled
        row.trueAbsoluteDopplerFcHz=s.dopplerHz(j,:).';
        row.centerRangeM=Rc;row.centerVelocityMps=vc;
        row.priorEpochSec=prior.epochSec;row.priorRangeM=prior.rangeM;row.priorVelocityMps=prior.velocityMps;
        row.cvModelRangeM=modelR;row.residualRangeAfterCvM=s.rangeM(j,:).'-modelR;
        row.residualDopplerAfterCvHz=residualFd(:,j);
        row.foldedResidualDoppler1kHz=mod(residualFd(:,j)+500,1000)-500;
        row.foldedResidualDoppler4kHz=mod(residualFd(:,j)+2000,4000)-2000;
        row.withinImagingWindow=mod(t+1e-10,cfg.stepSec)<cfg.cpiSec;
    end
    T=[T;row]; %#ok<AGROW>
end
trajectory=struct('timeSec',t,'rangeM',rho,'fdHz',fd,'position',position,'table',T);
if cfg.translationEnabled
    visible=mod(t+1e-10,cfg.stepSec)<cfg.cpiSec;
    residualFd(~visible,:)=NaN;trajectory.fdHz=residualFd;
end
end

function draw_rd_frame(fig,frame,cfg,prf,mode,traj)
clf(fig);tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
labels={'No range alignment','Shared envelope RCMC (3 passes)'};
methodNames=["noRCMC","commonRCMC"];
amplitudes={frame.rawAmplitude,frame.alignedAmplitude};colors=lines(4);
limit=min(prf/2,100*(ceil(max(abs(traj.fdHz),[],'all','omitnan')/100)+1));
for method=1:2
    nexttile;imagesc(frame.rangeAxisM,frame.dopplerAxisHz, ...
        20*log10(max(amplitudes{method},1e-3)));axis xy;
    xlim([-1 1]);ylim([-limit limit]);clim([-40 0]);hold on;
    m=frame.metrics(frame.metrics.method==methodNames(method),:);
    for j=1:4
        scatter(m.peakRangeOffsetM(j),m.peakDopplerHz(j),32,colors(j,:),'o','LineWidth',1.2);
        text(m.peakRangeOffsetM(j)+.025,m.peakDopplerHz(j),sprintf('P%d',j), ...
            'Color','w','FontWeight','bold','FontSize',9);
    end
    title(labels{method});xlabel('Range relative to motion model (m)');ylabel('Residual Doppler (Hz)');
    colorbar;
end
colormap(fig,'turbo');
sgtitle(sprintf('%s | PRF %.0f Hz | t = %.3f s | CPI %.0f ms | no autofocus / no slow-time DC removal', ...
    mode,prf,frame.timeSec,cfg.cpiSec*1000),'Interpreter','none','FontSize',11);
end

function plot_trajectories(traj,M,cfg,mode,out)
colors=lines(4);names={'P1 tip','P2 base','P3 base','P4 base'};
fig=figure('Visible','off','Color','w','Position',[60 60 1200 850]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile;hold on;sel=traj.timeSec<=min(1,cfg.durationSec);
for j=1:4
    xyz=reshape(traj.position(j,:,:),3,[]).';
    plot3(xyz(sel,1),xyz(sel,2),xyz(sel,3),'Color',colors(j,:),'LineWidth',1.2);
    scatter3(xyz(1,1),xyz(1,2),xyz(1,3),36,colors(j,:),'filled');
end
quiver3(-.95,0,0,.65,0,0,0,'r','LineWidth',1.6);text(-.95,0,.1,'LOS +X','Color','r');
axis equal;grid on;view(35,25);xlabel('LOS X (m)');ylabel('Y (m)');zlabel('Z (m)');
title('Geometric point trajectories (first cycle)');
nexttile;plot(traj.timeSec,traj.rangeM,'LineWidth',1.2);grid on;
xlabel('Time (s)');ylabel('Range relative to center (m)');title('Exact geometric range');
legend(names,'Location','best');
nexttile;plot(traj.timeSec,traj.fdHz,'LineWidth',1.1);hold on;yline(500,'k:');yline(-500,'k:');grid on;
xlabel('Time (s)');ylabel('True Doppler at 14.8 GHz (Hz)');title('Unaliased geometric Doppler; dotted lines = 1 kHz Nyquist');
if cfg.translationEnabled
    ylabel('Residual Doppler after CV (Hz)');title('True CV-residual Doppler in 32 ms windows; dotted: 1 kHz Nyquist');
end
nexttile;hold on;
for j=1:4
    x=M(M.method=="noRCMC"&M.pointIndex==j,:);
    for prf=cfg.prfHz
        y=x(x.prfHz==prf,:);marker='o';if prf==min(cfg.prfHz),marker='+';end
        plot(y.frameCenterSec,y.peakDopplerHz,marker,'Color',colors(j,:),'MarkerSize',4);
    end
end
grid on;xlabel('Window center (s)');ylabel('Measured point peak Doppler (Hz)');
title('Isolated-component diagnostic peaks: + lower PRF, o higher PRF');
fp=0;alpha=0;fn=0;amplitude=0;
if contains(lower(mode),"cone"),fp=cfg.precessionFrequencyHz;alpha=cfg.precessionHalfAngleDeg;end
if contains(lower(mode),"swing"),fn=cfg.swingFrequencyHz;amplitude=cfg.swingAngleMaxDeg;end
stateLabel='stationary center';
if cfg.translationEnabled,stateLabel=sprintf('v0 %.0f m/s, a %.0f m/s^2; CV prior / %.0f ms', ...
        cfg.initialRadialVelocityMps,cfg.radialAccelerationMps2,1000*cfg.priorUpdateSec);end
sgtitle(sprintf('%s | %s | spin %.1f Hz, coning %.1f Hz / %.0f deg, Z-swing %.1f Hz / +/-%.0f deg', ...
    mode,stateLabel,cfg.spinFrequencyHz,fp,alpha,fn,amplitude),'Interpreter','none');
exportgraphics(fig,fullfile(out,char(mode+"_trajectories.png")),'Resolution',150);close(fig);
end

function write_report(results,out)
cfg=results.configuration;S=results.summary;
path=fullfile(out,'REPORT.txt');fid=fopen(path,'w','n','UTF-8');guard=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n\n',cfg.experimentName);
if cfg.translationEnabled
    write_cv_report(fid,results);
    return;
end
fprintf(fid,'Fs %.0f MHz; duration %.3f s; CPI %.0f ms; step %.0f ms; no autofocus, no slow-time mean removal.\n\n', ...
    cfg.sampleRateHz/1e6,cfg.durationSec,cfg.cpiSec*1000,cfg.stepSec*1000);
fprintf(fid,'Spin factor %.3g Hz; precession %.3g Hz / %.3g deg; inertial-Z swing %.3g Hz, [%.3g, %.3g] deg. Swing starts at its minimum.\n\n', ...
    cfg.spinFrequencyHz,cfg.precessionFrequencyHz,cfg.precessionHalfAngleDeg, ...
    cfg.swingFrequencyHz,cfg.swingAngleMinDeg,cfg.swingAngleMaxDeg);
fprintf(fid,'Center fixed at [10050,0,0] m. Reference body-Y / LOS angle 45 deg; swing changes the instantaneous axis. Body midpoint is the rotation pivot.\n\n');
fprintf(fid,'The reference LO is evaluated every pulse from the center range, so its delay is correctly constant in this stationary case. Continuous delay is used; no 5 ns quantization. Beat convention: echo*conj(LO), -20 MHz, known offset removed before FFT.\n\n');
fprintf(fid,'Analytic low-IF difference phase is sampled at 200 MHz; no real-RF undersampling/mixing. Exact range/rate endpoints and cubic Hermite interpolation retain intrapulse motion. Exact internal checkpoints and an independent direct full-fast-time RF-product unit test validate the approximation.\n\n');
fprintf(fid,'Periodic cache reuse: %d. Distinct generated interval %.3f s. Integer default frequencies make an exact 1 s repeat; an in-RAM cropped-HRRP cycle represents the complete 5 s pulse train, and 1 kHz pulses are exact subsets of the 4 kHz samples. No pulse phase reset occurs at sliding windows. Disable ReusePeriodicCycle for a direct full-duration generation.\n\n', ...
    cfg.periodicCycleReuseEnabled,cfg.generatedDurationSec);
fprintf(fid,'| Mode | PRF Hz | Alias risk (upper band) | Max residual tail RCM cm | Median tail width/ideal: raw | aligned |\n|---|---:|---|---:|---:|---:|\n');
for i=1:height(S)
    fprintf(fid,'| %s | %.0f | %d | %.3f | %.3f | %.3f |\n',S.mode(i),S.prfHz(i), ...
        S.aliasRiskHighBand(i),S.maxTailRangePeakSpanAlignedCm(i),S.medianTailWidthRatioRaw(i),S.medianTailWidthRatioAligned(i));
end
fprintf(fid,'\nRCMC is fitted to the SUM scene and its common shifts are replayed to each point; isolated-point signals are used only for diagnostic widths/labels, never individual focusing. Residual migration over the 2.5 cm nominal range resolution, broadening, split peaks and alias flags must be read together. A narrow local -3 dB width alone does not establish focusing or 5 cm pair separation.\n\n');
fprintf(fid,'The RCMC comparison deliberately retains the legacy complex-profile integer shifts. Without autofocus these shifts can add common phase steps and change the apparent Doppler; no hidden point-reference phase reset is applied. Use the noRCMC panel and geometric tracks as the primary micro-Doppler reference. Residual rotation migration, nonlinear point phase and alignment-induced phase steps are distinct effects.\n\n');
fprintf(fid,'GIF panels: no RCMC / common three-pass envelope RCMC. Shared fixed dB scale relative to unit field amplitude, no frame-by-frame peak renormalization. Labels are peaks of isolated diagnostic components and may merge in the mixed scene.\n\n');
fprintf(fid,'The 100 ms hop means 10 animation frames/s and can give stroboscopic ambiguity for the threefold symmetric base. Dense 1 ms geometric trajectories are supplied; use them rather than animation repetition alone to infer rotation frequencies.\n\n');
fprintf(fid,'No RF, IF, HRRP cache or full range-Doppler cube is saved. CSV files contain trajectories and per-point/per-window measurements; summary.mat contains configuration, small tables and verification scalars only.\n');
end

function close_valid(fig)
if isgraphics(fig),close(fig);end
end

function value=opt(options,name,default)
if isfield(options,name),value=options.(name);else,value=default;end
end

function prior=make_prior(target,radar,t,cfg)
epoch=floor((t(:)+1e-10)/cfg.priorUpdateSec)*cfg.priorUpdateSec;
state=target.centerStateFcn(target,radar,epoch);
prior=struct('epochSec',epoch,'rangeM',vecnorm(state.position_I,2,2), ...
    'velocityMps',sum(state.position_I.*state.velocity_I,2)./vecnorm(state.position_I,2,2));
end

function write_cv_report(fid,results)
cfg=results.configuration;S=results.summary;
fprintf(fid,'Initial radial speed %.0f m/s, acceleration %.0f m/s^2, receding (+LOS). Project units: Mach=340 m/s, g=10 m/s^2.\n',cfg.initialRadialVelocityMps,cfg.radialAccelerationMps2);
fprintf(fid,'Fs %.0f MHz, duration %.3f s, CPI %.0f ms, hop %.0f ms; PRF 1 / 4 kHz. All %d high-PRF pulses are synthesized; no periodic reuse.\n',cfg.sampleRateHz/1e6,cfg.durationSec,cfg.cpiSec*1000,cfg.stepSec*1000,cfg.generatedPulses);
fprintf(fid,'Spin factor 1 Hz, precession 2 Hz / 10 deg, world-Z swing 1 Hz / +/-5 deg. Same four-point cone and initial attitude as the stationary baseline. No attenuation/noise changes are introduced.\n\n');
fprintf(fid,'Exact CENTER range and velocity are supplied at each %.0f ms update. Only their constant-velocity extrapolation is used between updates. The CV prior has no acceleration field; processing never differentiates successive velocity updates and uses no future update. Truth acceleration belongs to the echo generator and diagnostic geometry only.\n',1000*cfg.priorUpdateSec);
fprintf(fid,'Pulse LO delay solves tauLO=2*R_CV(t+tauLO)/c-bias/k. All pulses retain absolute slow-time phase. Conjugation of the analytic CV beat phase corrects intrapulse dispersion / range-Doppler coupling and known translation; known -20 MHz IF is digitally translated exactly, then the central 1024 of 32768 FFT bins are cropped. This is a model-centred frequency window without fractional-bin rounding.\n');
fprintf(fid,'The full CV phase correction includes the terms beyond the legacy first-order dispersion formula; its sign follows echo*conj(LO). It removes only the supplied CV model, with no data-derived phase autofocus.\n');
fprintf(fid,'No typical-point phase compensation, slow-time mean, acceleration search, point-specific range correction or point-specific phase fitting is applied. Shared three-pass envelope RCMC is an optional comparison, fitted only to the summed scene. Its legacy complex integer shifts can add common phase steps.\n\n');
fprintf(fid,'For a=%.0f m/s^2, the ideal CV-model error at 32 ms is %.3f cm, residual speed %.3f m/s, residual bulk Doppler %.3f Hz at 14.8 GHz (upper band %.3f Hz). At the actual 1/4 kHz last pulse the elapsed times are 31 / 31.75 ms, plus the receive-time delay.\n',cfg.radialAccelerationMps2,100*.5*cfg.radialAccelerationMps2*cfg.cpiSec^2,cfg.radialAccelerationMps2*cfg.cpiSec,2*14.8e9/3e8*cfg.radialAccelerationMps2*cfg.cpiSec,2*cfg.fH/3e8*cfg.radialAccelerationMps2*cfg.cpiSec);
fprintf(fid,'Ideal rectangular-aperture Doppler FWHM is %.5f Hz. It is a point-response bandwidth reference, not proof of 5 cm spatial resolution.\n\n',.885893/cfg.cpiSec);
fprintf(fid,'Mode, PRF, upper-band residual alias, max raw/aligned tail RCM cm, median raw/aligned tail width/ideal\n');
for i=1:height(S)
    fprintf(fid,'%s, %.0f, %d, %.4f / %.4f, %.4f / %.4f\n',S.mode(i),S.prfHz(i),S.aliasRiskHighBand(i),S.maxTailRangePeakSpanRawCm(i),S.maxTailRangePeakSpanAlignedCm(i),S.medianTailWidthRatioRaw(i),S.medianTailWidthRatioAligned(i));
end
fprintf(fid,'\nPer-point/per-window CSVs retain secondary-peak, amplitude-variation, nonlinear-phase and alias diagnostics; a narrow local FWHM alone must not be read as successful focusing. Trajectory CSVs distinguish absolute Doppler, micro-motion Doppler and CV-residual Doppler. Markers in images use isolated responses for diagnostics only, never for correction.\n');
for i=1:numel(results.verification)
    v=results.verification(i);
    fprintf(fid,'%s: peak ADC beat %.6f MHz, post-CV baseband %.6f MHz; median intrapulse range FWHM before [%.4f %.4f %.4f %.4f] cm, after [%.4f %.4f %.4f %.4f] cm; interpolation error %.3g m.\n',v.mode,v.maxAbsBeatHz/1e6,v.maxProcessedBeatHz/1e6,100*v.medianRangeWidthBeforeDispersionM,100*v.medianRangeWidthAfterDispersionM,v.rangeInterpolationErrorM);
end
fprintf(fid,'\nGIFs use a common fixed amplitude scale; 50 frames / 5 s. 100 ms sampling can create stroboscopic motion. No RF, IF, HRRP or RD cube is saved; only small tables, trajectory CSVs, plots, GIFs and configuration are retained.\n');
end
