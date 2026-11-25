function [cfg,sensor,base,cases]=isar_grid_setup(options)
%ISAR_GRID_SETUP Reproducible full-cone frequency/angle study, AF disabled.
if nargin<1,options=struct();end
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'functions')); ensure_md_zt_paths();
% Explicit numerical inputs match the preceding genuine-lib experiment;
% this runner does not require any previously generated result MAT file.
cfg=struct('sampleRateHz',200e6,'fL',11.8e9,'fH',17.8e9, ...
    'pulseWidthSec',100e-6,'propagationSpeed',3e8,'beatOffsetHz',20e6, ...
    'beatGuardHz',90e6,'rangeFftPoints',32768,'rangeBins',1024, ...
    'slowPadding',8,'rcmcSearchBins',16,'cpiSec',.032,'stepSec',.1, ...
    'durationSec',5,'startTimeSec',0,'centerRangeM',10050, ...
    'initialRadialVelocityMps',3400,'radialAccelerationMps2',100);
cfg.waveformBackend="lib";cfg.generationAndCompensationSeparated=true;
cfg.explicitConjugateMixing=true;cfg.fixedLoWithinCpi=true;
cfg.loUpdateUsesAccurateRangeOnly=true;cfg.libCacheExactStates=true;
cfg.modelSlowQuadraticPhaseRemoved=true;cfg.residualVideoPhaseRemoved=true;
cfg.fastMidPhasePreserved=false;cfg.interWindowPulseGenerationSkipped=true;
cfg.periodicCycleReuseEnabled=false;
cfg.libraryEchoGeneratorPath=string(which('radar.echo_generate_target_lfm'));
cfg.experimentName="锥体频率与小进动角参数扫描：固定本振、CV/CA、无特显点";
cfg.enableTypicalPointCompensation=false;cfg.phaseCompensationEnabled=false;
cfg.removeSlowTimeMean=false;cfg.preserveOriginalSlowPhase=false;
cfg.computeOnlyAlignedMetrics=true;cfg.batchPulses=8;
cfg.prfHz=[1000 4000];cfg.dopplerGuardFraction=.9;
cfg.rawEchoSaved=false;cfg.hrrfSaved=false;cfg.shiftPhaseRestorationEnabled=false;
cfg.frameCount=50;cfg.auditSlowStepSec=.000125;
sensor=struct('position_I',[0 0 0],'velocity_I',[0 0 0], ...
    'carrierFrequencyHz',14.8e9,'wavelengthM',cfg.propagationSpeed/14.8e9,'dopplerSign',-1);
motion=struct('type','accelerated','p0_I',[cfg.centerRangeM 0 0], ...
    'v0_I',[cfg.initialRadialVelocityMps 0 0],'a0_I',[cfg.radialAccelerationMps2 0 0], ...
    'attitude0YprRad',[-3*pi/4 0 0]);
spec=struct('Ntip',1,'radiusM',.3,'phase0Rad',pi/6,'axialSeparationM',1.5,'alpha0',[1;.3;.3;.3]);
base=scatterermodel.createLegacyDiskTarget(spec,motion);
family=opt(options,'Family',1);pairOnly=opt(options,'PairFrequencies',false);
cases=table();
for group=family(:).'
    if group==1
        rates=[1 3];angles=[1 3];swingHz=0;swingAmplitude=0;mode="coneSpin";
    else
        assert(group==2);
        rates=[2 3.5];angles=[5 7];mode="coneSwingSpin";
        if isfield(options,'SwingFrequencyHz')&&isfield(options,'SwingAmplitudeDeg')
            swingHz=options.SwingFrequencyHz;swingAmplitude=options.SwingAmplitudeDeg;
        else
            assert(isfield(options,'AllSwingLevels')&&options.AllSwingLevels, ...
                'Select swing amplitude/frequency, or explicitly request AllSwingLevels.');
            swingHz=[.5 1 1];swingAmplitude=[1 3 5];
        end
    end
    for swingIndex=1:numel(swingHz)
    for fs=rates
        for fp=rates
            if pairOnly&&fs~=fp,continue;end
            for alpha=angles
                id=string(sprintf('G%d_s%g_p%g_a%g',group,fs,fp,alpha));
                if group==2,id=id+string(sprintf('_w%g_f%g',swingAmplitude(swingIndex),swingHz(swingIndex)));end
                id=replace(id,'.','p');
                row=table(id,group,mode,fs,fp,alpha,swingAmplitude(swingIndex),swingHz(swingIndex), ...
                    'VariableNames',{'caseId','family','mode','spinFrequencyHz','precessionFrequencyHz', ...
                    'precessionHalfAngleDeg','swingAmplitudeDeg','swingFrequencyHz'});
                cases=[cases;row]; %#ok<AGROW>
            end
        end
    end
    end
end
end

function value=opt(s,name,fallback)
if isfield(s,name),value=s.(name);else,value=fallback;end
end
