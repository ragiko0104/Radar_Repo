function [H,refs,audit]=isarpath_synthesize_profiles(target,sensor,cfg,prior,t,p)
%ISARPATH_SYNTHESIZE_PROFILES Stream the actual original ISARpath frontend.
% The target RF is synthesized by radar.echo_generate_target_lfm at 1 GHz.
% Separate CV/CA centre-only reference RF uses prior parameters exclusively.
% Every row passes unchanged radar.lfm_local_reference_beat (fold, filtering,
% positive branch, resample, fast-time DC). Only small cropped HRRPs survive.
% Four diagnostic point channels are summed AFTER this real-linear frontend.
% A direct mixed-scene row checks this equivalence on the first pulse.
if isfield(cfg,'nativeWorkers')&&cfg.nativeWorkers>1&&numel(t)>1
    pool=gcp('nocreate');if isempty(pool),pool=parpool('Processes',cfg.nativeWorkers);end
    assert(~isa(pool,'parallel.ThreadPool'),'The library resample MEX requires process workers.');
    serial=cfg;serial.nativeWorkers=1;serial.workerThreadLimit=1;
    futures=parallel.FevalFuture.empty(0,numel(t));
    for j=1:numel(t)
        one=prior;one.verifySceneSum=j==1;
        futures(j)=parfeval(pool,@isarpath_synthesize_profiles,3,target,sensor,serial,one,t(j),p);
    end
    H=complex(zeros(numel(t),p.hrrfBins,5));refs=complex(zeros(numel(t),p.hrrfBins,2));
    audit=[];
    for j=1:numel(t)
        [index,h,r,a]=fetchNext(futures);H(index,:,:)=h;refs(index,:,:)=r;
        if isempty(audit),audit=a;
        else
            audit.coherentSumError=max(audit.coherentSumError,a.coherentSumError);
            audit.maximumRawBeatHz=max(audit.maximumRawBeatHz,a.maximumRawBeatHz);
        end
    end
    audit.generatedPulseCount=numel(t);return
end
if isfield(cfg,'workerThreadLimit'),maxNumCompThreads(cfg.workerThreadLimit);end
tf=(0:round(cfg.fs*cfg.T)-1)/cfg.fs;np=numel(t);hasTarget=~isempty(target);
verifySceneSum=~isfield(prior,'verifySceneSum')||prior.verifySceneSum;
nModel=2;
if isfield(cfg,'generatePriorWaveform')&&~cfg.generatePriorWaveform,nModel=0;end
H=complex(zeros(np,p.hrrfBins,5));refs=complex(zeros(np,p.hrrfBins,2));
audit=struct('coherentSumError',0,'referenceQuantizationSec',NaN,'generatedPulseCount',np, ...
    'rfSampleRateHz',cfg.fs,'ifSampleRateHz',cfg.fs_scope,'recordLengthSec',cfg.T, ...
    'motionCompensatedDuringGeneration',false,'maximumRawBeatHz',0);
if hasTarget
    provider=target;provider.groupId=(1:4)';provider.groupKeys=["P1";"P2";"P3";"P4"];provider.groupNames=provider.groupKeys;
    % Auto-gating asks only for centre position. Computing full 3-D attitude
    % for a million empty-tail samples would change cost, not the waveform.
    provider.centerStateFcn=@(value,~,query) struct('position_I',value.payload.motion.positionFcn(query));
end
cached=[];cacheFirstTime=0;
request=struct('pulseTimeSec',0,'fastTimeSec',tf,'f0Hz',cfg.fL,'f1Hz',cfg.fH, ...
    'fLoHz',0,'pulseWidthSec',cfg.Tr,'deltaTauSec',0,'propagationSpeed',cfg.propagationSpeed, ...
    'waveformType','LFM','pathLossMode','none','targetAmplitude',1);
options=struct('Mode','fastTimeExact','ReturnGroups',true,'ReturnPerTarget',false, ...
    'AutoFastTimeGate',true,'FastTimeChunkSize',32768,'ScattererChunkSize',4);
beatRequest=struct('fastTimeSec',tf,'referenceStartTimeSec',prior.loStartSec, ...
    'f0Hz',cfg.fL,'f1Hz',cfg.fH,'pulseWidthSec',cfg.Tr,'referenceDelaySec',0,'outputSampleRateHz',cfg.fs_scope);
beatOptions=struct('Mode','ANALYTIC_IF','RemoveMirrorFrequency',true, ...
    'RemoveSamplingAlias',true,'RemoveTimeOverlapCrosstalk',true,'RemoveDc',true, ...
    'LowpassCenterHz',cfg.fs_scope/4,'LowpassBandwidthHz',cfg.fs_scope/2, ...
    'FilterSuppressionDb',120,'ReturnPerTarget',false);
for j=1:np
    nch=nModel+4*hasTarget+double(hasTarget&&j==1&&verifySceneSum);
    rf=complex(zeros(nch,numel(tf)));first=1;
    if hasTarget
        useCache=~isfield(cfg,'cacheExactGeometry')||cfg.cacheExactGeometry;
        if useCache
            centre=target.payload.motion.positionFcn(t(j)+[0;cfg.T]);
            distances=vecnorm(centre-sensor.position_I,2,2);
            firstSample=max(0,floor(2*(min(distances)-target.boundingRadiusM)/cfg.propagationSpeed*cfg.fs)-2);
            lastSample=min(numel(tf)-1,ceil((2*(max(distances)+target.boundingRadiusM)/cfg.propagationSpeed+cfg.Tr)*cfg.fs)+2);
            cacheTimes=t(j)+tf(firstSample+1:lastSample+1);cacheFirstTime=cacheTimes(1);
            cached=scatterermodel.evaluateTarget(target,sensor,cacheTimes);
            queryFast=tf(firstSample+1:lastSample+1);
            lo=round(prior.loStartSec*cfg.fs)/cfg.fs;
            delay=2*cached.rangeM/cfg.propagationSpeed;
            echoTime=queryFast-delay;
            rawFrequency=cfg.kr*(delay-lo)+2*cached.rangeRateMps/cfg.propagationSpeed.*(cfg.fL+cfg.kr*echoTime);
            active=echoTime>=0&echoTime<cfg.Tr&queryFast>=lo&queryFast<lo+cfg.Tr;
            audit.maximumRawBeatHz=max(audit.maximumRawBeatHz,max(abs(rawFrequency(active))));
            provider.evaluateFcn=@cached_state;
        else
            provider.evaluateFcn=target.evaluateFcn;
        end
        request.pulseTimeSec=t(j);
        [echo,generation]=radar.echo_generate_target_lfm(provider,sensor,request,options);
        assert(strcmp(generation.mode,'fastTimeExact'));
        for point=1:4,rf(point,:)=echo.groups(point).signal;end
        if j==1&&verifySceneSum,rf(end,:)=echo.total;end
        clear echo
        first=5;
    end
    % Prior RF generation is separate from all compensation. Model time is
    % absolute received time, matching fastTimeExact; no truth acceleration
    % enters the CV reference. The reference has one point at the body origin.
    for model=1:nModel
        acc=0;if model==2,acc=prior.accelerationMps2;end
        dt=t(j)-prior.epochSec+tf;
        R=prior.rangeM+prior.velocityMps*dt+.5*acc*dt.^2;
        u=tf-2*R/cfg.propagationSpeed;
        valid=u>=0 & u<cfg.Tr;
        rf(first+model-1,valid)=exp(1j*2*pi*(cfg.fL*u(valid)+.5*cfg.kr*u(valid).^2));
        vr=prior.velocityMps+acc*dt;
        lo=round(prior.loStartSec*cfg.fs)/cfg.fs;
        f=cfg.kr*(2*R/cfg.propagationSpeed-lo)+2*vr/cfg.propagationSpeed.*(cfg.fL+cfg.kr*u);
        overlap=valid & tf>=lo & tf<lo+cfg.Tr;
        audit.maximumRawBeatHz=max(audit.maximumRawBeatHz,max(abs(f(overlap))));
    end
    [beat,bi]=radar.lfm_local_reference_beat({rf},beatRequest,beatOptions);clear rf
    audit.referenceQuantizationSec=bi.referenceQuantizationErrorSec;
    % Each diagnostic row is an independent channel with the SAME physical
    % target's LO. It must not be passed as a separate target (1/Ntarget).
    if hasTarget&&j==1&&verifySceneSum
        summed=sum(beat.total(1:4,:),1);direct=beat.total(end,:);
        audit.coherentSumError=norm(summed-direct)/max(norm(direct),eps);
        assert(audit.coherentSumError<1e-11);
    end
    rows=repmat(t(j),nch,1);
    segment=isarpath_extract_segment(beat.total,cfg,prior,rows,p);clear beat
    profiles=isarpath_range_profiles(segment,cfg,prior,rows,p);clear segment
    if hasTarget
        H(j,:,2:5)=permute(profiles(1:4,:),[3 2 1]);H(j,:,1)=sum(H(j,:,2:5),3);
    end
    if nModel>0,refs(j,:,:)=permute(profiles(first:first+1,:),[3 2 1]);end
end

    function state=cached_state(~,~,query,indices)
        column=round((query(:).'-cacheFirstTime)*cfg.fs)+1;
        assert(all(column>=1&column<=size(cached.rangeM,2)));
        assert(max(abs(cacheFirstTime+(column-1)/cfg.fs-query(:).'))<2e-14);
        state=struct('rangeM',cached.rangeM(indices,column),'alpha',cached.alpha(indices,column), ...
            'position_I',cached.position_I(indices,:,column),'velocity_I',cached.velocity_I(indices,:,column), ...
            'rangeRateMps',cached.rangeRateMps(indices,column),'dopplerHz',cached.dopplerHz(indices,column));
    end
end
