function frame=isar_long_rd_frame(pointHrrf,cfg,prf,target,radar,pulseTime,sceneHrrf)
%ISAR_LONG_RD_FRAME Sliding range-Doppler, with/without shared envelope RCMC.
% Optional legacy variance-point autofocus; no slow-time DC/MTI removal.
H=double(pointHrrf);np=size(H,1);nr=size(H,2);ns=size(H,3);
scene=sum(H,3);externalScene=nargin>=7;
if externalScene
    assert(isequal(size(sceneHrrf),[np nr]));scene=double(sceneHrrf);
end
sceneError=norm(sum(H,3)-scene,'fro')/max(norm(scene,'fro'),eps);
assert(sceneError<1e-12,'Diagnostic points do not sum to received scene.');
work=scene;shifts=zeros(np,1);
% External-test style: three iterations, first-pulse then mean envelope.
% Keep the current ISARpath +/-16-bin search, and replay the same shifts
% onto every point. This never uses the isolated points to fit alignment.
for iteration=1:3
    magnitude=abs(work);
    if iteration==1,reference=magnitude(1,:);else,reference=mean(magnitude,1);end
    spectrum=fft([magnitude zeros(np,nr)],2*nr,2);
    refSpectrum=fft([reference zeros(1,nr)],2*nr);
    correlation=real(ifft(spectrum.*conj(refSpectrum),2*nr,2));
    candidates=-cfg.rcmcSearchBins:cfg.rcmcSearchBins;
    [~,best]=max(correlation(:,mod(candidates,2*nr)+1),[],2);
    step=reshape(candidates(best),[],1);
    work=shift_rows(work,step);shifts=shifts+step;
end
aligned=complex(zeros(size(H)));
for j=1:ns,aligned(:,:,j)=shift_rows(H(:,:,j),shifts);end
replayError=norm(sum(aligned,3)-work,'fro')/max(norm(work,'fro'),eps);
assert(replayError<1e-12,'Common RCMC point replay is not coherent.');
doAutofocus=isfield(cfg,'enableTypicalPointCompensation')&&cfg.enableTypicalPointCompensation;
phaseCorrection=ones(np,1);referenceBin=NaN;minimumVariance=NaN;referencePointFractions=nan(ns,1);
if doAutofocus
    criterion=filter([.5 .5],1,var(abs(work),0,1)./mean(abs(work).^2,1));
    [minimumVariance,referenceBin]=min(criterion);
    referencePointPower=reshape(sum(abs(aligned(:,referenceBin,:)).^2,1),ns,1);
    referencePointFractions=referencePointPower/sum(referencePointPower);
    phaseCorrection=conj(work(:,referenceBin)./abs(work(:,referenceBin)));
    corrected=feval('radar.phase_cmp',work,.12);
    assert(norm(corrected-work.*phaseCorrection,'fro')/norm(work,'fro')<1e-12);
    work=corrected;aligned=aligned.*phaseCorrection;
    replayError=norm(sum(aligned,3)-work,'fro')/max(norm(work,'fro'),eps);
    assert(replayError<1e-12&&all(isfinite(phaseCorrection)));
end
nDoppler=cfg.slowPadding*np;
alignedOnly=isfield(cfg,'computeOnlyAlignedMetrics')&&cfg.computeOnlyAlignedMetrics;
raw=[];
if ~alignedOnly,raw=fftshift(fft(scene,nDoppler,1),1)/np;end
compensated=fftshift(fft(work,nDoppler,1),1)/np;
fd=(-nDoppler/2:nDoppler/2-1)'*prf/nDoppler;
k=(cfg.fH-cfg.fL)/cfg.pulseWidthSec;
rangeDescending=-(-nr/2:nr/2-1)*(cfg.sampleRateHz/cfg.rangeFftPoints)* ...
    cfg.propagationSpeed/(2*k);
echoTime=pulseTime(:)+2*cfg.centerRangeM/cfg.propagationSpeed+cfg.pulseWidthSec/2;
fixedLoDiagnostic=isfield(cfg,'frameFixedLo');
cvDiagnostic=isfield(cfg,'frameMotionPrior');
if cvDiagnostic
    prior=cfg.frameMotionPrior;k=(cfg.fH-cfg.fL)/cfg.pulseWidthSec;
    acceleration=0;
    if isfield(prior,'accelerationMps2'),acceleration=prior.accelerationMps2;end
    dt=pulseTime(:)-prior.epochSec;
    Rtx=prior.rangeM+prior.velocityMps*dt+.5*acceleration*dt.^2;
    Vtx=prior.velocityMps+acceleration*dt;
    A=1-2*Vtx/cfg.propagationSpeed;B=2*Rtx/cfg.propagationSpeed-cfg.beatOffsetHz/k;
    delay=2*B./(A+sqrt(A.^2-4*acceleration/cfg.propagationSpeed*B));
    echoTime=pulseTime(:)+delay+cfg.pulseWidthSec/2;
end
if fixedLoDiagnostic
    fixed=cfg.frameFixedLo;
    echoTime=pulseTime(:)+fixed.loDelaySec+(round(cfg.sampleRateHz*cfg.pulseWidthSec)-1)/(2*cfg.sampleRateHz);
end
truth=scatterermodel.evaluateTarget(target,radar,echoTime);
diagnosticRanges=truth.rangeM;diagnosticDoppler=truth.dopplerHz;
if cvDiagnostic
    dt=echoTime-prior.epochSec;
    modelRange=prior.rangeM+prior.velocityMps*dt+.5*acceleration*dt.^2;
    modelVelocity=prior.velocityMps+acceleration*dt;
    diagnosticRanges=truth.rangeM-modelRange.';
    diagnosticDoppler=truth.dopplerHz+2*modelVelocity.'/radar.wavelengthM;
end
if fixedLoDiagnostic
    uc=(round(cfg.sampleRateHz*cfg.pulseWidthSec)-1)/(2*cfg.sampleRateHz);
    delta=2*truth.rangeM/cfg.propagationSpeed-fixed.loDelaySec;
    diagnosticDoppler=-2*truth.rangeRateMps/cfg.propagationSpeed.*(cfg.fL+k*(uc-delta))-fixed.commonLinearDopplerHz;
    if isfield(fixed,'modelSlowPhaseRemoved')&&fixed.modelSlowPhaseRemoved
        deltaModel=2*modelRange.'/cfg.propagationSpeed-fixed.loDelaySec;
        diagnosticDoppler=diagnosticDoppler+2*modelVelocity.'/cfg.propagationSpeed.*(cfg.fL+k*(uc-deltaModel));
    end
end
rows=struct([]);methods=1:2;if alignedOnly,methods=2;end
for method=methods
    if method==1,source=H;name="noRCMC";else,source=aligned;name="commonRCMC";end
    for j=1:ns
        A=abs(fftshift(fft(source(:,:,j),nDoppler,1),1));
        [~,idx]=max(A(:));[dopplerBin,column]=ind2sub(size(A),idx);
        s=source(:,column,j);
        nFine=64*np;line=abs(fftshift(fft(s,nFine)));
        [finePeak,peak]=max(line);line=circshift(line,nFine/2+1-peak);line=line/max(line);
        df=prf/nFine;relativeFrequency=(-nFine/2:nFine/2-1)'*df;
        m=signal.resolution_measure_3db(line,relativeFrequency,'linear',struct(),false);
        inside=relativeFrequency>=m.x_left&relativeFrequency<=m.x_right;
        localPeak=false(size(line));
        localPeak(2:end-1)=line(2:end-1)>=line(1:end-2)&line(2:end-1)>line(3:end);
        other=line(localPeak&~inside);
        otherDb=-Inf;if ~isempty(other),otherDb=20*log10(max(other));end
        [radius,order]=sort(abs(relativeFrequency));power=line.^2;
        cumulative=cumsum(power(order))/sum(power);
        energy90=2*radius(find(cumulative>=.9,1));
        [~,rangePeaks]=max(abs(source(:,:,j)),[],2);
        rangeValues=rangeDescending(rangePeaks);
        phase=unwrap(angle(s));fitPhase=polyval(polyfit(pulseTime(:),phase,1),pulseTime(:));
        trueFd=diagnosticDoppler(j,:);
        inputAlias=max(abs(trueFd))*cfg.fH/((cfg.fH+cfg.fL)/2)>=prf/2;
        if fixedLoDiagnostic,inputAlias=fixed.slowDopplerBoundHz(j)>=prf/2;end
        row=struct('method',name,'pointIndex',j,'peakRangeOffsetM',rangeDescending(column), ...
            'peakDopplerHz',fd(dopplerBin),'doppler3dBWidthHz',m.width, ...
            'idealDoppler3dBWidthHz',.885893/(np/prf), ...
            'widthToIdealRatio',m.width/(.885893/(np/prf)), ...
            'largestOtherPeakDb',otherDb,'mainlobeEnergyFraction',sum(line(inside).^2)/sum(line.^2), ...
            'energy90WidthAboutPeakHz',energy90, ...
            'coherentPeakPowerFraction',finePeak^2/(np*sum(abs(s).^2)), ...
            'residualHrrfPeakSpanM',max(rangeValues)-min(rangeValues), ...
            'geometricRangeSpanM',max(diagnosticRanges(j,:))-min(diagnosticRanges(j,:)), ...
            'fixedGateAmplitudeCV',std(abs(s))/max(mean(abs(s)),eps), ...
            'nonlinearPhaseRmsRad',sqrt(mean((phase-fitPhase).^2)), ...
            'truthDopplerMidHz',interp1(echoTime,trueFd,mean(echoTime),'linear'), ...
            'maxAbsTruthDopplerHz',max(abs(trueFd)), ...
            'frameAliasRiskAtHighBand',inputAlias,'inputPhaseAliasRisk',inputAlias, ...
            'phaseAutofocusApplied',doAutofocus&&method==2, ...
            'commonShiftSpanBins',max(shifts)-min(shifts));
        rows=[rows;row]; %#ok<AGROW>
    end
end
frame=struct('rawAmplitude',single(abs(raw(:,end:-1:1))), ...
    'alignedAmplitude',single(abs(compensated(:,end:-1:1))), ...
    'dopplerAxisHz',fd,'rangeAxisM',rangeDescending(end:-1:1), ...
    'metrics',struct2table(rows),'rcmcShiftBins',shifts, ...
    'replayError',replayError,'timeSec',mean(echoTime), ...
    'phaseCompensationEnabled',doAutofocus,'slowMeanRemovalEnabled',false, ...
    'phaseCorrection',phaseCorrection,'phaseReferenceBin',referenceBin, ...
    'phaseReferenceMinimumVariance',minimumVariance,'referencePointPowerFractions',referencePointFractions);
frame.externalSceneUsed=externalScene;frame.sceneDecompositionError=sceneError;
end

function out=shift_rows(H,shift)
[np,nr]=size(H);indices=mod((0:nr-1)+shift,nr)+1;
out=H(sub2ind([np nr],repmat((1:np)',1,nr),indices));
end
