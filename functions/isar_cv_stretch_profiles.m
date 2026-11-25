function [profiles,info]=isar_cv_stretch_profiles(target,radar,pulseTime,cfg,prior)
%ISAR_CV_STRETCH_PROFILES Streaming low-IF synthesis and CV-only processing.
% PRIOR supplies epochSec, rangeM and velocityMps, scalar or one per pulse.
% No acceleration is accepted by this processor. Target truth is used ONLY
% to synthesize the received waveform and verify interpolation / beat bounds.
% The pulse LO follows the CV range. A conjugate CV phase removes intrapulse
% dispersion, range-Doppler coupling and known translation coherently. Exact
% digital frequency translation followed by FFT cropping implements the
% model-centred frequency window without rounding a fractional FFT bin.
t=double(pulseTime(:));np=numel(t);c=cfg.propagationSpeed;
epoch=expand(prior.epochSec,np);Rp=expand(prior.rangeM,np);vp=expand(prior.velocityMps,np);
assert(~isfield(prior,'accelerationMps2'),'Acceleration must not enter the processor.');
fs=cfg.sampleRateHz;Tr=cfg.pulseWidthSec;k=(cfg.fH-cfg.fL)/Tr;
nFast=round(fs*Tr);nfft=cfg.rangeFftPoints;nr=cfg.rangeBins;
assert(nfft>=nFast&&mod(nr,2)==0&&nr<=nfft);
Rtx=Rp+vp.*(t-epoch);d0=cfg.beatOffsetHz/k;
% Solve tauLO=2*R_CV(t+tauLO)/c-d0 using ONLY the supplied CV prior.
tauLo=(2*Rtx/c-d0)./(1-2*vp/c);
u=(0:nFast-1)/fs;span=u(end);z=u/span;
h00=2*z.^3-3*z.^2+1;h10=z.^3-2*z.^2+z;
h01=-2*z.^3+3*z.^2;h11=z.^3-z.^2;
start=t+tauLo;state=scatterermodel.evaluateTarget(target,radar,[start;start+span]);
ns=size(state.rangeM,1);reference=Rp+vp.*(start-epoch);
R0=state.rangeM(:,1:np).'-reference;R1=state.rangeM(:,np+1:end).'-reference;
V0=state.rangeRateMps(:,1:np).';V1=state.rangeRateMps(:,np+1:end).';
zc=[.113 .317 .573 .829];uc=zc*span;tc=start+uc;
exact=scatterermodel.evaluateTarget(target,radar,tc(:));
rangeError=0;maxBeat=0;maxProcessedBeat=0;
profiles=complex(zeros(np,nr,ns,'single'));
crop=nfft/2+1+(-nr/2:nr/2-1);
dModel=d0+2*vp/c.*u;
modelCycles=-dModel.*(cfg.fL+k*u)+.5*k*dModel.^2;
nominalCycles=-d0*(cfg.fL+k*u)+.5*k*d0^2;
beforeWidth=nan(1,ns);afterWidth=beforeWidth;
for j=1:ns
    rho=R0(:,j)*h00+span*V0(:,j)*h10+R1(:,j)*h01+span*V1(:,j)*h11;
    delta=d0+2*rho/c;
    actualCycles=-delta.*(cfg.fL+k*u)+.5*k*delta.^2;
    gate=(u-delta>=0)&(u-delta<=Tr);
    % Difference of phases, not separately sampled GHz RF signals. The known
    % negative IF is shifted exactly to zero before the cropped range FFT.
    cycles=actualCycles-modelCycles+nominalCycles+cfg.beatOffsetHz*u;
    beat=target.payload.alpha0(j)*double(gate).*exp(1j*2*pi*cycles);
    spectrum=fftshift(fft(beat,nfft,2),2)/nFast;
    profiles(:,:,j)=single(spectrum(:,crop));
    raw=target.payload.alpha0(j)*double(gate(1,:)).*exp(1j*2*pi*actualCycles(1,:));
    beforeWidth(j)=range_width(raw,fs,nfft,k,c);
    afterWidth(j)=range_width(beat(1,:),fs,nfft,k,c);
    rc=R0(:,j)*(2*zc.^3-3*zc.^2+1)+span*V0(:,j)*(zc.^3-2*zc.^2+zc)+ ...
        R1(:,j)*(-2*zc.^3+3*zc.^2)+span*V1(:,j)*(zc.^3-zc.^2);
    truth=reshape(exact.rangeM(j,:),np,numel(zc))-reference;
    rangeError=max(rangeError,max(abs(rc-truth),[],'all'));
    rr=[R0(:,j),truth,R1(:,j)];
    vv=[V0(:,j),reshape(exact.rangeRateMps(j,:),np,numel(zc)),V1(:,j)];
    uf=[0 uc span];dd=d0+2*rr/c;
    fb=-k*dd-(2*vv/c).*(cfg.fL+k*(uf-dd));
    dm=d0+2*vp/c.*uf;
    fmodel=-k*dm-(2*vp/c).*(cfg.fL+k*(uf-dm));
    maxBeat=max(maxBeat,max(abs(fb),[],'all'));
    maxProcessedBeat=max(maxProcessedBeat,max(abs(fb-fmodel),[],'all'));
end
assert(rangeError<1e-8,'Intrapulse range interpolation exceeded 10 nm.');
assert(maxBeat<cfg.beatGuardHz&&cfg.beatGuardHz<fs/2,'Raw beat exceeds ADC guard.');
assert(maxProcessedBeat<(nr/2-8)*fs/nfft,'Model-centred range window loses target.');
info=struct('maxRangeInterpolationErrorM',rangeError,'maxAbsBeatHz',maxBeat, ...
    'maxProcessedBeatHz',maxProcessedBeat,'loDelaySec',tauLo, ...
    'rangeWidthBeforeDispersionM',beforeWidth,'rangeWidthAfterDispersionM',afterWidth, ...
    'frequencyWindowCenterBeforeTranslationHz',-cfg.beatOffsetHz, ...
    'knownTranslationPhaseRemoved',true,'accelerationPriorUsed',false);
end

function x=expand(x,n)
x=double(x(:));if isscalar(x),x=repmat(x,n,1);end
assert(numel(x)==n&&all(isfinite(x)));
end

function w=range_width(x,fs,nfft,k,c)
nfft=4*nfft; % Resolve the 2.2 cm mainlobe beyond the 1.53 cm display grid.
a=abs(fftshift(fft(x,nfft)));[~,p]=max(a);
a=circshift(a,nfft/2+1-p);a=a/max(a);
f=(-nfft/2:nfft/2-1)*fs/nfft*c/(2*k);
m=signal.resolution_measure_3db(a,f,'linear',struct(),false);w=m.width;
end
