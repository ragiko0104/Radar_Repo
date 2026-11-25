function [profiles,info]=isar_streamed_stretch_profiles(target,radar,pulseTime,cfg)
%ISAR_STREAMED_STRETCH_PROFILES Low-IF equivalent, continuous-delay dechirp.
% No sampled-real-RF mixing: analytic RF phase DIFFERENCE is evaluated at
% cfg.sampleRateHz (200 MHz). Negative beat offset preserves the native
% echo*conj(LO) Doppler sign. The known IF offset is removed in fast time.
% Output: pulse x cropped frequency bin x isolated point (complex single).
% Cubic Hermite range interpolation retains intra-pulse motion, using exact
% range/range-rate endpoints. Exact internal checkpoints bound its error.
% Only this caller-sized block of IF samples exists; no raw echo is saved.
t=double(pulseTime(:));np=numel(t);c=cfg.propagationSpeed;
fs=cfg.sampleRateHz;Tr=cfg.pulseWidthSec;k=(cfg.fH-cfg.fL)/Tr;
nFast=round(fs*Tr);nfft=cfg.rangeFftPoints;nr=cfg.rangeBins;
assert(nfft>=nFast&&mod(nr,2)==0&&nr<=nfft);
assert(all(radar.position_I==[0 0 0]),'This experiment uses a fixed radar at the origin.');
centre=target.centerStateFcn(target,radar,t);
Rc=vecnorm(centre.position_I,2,2);
assert(max(vecnorm(centre.velocity_I,2,2))<1e-12, ...
    'This independent experiment requires no target translation.');
tauLo=2*Rc/c-cfg.beatOffsetHz/k;
u=(0:nFast-1)/fs;span=u(end);z=u/span;
h00=2*z.^3-3*z.^2+1;h10=z.^3-2*z.^2+z;
h01=-2*z.^3+3*z.^2;h11=z.^3-z.^2;
start=t+tauLo;finish=start+span;
state=scatterermodel.evaluateTarget(target,radar,[start;finish]);
ns=size(state.rangeM,1);
R0=state.rangeM(:,1:np).'-Rc;R1=state.rangeM(:,np+1:end).'-Rc;
V0=state.rangeRateMps(:,1:np).';V1=state.rangeRateMps(:,np+1:end).';
% Continuous geometry checks are independent of the interpolation itself.
zc=[.113 .317 .573 .829];uc=zc*span;
tc=start+uc;exact=scatterermodel.evaluateTarget(target,radar,tc(:));
rangeError=0;maxBeat=0;
profiles=complex(zeros(np,nr,ns,'single'));
crop=nfft/2+1+(-nr/2:nr/2-1);
for j=1:ns
    rho=R0(:,j)*h00+span*V0(:,j)*h10+R1(:,j)*h01+span*V1(:,j)*h11;
    dTau=2*rho/c+cfg.beatOffsetHz/k;
    % Stable phase subtraction: never subtract two enormous RF phases.
    cycles=-dTau.*(cfg.fL+k*u)+.5*k*dTau.^2;
    gate=(u-dTau>=0)&(u-dTau<=Tr);
    alpha=target.payload.alpha0(j);
    IF=alpha*double(gate).*exp(1j*2*pi*cycles);
    % No fast-time DC rejection, no slow-time mean subtraction or autofocus.
    baseband=IF.*exp(1j*2*pi*cfg.beatOffsetHz*u);
    spectrum=fftshift(fft(baseband,nfft,2),2)/nFast;
    profiles(:,:,j)=single(spectrum(:,crop));
    zz=zc;
    rc=R0(:,j)*(2*zz.^3-3*zz.^2+1)+span*V0(:,j)*(zz.^3-2*zz.^2+zz)+ ...
        R1(:,j)*(-2*zz.^3+3*zz.^2)+span*V1(:,j)*(zz.^3-zz.^2);
    truth=reshape(exact.rangeM(j,:),np,numel(zc))-Rc;
    rangeError=max(rangeError,max(abs(rc-truth),[],'all'));
    rr=[state.rangeM(j,1:np).',reshape(exact.rangeM(j,:),np,numel(zc)), ...
        state.rangeM(j,np+1:end).'];
    vv=[V0(:,j),reshape(exact.rangeRateMps(j,:),np,numel(zc)),V1(:,j)];
    dd=2*(rr-Rc)/c+cfg.beatOffsetHz/k;
    fb=-k*dd-(2*vv/c).*(cfg.fL+k*([0 uc span]-dd));
    maxBeat=max(maxBeat,max(abs(fb),[],'all'));
end
assert(rangeError<1e-8,'Intrapulse interpolation exceeded 10 nm: %.3g m.',rangeError);
assert(maxBeat<cfg.beatGuardHz&&cfg.beatGuardHz<fs/2, ...
    'Beat exceeds the configured anti-alias guard.');
info=struct('maxRangeInterpolationErrorM',rangeError, ...
    'maxPhaseErrorBoundRad',4*pi*cfg.fH/c*rangeError, ...
    'maxAbsBeatHz',maxBeat,'loDelaySec',tauLo,'centerRangeM',Rc, ...
    'nFast',nFast,'phaseConvention',"echo*conj(LO); negative IF offset", ...
    'fastTimeModel',"200 MHz analytic difference phase; exact-endpoint cubic Hermite range");
end
