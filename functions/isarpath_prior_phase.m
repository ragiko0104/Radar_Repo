function reference=isarpath_prior_phase(cfg,prior,t,p,shifts,mode)
%ISARPATH_PRIOR_PHASE Kinematic phase in the EXACT original processing gauge.
% Evaluate the known positive-IF phase at the prior receive-aperture midpoint,
% then subtract the original dispersion phase and moving FFT-window phase.
% Unlike taking phase from a synthetic HRRP gate, this definition cannot
% acquire pi jumps when the mismatched CV reference crosses a sinc zero.
% No scatterer/micro-motion truth or measured scene phase enters this helper.
acc=0;if string(mode)=="MODEL_CA",acc=prior.accelerationMps2;end
dt=t(:)-prior.epochSec;lo=round(prior.loStartSec*cfg.fs)/cfg.fs;
starts=round((prior.loStartSec-p.timeMarginSamples/cfg.fs_scope)*cfg.fs_scope)+ones(size(dt));
if p.timeWindowMode=="MOVING_PULSE"
    starts=round((prior.loStartSec+2*prior.velocityMps*dt/cfg.propagationSpeed- ...
        p.timeMarginSamples/cfg.fs_scope)*cfg.fs_scope)+1;
end
segmentStart=(starts-1)/cfg.fs_scope;
echoStart=lo+zeros(size(dt));echoStop=echoStart+cfg.Tr;
for iteration=1:4
    echoStart=2*range_at(dt+echoStart)/cfg.propagationSpeed;
    echoStop=cfg.Tr+2*range_at(dt+echoStop)/cfg.propagationSpeed;
end
left=max(max(echoStart,lo),segmentStart);
right=min(min(echoStop-1/cfg.fs,lo+cfg.Tr-1/cfg.fs), ...
    segmentStart+(p.timeWindowSamples-1)/cfg.fs_scope);
assert(all(right>left),'Prior-predicted echo has no LO overlap.');
fast=(left+right)/2;local=fast-segmentStart;
R=range_at(dt+fast);velocity=prior.velocityMps+acc*(dt+fast);
D=2*R/cfg.propagationSpeed-lo;u=fast-lo;
positiveCycles=D.*(cfg.fL+cfg.kr*u)-.5*cfg.kr*D.^2;
dispersionCycles=zeros(size(dt));dispersionHz=zeros(size(dt));
if p.dispersionCompensationEnabled
    dispersionCycles=2*prior.velocityMps/cfg.propagationSpeed*(cfg.kr*local.^2+cfg.fL*local);
    dispersionHz=2*prior.velocityMps/cfg.propagationSpeed*(2*cfg.kr*local+cfg.fL);
end
frequency=prior.beatOffsetHz+2*cfg.kr*prior.velocityMps*dt/cfg.propagationSpeed;
centres=round(frequency/p.frequencyBinSpacingHz)+1+p.cropCenterOffsetBinsEffective;
fftFrequency=(centres-1+shifts(:))*p.frequencyBinSpacingHz;
instantaneousHz=cfg.kr*D+2*velocity/cfg.propagationSpeed.*(cfg.fL+cfg.kr*(fast-2*R/cfg.propagationSpeed))-dispersionHz;
% A fixed virtual sub-bin range gate is centred on the model's FIRST pulse.
% Its offset remains fixed for the entire CPI, just as a real HRRP column.
fixedGateOffsetHz=instantaneousHz(1)-fftFrequency(1);
cycles=positiveCycles-dispersionCycles-(fftFrequency+fixedGateOffsetHz).*local;
correction=exp(-1j*2*pi*(cycles-cycles(1)));
reference=struct('phaseCorrection',correction,'phaseCycles',cycles, ...
    'fixedGateOffsetHz',fixedGateOffsetHz,'apertureMidpointSec',fast, ...
    'method',"prior polynomial + positive IF + dispersion + crop/alignment coordinates", ...
    'accelerationUsedMps2',acc);

    function R=range_at(x)
        R=prior.rangeM+prior.velocityMps*x+.5*acc*x.^2;
    end
end
