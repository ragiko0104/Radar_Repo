function test_isar_micro_motion()
% Analytic state parity, centre-motion preservation, and physical scaling.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'functions')); ensure_md_zt_paths();
spec=struct('Ntip',1,'radiusM',.3,'phase0Rad',pi/6, ...
    'axialSeparationM',1.5,'alpha0',[1;.3;.3;.3]);
motion=struct('type','accelerated','referenceTimeSec',.1, ...
    'p0_I',[10050,0,0],'v0_I',[3468,0,0],'a0_I',[-150,0,0], ...
    'attitude0YprRad',[-3*pi/4,0,0]);
target=scatterermodel.createLegacyDiskTarget(spec,motion);
radar=struct('position_I',[0 0 0],'velocity_I',[0 0 0], ...
    'carrierFrequencyHz',14.8e9,'wavelengthM',3e8/14.8e9,'dopplerSign',-1);
cfg=struct('bodyRotationMode',"spin",'spinFrequencyHz',1, ...
    'precessionFrequencyHz',1,'precessionHalfAngleDeg',10, ...
    'startTimeSec',.2,'propagationSpeed',3e8,'Tr',100e-6, ...
    'lambda',3e8/14.8e9,'Ta',1e-3);
t=[.2;.200001;.2155;.231;.4];
for mode=["spin","coneSpin","swingSpin","coneSwingSpin"]
    cfg.bodyRotationMode=mode;
    for fs=[.5 1 2 3 3.5]
        cfg.spinFrequencyHz=fs;cfg.precessionFrequencyHz=4-fs;
        cfg.swingFrequencyHz=.5;
        fast=isar_apply_micro_motion(target,cfg);
        generic=fast;generic.evaluateFcn=fast.metadata.isarOriginalEvaluateFcn;
        for selection={1:4,[4 2],1}
            a=scatterermodel.evaluateTarget(fast,radar,t,selection{1});
            b=scatterermodel.evaluateTarget(generic,radar,t,selection{1});
            for field={'rangeM','position_I','velocity_I','rangeRateMps','dopplerHz','alpha'}
                difference=max(abs(a.(field{1})-b.(field{1})),[],'all');
                assert(difference<2e-8,'Vectorized state mismatch: %s',field{1});
            end
        end
        oldKin=scatterermodel.evaluateRigidMotion(struct(),motion,t);
        assert(max(abs(a.centerPosition_I-oldKin.position_I),[],'all')<1e-10);
        assert(max(abs(a.centerVelocity_I-oldKin.velocity_I),[],'all')<1e-10);
        if mode=="spin"
            old=scatterermodel.withConstantBodyRate(target,[0,2*pi*fs,0],cfg.startTimeSec);
            oldState=scatterermodel.evaluateTarget(old,radar,t);
            newState=scatterermodel.evaluateTarget(fast,radar,t);
            assert(max(abs(oldState.rangeM-newState.rangeM),[],'all')<1e-9);
            % Pure spin does not move the tip relative to the body origin.
            tip=reshape(newState.position_I(1,:,:),3,[]).'-newState.centerPosition_I;
            assert(max(abs(tip-tip(1,:)),[],'all')<1e-9);
            g=isar_cross_range_geometry(fast,radar,cfg,.2+(0:31)'*.001);
            assert(abs(g.effectiveRateRadps-2*pi*fs/sqrt(2))<1e-10);
            assert(g.maximumPlaneDriftDeg<2e-6);
            assert(abs(g.nominalResolutionM*fs-.0712877150673384)<1e-10);
        end
    end
end
% Changing LOS contributes apparent rotation with the correct sign.
cfg.bodyRotationMode="spin";cfg.spinFrequencyHz=1;
moving=target;moving.payload.motion.v0_I=[3468,100,0];
fast=isar_apply_micro_motion(moving,cfg);
g=isar_cross_range_geometry(fast,radar,cfg,.2+(0:31)'*.001);
centre=fast.centerStateFcn(fast,radar,g.referenceTimeSec);
n=centre.position_I/norm(centre.position_I);
nDot=(centre.velocity_I-dot(centre.velocity_I,n)*n)/norm(centre.position_I);
expected=norm(cross(n,centre.angularVelocity_I-cross(n,nDot)));
assert(abs(g.effectiveRateRadps-expected)<1e-12);
% Each target gets its own scale; a shared angular-speed norm is incorrect.
twoTargets=[target,target];
twoTargets(2).payload.motion.attitude0YprRad=[0 0 0];
twoTargets=isar_apply_micro_motion(twoTargets,cfg);
twoGeometry=isar_cross_range_geometry(twoTargets,radar,cfg,.2+(0:31)'*.001);
assert(abs(twoGeometry(1).effectiveRateRadps-2*pi/sqrt(2))<1e-10);
assert(abs(twoGeometry(2).effectiveRateRadps-2*pi)<1e-10);
% Library translation aliases and custom paths survive attitude replacement.
alias=target;
alias.payload.motion=struct('motionType','accelerated','t0Sec',.1, ...
    'position0_I',[10050 0 0],'velocityMps',[3468 2 -1], ...
    'accelerationMps2',[-150 1 0],'attitude0YprRad',[-3*pi/4 0 0]);
custom=target;
custom.payload.motion=struct('type','custom', ...
    'positionFcn',@(x) [10050+3*x(:).^2,2*sin(x(:)),zeros(numel(x),1)], ...
    'velocityFcn',@(x) [6*x(:),2*cos(x(:)),zeros(numel(x),1)], ...
    'accelerationFcn',@(x) [6*ones(numel(x),1),-2*sin(x(:)),zeros(numel(x),1)], ...
    'attitudeFcn',@(x) repmat(eye(3),1,1,numel(x)), ...
    'omegaBodyFcn',@(x) zeros(numel(x),3));
for original=[alias,custom]
    converted=isar_apply_micro_motion(original,cfg);
    before=scatterermodel.evaluateRigidMotion(struct(),original.payload.motion,t);
    after=converted.centerStateFcn(converted,radar,t);
    assert(max(abs(before.position_I-after.position_I),[],'all')<1e-10);
    assert(max(abs(before.velocity_I-after.velocity_I),[],'all')<1e-10);
    assert(max(abs(before.acceleration_I-after.acceleration_I),[],'all')<1e-10);
end
% A spin axis along LOS cannot be calibrated as a spatial cross-range axis.
degenerate=target;degenerate.payload.motion.attitude0YprRad=[-pi/2 0 0];
degenerate=isar_apply_micro_motion(degenerate,cfg);failed=false;
try,isar_cross_range_geometry(degenerate,radar,cfg,.2+(0:31)'*.001);
catch problem,failed=strcmp(problem.identifier,'RadarProject:ISARMotion:UnobservableAspect');end
assert(failed);
% Invalid controls fail before allocating echo arrays.
for bad=[-.1,Inf,NaN]
    cfg.spinFrequencyHz=bad;failed=false;
    try,isar_apply_micro_motion(target,cfg);catch,failed=true;end
    assert(failed);
end
cfg.spinFrequencyHz=1;cfg.bodyRotationMode="notAMode";failed=false;
try,isar_apply_micro_motion(target,cfg);catch,failed=true;end
assert(failed);
% Exercise the real fullpath configuration block without synthesizing its
% ten-target scene; the public caller override must reach cfg_ISARpath.
ISAR_MOTION_OPTIONS=struct('Mode',"coneSpin",'SpinFrequencyHz',1, ...
    'PrecessionFrequencyHz',2,'PrecessionHalfAngleDeg',10); %#ok<NASGU>
fullpathScript=which('multiuav_iasr_fullpath.m');
assert(~isempty(fullpathScript),'RV_ISAR_Path/code/multiuav_iasr_fullpath.m not on path.');
source=fileread(fullpathScript);
configBlock=extractBefore(string(source),"%% 其余仿真参数定义");
eval(char(configBlock));
assert(cfg_ISARpath.bodyRotationMode=="coneSpin");
assert(cfg_ISARpath.spinFrequencyHz==1&&cfg_ISARpath.precessionFrequencyHz==2);
assert(cfg_ISARpath.precessionHalfAngleDeg==10);
assert(cfg_ISARpath.Na==32&&cfg_ISARpath.Ta==.001);
fprintf('ISAR micro-motion analytic parity and scaling tests passed.\n');
end
