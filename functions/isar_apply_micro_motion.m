function [targetsOut, settings] = isar_apply_micro_motion(targets, cfg)
%ISAR_APPLY_MICRO_MOTION Preserve translation and select spin/coning/swing.
% Modes: spin, coneSpin, swingSpin, coneSwingSpin. Nonnegative Hz rates.
% Swing is the library's inertial-Z periodic rigid-body swing, not a
% time-varying coning half angle. Optional swingFrequencyHz=1,
% swingAngleMinDeg=-5, swingAngleMaxDeg=5. Its initial angle is the minimum.
% coneSpin uses C(t)=C0*Ry(wp*t)*Rz(alpha)*Ry(ws*t), exactly the convention
% of scatterermodel.createConicalMicroMotion. ws is the final Euler-factor
% rate: omega dot the symmetry axis is ws+wp*cos(alpha), not simply ws.
% Rotation is about the body origin (the original cone midpoint), not a
% fixed tip pivot. The tip is on the spin axis and stays fixed in spin mode.
% Legacy callers with no new mode, or constantBodyRate, retain omegaBody_B.

mode = string(field_value(cfg,'bodyRotationMode',"constantBodyRate"));
if ~isscalar(mode) || ~ismember(lower(mode), ...
        ["spin","conespin","swingspin","coneswingspin","constantbodyrate"])
    error('RadarProject:ISARMotion:InvalidMode', ...
        'Use spin, coneSpin, swingSpin or coneSwingSpin.');
end
reference = field_value(cfg,'startTimeSec',0);
validateattributes(reference,{'numeric'},{'real','finite','scalar'});
settings = struct('mode',mode,'spinFrequencyHz',NaN, ...
    'precessionFrequencyHz',0,'precessionHalfAngleDeg',0, ...
    'rotationReferenceTimeSec',reference,'pivot',"bodyOrigin");
if strcmpi(mode,"constantBodyRate")
    targetsOut = scatterermodel.withConstantBodyRate(targets,cfg.omegaBody_B,reference);
    settings.spinFrequencyHz = norm(cfg.omegaBody_B)/(2*pi);
    return;
end
fs = field_value(cfg,'spinFrequencyHz',1);
fp = field_value(cfg,'precessionFrequencyHz',1);
angleDeg = field_value(cfg,'precessionHalfAngleDeg',10);
swingHz = field_value(cfg,'swingFrequencyHz',1);
swingMin = field_value(cfg,'swingAngleMinDeg',-5);
swingMax = field_value(cfg,'swingAngleMaxDeg',5);
validateattributes(fs,{'numeric'},{'real','finite','scalar','nonnegative'});
validateattributes(fp,{'numeric'},{'real','finite','scalar','nonnegative'});
validateattributes(angleDeg,{'numeric'},{'real','finite','scalar','>=',0,'<',90});
validateattributes(swingHz,{'numeric'},{'real','finite','scalar','nonnegative'});
validateattributes(swingMin,{'numeric'},{'real','finite','scalar','>',-90,'<',90});
validateattributes(swingMax,{'numeric'},{'real','finite','scalar','>=',swingMin,'<',90});
if ~contains(lower(mode),"cone"), fp=0; angleDeg=0; end
if ~contains(lower(mode),"swing"), swingHz=0; swingMin=0; swingMax=0; end
settings.spinFrequencyHz=fs;
settings.precessionFrequencyHz=fp;
settings.precessionHalfAngleDeg=angleDeg;
settings.axialRotationFrequencyHz=fs+fp*cosd(angleDeg);
settings.swingFrequencyHz=swingHz;
settings.swingAngleMinDeg=swingMin;
settings.swingAngleMaxDeg=swingMax;
settings.axialRotationFrequencyNote= ...
    "spin+precession projection only; swing adds a time-dependent projection";
targetsOut=targets;
for k=1:numel(targets)
    target=scatterermodel.validateTarget(targets(k));
    if strcmpi(target.kind,'projectedRange')
        error('RadarProject:ISARMotion:ProjectedRange','Micro-motion requires 3-D geometry.');
    end
    old=target.payload.motion;
    % This entry point replaces attitude only, as withConstantBodyRate did.
    % Preserve the full original analytic centre trajectory, including custom
    % motion. Common polynomial translations take a vectorized fast path.
    initial=scatterermodel.evaluateRigidMotion(struct(),old,reference);
    motionConfig=struct('Mode',mode,'PositionM',initial.position_I, ...
        'InitialBodyToInertial',initial.C_IB,'ReferenceTimeSec',reference, ...
        'SpinRateRadps',2*pi*fs,'ConingRateRadps',2*pi*fp, ...
        'ConingHalfAngleRad',deg2rad(angleDeg), ...
        'SwingRateRadps',2*pi*swingHz, ...
        'SwingAngleMinRad',deg2rad(swingMin),'SwingAngleMaxRad',deg2rad(swingMax));
    [motion,truth]=scatterermodel.createConicalMicroMotion(motionConfig);
    motion.positionFcn=@(t) original_translation(old,t,'position_I');
    motion.velocityFcn=@(t) original_translation(old,t,'velocity_I');
    motion.accelerationFcn=@(t) original_translation(old,t,'acceleration_I');
    target.payload.motion=motion;
    target.metadata.isarMicroMotion=settings;
    target.metadata.isarMicroMotion.truth=truth;
    % The generic rigid evaluator loops over every 1 GHz sample. Evaluate
    % the same analytic state in arrays; UAV evaluators keep their rotors.
    if ismember(lower(string(target.kind)),["legacydisk","rigidpointcloud","point"])
        target.metadata.isarOriginalEvaluateFcn=target.evaluateFcn;
        target.evaluateFcn=@isar_rigid_micro_motion_state;
        target.centerStateFcn=@(targetValue,~,t) ...
            isar_rigid_micro_motion_state(targetValue,struct(),t,[]);
    end
    targetsOut(k)=scatterermodel.validateTarget(target);
end
end

function value=original_translation(motion,t,name)
t=double(t(:));
mode=lower(string(first_field(motion,{'type','mode','motionType'},'hover')));
vibration=field_value(motion,'vibration',struct());
hasVibration=isfield(vibration,'enabled')&&vibration.enabled;
if ismember(mode,["hover","stationary","fixed","linear","constantvelocity", ...
        "constant_velocity","accelerated","constantacceleration", ...
        "constant_acceleration","constantbodyrate","constant_body_rate","bodyrate","body_rate"]) && ~hasVibration
        dt=t-first_field(motion,{'referenceTimeSec','timeReferenceSec','t0Sec'},0);
        p=first_field(motion,{'p0_I','position0_I','initialPosition_I','initialPositionM','positionM'},[0 0 0]);
        v=first_field(motion,{'v0_I','velocity0_I','initialVelocity_I','velocityMps','velocity_I'},[0 0 0]);
        a=first_field(motion,{'a0_I','acceleration0_I','accelerationMps2','acceleration_I'},[0 0 0]);
        p=reshape(p,1,3);v=reshape(v,1,3);a=reshape(a,1,3);
        if ismember(mode,["hover","stationary","fixed"]),v=[0 0 0];a=v;end
        if ismember(mode,["linear","constantvelocity","constant_velocity"]),a=[0 0 0];end
        switch name
            case 'position_I',value=p+dt.*v+.5*dt.^2.*a;
            case 'velocity_I',value=v+dt.*a;
            case 'acceleration_I',value=repmat(a,numel(t),1);
        end
        return;
end
kin=scatterermodel.evaluateRigidMotion(struct(),motion,t);
value=kin.(name);
end

function value=first_field(s,names,fallback)
value=fallback;
for k=1:numel(names)
    if isfield(s,names{k})&&~isempty(s.(names{k})),value=s.(names{k});return;end
end
end

function value=field_value(s,name,fallback)
if isfield(s,name)&&~isempty(s.(name)),value=s.(name);else,value=fallback;end
end
