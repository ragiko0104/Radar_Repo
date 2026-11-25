%TEST_LEGACY_FOUR_SCATTERER_UNIFICATION Verify lib legacyDisk geometry and LOS projection.

clear;
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(thisFile));
run(fullfile(repoRoot,'startup_project.m'));

baseSpec = struct('Ntip',1,'radiusM',0.3,'phase0Rad',pi/6, ...
    'axialSeparationM',1.5,'alpha0',[1;0.3;0.3;0.3]);
motion = struct('type','accelerated','referenceTimeSec',0, ...
    'p0_I',[10050,0,0],'v0_I',[3468,0,0],'a0_I',[-150,0,0], ...
    'attitude0YprRad',[-3*pi/4,0,0]);
radarState = struct('position_I',[0,0,0],'velocity_I',[0,0,0], ...
    'carrierFrequencyHz',14.8e9,'dopplerSign',-1);
pulseTime = (0:31).'*1e-3;

for nTip = [1,0]
    spec = baseSpec;
    spec.Ntip = nTip;
    [target,sharedInfo] = scatterermodel.createLegacyDiskTarget(spec,motion);
    target = scatterermodel.withConstantBodyRate(target,[0,2*pi,0],0);
    radarPoints = target.payload.positionTemplate_B.';
    assert(isequal(size(radarPoints),[3,4]));
    assert(max(abs(sharedInfo.alpha0-baseSpec.alpha0),[],'all') < 1e-12);

    if nTip == 1
        expectedDefault = [0,0.3*cos(pi/6),0.3*cos(5*pi/6),0.3*cos(3*pi/2); ...
            -0.75,0.75,0.75,0.75; ...
            0,0.3*sin(pi/6),0.3*sin(5*pi/6),0.3*sin(3*pi/2)];
        assert(max(abs(radarPoints-expectedDefault),[],'all') < 1e-12);
    end

    state = scatterermodel.evaluateTarget(target,radarState,pulseTime,1:4);
    initialRotation = state.kinematics.C_IB(:,:,1);
    expectedRotation = [-sqrt(0.5),sqrt(0.5),0; ...
        -sqrt(0.5),-sqrt(0.5),0; 0,0,1];
    assert(max(abs(initialRotation-expectedRotation),[],'all') < 1e-12);
    bodyAxisLosAngles = acosd(max(-1,min(1,[1,0,0]*initialRotation)));
    assert(max(abs(bodyAxisLosAngles-[135,45,90])) < 1e-10);

    oldLos = [0;1;0];
    oldBodyY = [sqrt(0.5);sqrt(0.5);0];
    oldBodyZ = [0;0;1];
    oldBodyX = cross(oldBodyY,oldBodyZ);
    oldRotation = [oldBodyX,oldBodyY,oldBodyZ];
    oldProjected = radar.ext_target_rot_pos_cal_opt( ...
        radarPoints,1e-3,32,[0;2*pi;0],oldLos,oldRotation);
    newProjected = zeros(32,4);
    for pulseIndex = 1:32
        newProjected(pulseIndex,:) = [1,0,0]* ...
            state.kinematics.C_IB(:,:,pulseIndex)*radarPoints;
    end
    assert(max(abs(newProjected-oldProjected),[],'all') < 1e-12, ...
        'Refactored +x-LOS projection must equal the initial +y-LOS code.');

    pairIndex = nchoosek(1:4,2);
    pairRms = zeros(size(pairIndex,1),1);
    for pairNumber = 1:size(pairIndex,1)
        difference = state.rangeM(pairIndex(pairNumber,1),:)- ...
            state.rangeM(pairIndex(pairNumber,2),:);
        pairRms(pairNumber) = sqrt(mean(difference.^2));
    end
    assert(all(pairRms > 1e-3),'Each scatterer pair must have a distinct range track.');
end

fprintf('Legacy four-scatterer lib unification test passed.\n');
