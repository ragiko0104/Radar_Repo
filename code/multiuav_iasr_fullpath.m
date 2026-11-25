% 测速测距 + ISAR 成像全链路仿真（基于 matlab_research_lib）
% created by lms in 2026/07/13
if isempty(which('radar.echo_generate_target_lfm'))
    run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'startup_project.m'));
end
%% 常量定义
c = 3e8;
Z = 50;                                         % 50Ω匹配负载
q = 1.6e-19;                                    % 单位电荷
lc = 1550e-9;       nu = c/lc;                  % 使用的光波长λ= 1550nm
h = 6.62e-34;                                   % Planck常数h = 6.62e-34 J·S
kB = 1.38e-23;                                  % Boltzmann常数kB = 1.38e-23 J/K
T_exp = 300;                                    % 实验室温度T_exp = 27℃ = 300K
Mach = 340;
g = 10;
win_G = 1;
%% 开关定义（允许调用方在工作区预先覆盖）
if ~exist('run_RVpath', 'var') || isempty(run_RVpath)
    run_RVpath = false;
end
if ~exist('run_ISARpath', 'var') || isempty(run_ISARpath)
    run_ISARpath = true;
end
if ~exist('UAV_PLOT_GEOMETRY_ENABLED', 'var') || isempty(UAV_PLOT_GEOMETRY_ENABLED)
    UAV_PLOT_GEOMETRY_ENABLED = true;
end
%% 目标运动参数定义
%% 测量波形参数定义
% fs = 64e9;
cfg_RVpath = struct();
cfg_RVpath.fs = 200e6;
cfg_RVpath.Tr = 50e-6;
cfg_RVpath.T_space = 5e-6;
cfg_RVpath.T = 1.2e-3;
cfg_RVpath.Ta = cfg_RVpath.T;
cfg_RVpath.fc = 12.5e9;
cfg_RVpath.f_space = 5e6;
cfg_RVpath.B = [5, 7.5, 10]*1e6;
cfg_RVpath.f1_up = 10e6;
cfg_RVpath.Na = 1;
if cfg_RVpath.Ta < cfg_RVpath.T
    error("In RVpath, Ta should be at least the same with T or larger\n");
end
% R/V cross-band consistency gate. The legacy optimized chain used fixed
% gates c/min(B)/5 and 2/5*c/(2*Tr*fc), which are 12 m and 96 m/s for the
% current waveform. The new default estimates each refined peak's effective
% SINR from its local CFAR training cells, converts it to a delay-precision standard
% deviation, and propagates the three delay uncertainties through the two
% Cal_rv_mat solutions before matching them.
cfg_RVpath.rvConsistency = struct( ...
    'mode','snrPrecision', ...                 % snrPrecision / legacyResolution
    'sigmaMultiplier',2, ...                  % two-sided consistency gate in sigma
    'snrFloorDb',0, ...                       % numerical lower bound
    'snrCeilingDb',60, ...                    % finite gate for noiseless simulation
    'delayStdFloorBins',0.05, ...              % peak-fit/model precision floor
    'noiseEstimator','localCfar', ...          % localCfar effective SINR / robustOffTarget
    'biasCorrection','rankMedian', ...         % none / rankMedian
    'maxPeakRankDifference',0, ...             % equal-range-order association
    'enforceOneToOnePeakUse',true);            % prevent combinatorial duplicates
%% 成像波形参数定义
cfg_ISARpath = struct();
cfg_ISARpath.fs = 1e9;
cfg_ISARpath.Tr = 100e-6;
cfg_ISARpath.T = 1e-3;
cfg_ISARpath.Ta = 1e-3;
cfg_ISARpath.fL = 11.8e9;
cfg_ISARpath.fH = 17.8e9;
cfg_ISARpath.Na = 32;
if exist('FULLPATH_OVERRIDES', 'var') && isstruct(FULLPATH_OVERRIDES) && ...
        isfield(FULLPATH_OVERRIDES, 'isarNa') && ~isempty(FULLPATH_OVERRIDES.isarNa)
    cfg_ISARpath.Na = double(FULLPATH_OVERRIDES.isarNa);
end
cfg_ISARpath.fs_scope = 200e6;                  % FPGA/示波器侧去斜信号采样率
cfg_ISARpath.echoMode = "fastTimeExact";        % 在快时间内更新动态散射点位置
if exist('FULLPATH_OVERRIDES', 'var') && isstruct(FULLPATH_OVERRIDES) && ...
        isfield(FULLPATH_OVERRIDES, 'isarEchoMode') && ~isempty(FULLPATH_OVERRIDES.isarEchoMode)
    cfg_ISARpath.echoMode = string(FULLPATH_OVERRIDES.isarEchoMode);
end
cfg_ISARpath.propagationSpeed = c;
cfg_ISARpath.B = abs(cfg_ISARpath.fH-cfg_ISARpath.fL);
cfg_ISARpath.kr = (cfg_ISARpath.fH-cfg_ISARpath.fL)/cfg_ISARpath.Tr;
cfg_ISARpath.fc = (cfg_ISARpath.fH+cfg_ISARpath.fL)/2;
cfg_ISARpath.lambda = cfg_ISARpath.propagationSpeed/cfg_ISARpath.fc;
cfg_ISARpath.startTimeSec = 0;                  % 成像孔径从绝对时刻 t=0 s 开始
cfg_ISARpath.bodyRotationMode = "constantBodyRate";
cfg_ISARpath.omegaBody_B = [0, 2*pi, 0];       % legacy：绕机体系 +y 轴以 2*pi rad/s 匀速转动
if cfg_ISARpath.Ta < cfg_ISARpath.T
    error("In ISARpath, Ta should be at least the same with T or larger\n");
end
%% 其余仿真参数定义
% R0 = (2:11)'*100;
% v_Mach = ones(10, 1);
% a_g = zeros(10, 1);
R0 = [10050;20050;30050;40050;50050;65050;75050;85050;100050;120050];
v_Mach = ones(10, 1)*10.2*sqrt(2);
a_g = [-15, -12, -11, -10.8, -10.6, -10.3, -10.2, -10.1, -10.1, -10.1]'*sqrt(2);
if exist('FULLPATH_OVERRIDES', 'var') && isstruct(FULLPATH_OVERRIDES) && ...
        isfield(FULLPATH_OVERRIDES, 'targetIndices') && ~isempty(FULLPATH_OVERRIDES.targetIndices)
    idxKeep = FULLPATH_OVERRIDES.targetIndices(:);
    R0 = R0(idxKeep);
    v_Mach = v_Mach(idxKeep);
    a_g = a_g(idxKeep);
end
n_Targets = numel(R0);
%% 目标散射体模型加载
% 每个物理目标都将在 motion 定义后转换为统一 scatterermodel target。
% targetModelTypes 可以是单个字符串或 n_Targets 个字符串，从而允许异构场景：
%   "uav"             - 参数化无人机刚体+旋翼散射点；
%   "legacyDisk"      - 旧 nPoint/r/theta0 圆周四小球或多小球模型；
%   "rigidPointCloud" - 调用方给定的三维刚体点云；
%   "point"           - 单散射点目标；
%   "projectedRange"  - 已有逐脉冲 tgt_posr，仅支持 stop-and-go。
if ~exist('targetModelTypes','var') || isempty(targetModelTypes)
    targetModelTypes = repmat("uav", n_Targets, 1);
end

% 可选的参数化无人机散射点模板：
%   cfg_quad_small          - Phantom 4 类小型四旋翼，默认 6000 RPM，适合小型消费级无人机目标。
%   cfg_quad_heavy          - M350 类重型四旋翼，默认 3600 RPM，适合较大机身和桨叶的工业级无人机目标。
%   cfg_hex_heavy           - M600 类重型六旋翼，默认 3200 RPM，适合六旋翼载荷无人机目标。
%   cfg_hybrid_vtol_4plus1  - TX25A-like 复合翼 4+1 VTOL 参数模型，不代表精确厂商 CAD；
%       "hover"             - 四个升力旋翼高速旋转，机头拉进旋翼停止，表示悬停模式；
%       "transition"        - 升力旋翼由 2400 RPM 降至 0，拉进旋翼由 0 升至 3000 RPM，表示转换模式；
%       "cruise"            - 升力旋翼停止，机头拉进旋翼以 3000 RPM 旋转，表示固定翼巡航模式。
% 前三种模板默认采用恒定 RPM；目标的 hover/linear/accelerated/yawScan 等刚体运动由 motion 另行定义。
% 复合翼模式选择示例：uavConfig = uavmodel.configs.cfg_hybrid_vtol_4plus1("transition");
uavConfig = uavmodel.configs.cfg_quad_small();
model = uavmodel.createModel(uavConfig);

% legacy 与通用刚体点云的默认示例参数；实际研究中可替换这些模板和复系数。
% Geometry and deterministic complex field reflectivity share one spec.
% Ntip~=0 selects the cone; Ntip==0 selects four body x-z-plane spheres.
legacyDiskSpec = struct('Ntip',1,'radiusM',0.3,'phase0Rad',pi/6, ...
    'axialSeparationM',1.5,'alpha0',[1;0.3;0.3;0.3]);
% Equivalent to the first-generation R0_body_to_world convention after
% mapping its +y LOS onto this script's +x LOS while preserving radar up:
% body-x/body-y/body-z make [135,45,90] deg with LOS, respectively.
legacyDiskInitialAttitudeYprRad = [-3*pi/4,0,0];
rigidPointCloud_B = [0.20,0,0; 0,0.20,0; -0.20,0,0; 0,-0.20,0];
rigidPointCloudAlpha0 = legacyDiskSpec.alpha0;
singlePointAlpha0 = 1;
projectedRangeOffsetM = zeros(cfg_RVpath.Na, 1);
projectedRangePulseTimeSec = (0:cfg_RVpath.Na-1).' .* cfg_RVpath.T;

% 仅当场景使用 UAV 且允许绘图时绘制 t=0 s 机体系散射点几何。
if any(strcmpi(targetModelTypes, "uav")) && UAV_PLOT_GEOMETRY_ENABLED
    uavmodel.plotGeometry(model, 0, struct( ...
        'ShowLabels', true, 'ShowFrame', true, 'ShowRotorNormals', true));
end

%% 目标运动模型与雷达状态定义
% 当前令惯性系 I 与雷达坐标系 R 的坐标轴平行，雷达初始位置为原点。
% LOS 从雷达指向目标，默认沿雷达系 +x；距离率为正表示目标远离雷达。
losHat_R = [1, 0, 0];
losHat_R = losHat_R ./ norm(losHat_R);

% 雷达运动设置入口。方向为 [0,0,0] 且速度为 0 时表示雷达静止；
% 若设置运动雷达，请同时给出非零方向和非负速度，当前假设雷达做匀速直线运动。
radarPosition0_I = [0, 0, 0];
radarReferenceTimeSec = 0;
radarMotionDirection_R = [0, 0, 0];
radarSpeedMps = 0;
if norm(radarMotionDirection_R) == 0
    if radarSpeedMps ~= 0
        error('非零雷达速度必须配合非零 radarMotionDirection_R。');
    end
    radarMotionDirectionHat_R = [0, 0, 0];
else
    radarMotionDirectionHat_R = radarMotionDirection_R ./ norm(radarMotionDirection_R);
end
radarVelocity_I = radarSpeedMps .* radarMotionDirectionHat_R;

radarState = struct();
radarState.position_I = @(tQuery) radarPosition0_I + ...
    (double(tQuery(:)) - radarReferenceTimeSec) .* radarVelocity_I;
radarState.velocity_I = radarVelocity_I;
radarState.carrierFrequencyHz = cfg_RVpath.fc;
radarState.wavelengthM = c ./ cfg_RVpath.fc;
radarState.dopplerSign = -1;  % fD=-2*rangeRate/lambda，接近雷达时多普勒为正。
radarState.polarization = 'VV';

% 通过两个角度设置目标正速度/正加速度所对应的公共运动方向：
%   targetMotionLosAngleDeg       - 与 +LOS（雷达系 +x）的夹角；
%       0 deg 表示沿 +x 远离，90 deg 表示横向运动，180 deg 表示沿 -x 接近。
%   targetMotionAzimuthDeg        - 横向分量绕 LOS 的方位角；
%       0 deg 指向雷达系 +y，90 deg 指向雷达系 +z。
targetMotionLosAngleDeg = 45;
targetMotionAzimuthDeg = 0;
targetMotionLosAngleRad = targetMotionLosAngleDeg .* pi ./ 180;
targetMotionAzimuthRad = targetMotionAzimuthDeg .* pi ./ 180;

losSideHat_R = [0, 1, 0];
losUpHat_R = [0, 0, 1];
targetMotionDirectionHat_R = cos(targetMotionLosAngleRad) .* losHat_R + ...
    sin(targetMotionLosAngleRad) .* (...
    cos(targetMotionAzimuthRad) .* losSideHat_R + ...
    sin(targetMotionAzimuthRad) .* losUpHat_R);
targetMotionDirectionHat_R = targetMotionDirectionHat_R ./ ...
    norm(targetMotionDirectionHat_R);

% 换算后的单位运动矢量及其与雷达系 x/y/z 三轴的夹角，顺序为 [x,y,z]。
targetMotionAngleXYZRad = acos(max(-1, min(1, targetMotionDirectionHat_R)));
targetMotionAngleXYZDeg = targetMotionAngleXYZRad .* 180 ./ pi;
fprintf('目标运动单位方向矢量（雷达系）：[%.6f, %.6f, %.6f]\n', ...
    targetMotionDirectionHat_R(1), targetMotionDirectionHat_R(2), ...
    targetMotionDirectionHat_R(3));
fprintf('目标运动方向与雷达系 x/y/z 轴夹角：[%.3f, %.3f, %.3f] deg\n', ...
    targetMotionAngleXYZDeg(1), targetMotionAngleXYZDeg(2), ...
    targetMotionAngleXYZDeg(3));

% 将以 Mach/g 表示的理论真值转换到 SI，并沿公共运动方向构造三维状态。
targetSpeed0Mps = v_Mach(:) .* Mach;
targetAccelerationMps2 = a_g(:) .* g;
targetPosition0_I = radarPosition0_I + R0(:) .* losHat_R;
initialLosProjection = dot(targetMotionDirectionHat_R,losHat_R);
v_radial0 = targetSpeed0Mps .* initialLosProjection - ...
    radarVelocity_I * losHat_R.';
a_radial = targetAccelerationMps2 .* initialLosProjection;

% The monostatic radar remains fixed.  Preserve only the t=0 LOS
% projections of target velocity and acceleration in the target-centre
% motion; do not introduce a transverse centre trajectory.
targetVelocity0_I = v_radial0 .* losHat_R;
targetAcceleration0_I = a_radial .* losHat_R;

% Echo truth contains radial acceleration.  The downstream range predictor
% intentionally remains linear, R + v_r*t, and leaves this residual to RCMC.
motionPulseTimeSec = (0:cfg_ISARpath.Na-1).' .* cfg_ISARpath.Ta;
motionTrueCenterRangeM = R0(:).' + ...
    motionPulseTimeSec .* v_radial0(:).' + ...
    0.5 .* motionPulseTimeSec.^2 .* a_radial(:).';
motionLinearPredictedRangeM = R0(:).' + ...
    motionPulseTimeSec .* v_radial0(:).';
initialLosMotionDiagnostics = struct( ...
    'model',"fixedRadar_initialLos_radialAccelerationTruth_linearPrediction", ...
    'initialLosHat',losHat_R, ...
    'projection',initialLosProjection, ...
    'radialVelocityMps',v_radial0, ...
    'radialAccelerationMps2',a_radial, ...
    'pulseTimeSec',motionPulseTimeSec, ...
    'trueCenterRangeM',motionTrueCenterRangeM, ...
    'linearPredictedRangeM',motionLinearPredictedRangeM, ...
    'accelerationResidualM',motionTrueCenterRangeM-motionLinearPredictedRangeM, ...
    'maximumTransverseVelocityMps',max(vecnorm( ...
        targetVelocity0_I-(targetVelocity0_I*losHat_R.').*losHat_R,2,2)), ...
    'maximumTransverseAccelerationMps2',max(vecnorm( ...
        targetAcceleration0_I-(targetAcceleration0_I*losHat_R.').*losHat_R,2,2)));

% 每个目标使用同一个散射点结构模板，但具有独立的质心运动真值。
motion = repmat(struct('type', 'accelerated', 'referenceTimeSec', 0, ...
    'p0_I', [0, 0, 0], 'v0_I', [0, 0, 0], 'a0_I', [0, 0, 0], ...
    'attitude0YprRad', [0, 0, 0]), n_Targets, 1);
for targetIndex = 1:n_Targets
    motion(targetIndex).p0_I = targetPosition0_I(targetIndex, :);
    motion(targetIndex).v0_I = targetVelocity0_I(targetIndex, :);
    motion(targetIndex).a0_I = targetAcceleration0_I(targetIndex, :);
end

%% 将几何模板和运动状态组合成统一散射体 target 数组
targetModelTypes = string(targetModelTypes(:));
if isscalar(targetModelTypes)
    targetModelTypes = repmat(targetModelTypes, n_Targets, 1);
elseif numel(targetModelTypes) ~= n_Targets
    error('targetModelTypes 必须是标量或包含 n_Targets 个元素。');
end
targetCells = cell(n_Targets, 1);
for targetIndex = 1:n_Targets
    targetName = sprintf('target%02d_%s', targetIndex, targetModelTypes(targetIndex));
    switch lower(targetModelTypes(targetIndex))
        case "uav"
            targetCells{targetIndex} = scatterermodel.createUavTarget( ...
                model, motion(targetIndex), struct('Name',targetName, ...
                'Grouping','bodyRotor'));
        case "legacydisk"
            legacyMotion = motion(targetIndex);
            legacyMotion.attitude0YprRad = legacyDiskInitialAttitudeYprRad;
            targetCells{targetIndex} = scatterermodel.createLegacyDiskTarget( ...
                legacyDiskSpec,legacyMotion,struct('Name',targetName));
        case "rigidpointcloud"
            targetCells{targetIndex} = scatterermodel.createRigidPointCloudTarget( ...
                rigidPointCloud_B,rigidPointCloudAlpha0,motion(targetIndex), ...
                struct('Name',targetName));
        case "point"
            targetCells{targetIndex} = scatterermodel.createPointTarget( ...
                motion(targetIndex),singlePointAlpha0,struct('Name',targetName));
        case "projectedrange"
            targetCells{targetIndex} = scatterermodel.createProjectedRangeTarget( ...
                projectedRangeOffsetM,singlePointAlpha0,motion(targetIndex), ...
                projectedRangePulseTimeSec,struct('Name',targetName));
        otherwise
            error('不支持的 targetModelTypes(%d)="%s"。', ...
                targetIndex,targetModelTypes(targetIndex));
    end
end
targets = vertcat(targetCells{:});
%% 运动参数测量
if run_RVpath
    multiuav_rvpath;
end
%% 运动参数传递
if run_RVpath
    R_isar = rvMeasurement.calibratedRangeM;
    v_isar = rvMeasurement.calibratedVelocityMps;
else
    R_isar = R0;
    v_isar = v_radial0;
end

%% ISAR 成像
if run_ISARpath
    multiuav_isarpath;
end
