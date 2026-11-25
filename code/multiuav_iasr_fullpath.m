% 测速测距 + ISAR 成像全链路执行脚本（参数与开关已独立至 config/rv_isar_default_config.m）
% 调用方式：
%   1. 推荐：通过 main.m 或外部脚本构造 simCfg = rv_isar_default_config()，
%      设置 simCfg.switches.runMode = "rv" | "isar" | "both" 后再运行本脚本；
%   2. 兼容：直接运行本脚本将自动加载默认配置，并兼容旧工作区覆盖变量。
if isempty(which('radar.echo_generate_target_lfm')) || isempty(which('rv_isar_default_config'))
    run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'startup_project.m'));
end

%% 加载并合并外部配置
if exist('simCfg', 'var') && ~isempty(simCfg)
    simCfg = rv_isar_default_config(simCfg);
else
    simCfg = rv_isar_default_config();
end

%% 常量导出
c     = simCfg.constants.c;
Z     = simCfg.constants.Z;
q     = simCfg.constants.q;
lc    = simCfg.constants.lc;
nu    = c / lc;
h     = simCfg.constants.h;
kB    = simCfg.constants.kB;
T_exp = simCfg.constants.T_exp;
Mach  = simCfg.constants.Mach;
g     = simCfg.constants.g;
win_G = simCfg.constants.win_G;

%% 运行模式与开关解析（支持外部 simCfg.switches.runMode 及旧工作区布尔开关）
switch lower(string(simCfg.switches.runMode))
    case "rv"
        defaultRunRV   = true;
        defaultRunISAR = false;
    case "isar"
        defaultRunRV   = false;
        defaultRunISAR = true;
    case "both"
        defaultRunRV   = true;
        defaultRunISAR = true;
    otherwise
        error('RV_ISAR_Path:InvalidRunMode', ...
            'simCfg.switches.runMode 必须为 "rv"、"isar" 或 "both"，当前为 "%s"。', ...
            string(simCfg.switches.runMode));
end

if ~exist('run_RVpath', 'var') || isempty(run_RVpath)
    run_RVpath = defaultRunRV;
else
    run_RVpath = logical(run_RVpath);
end
if ~exist('run_ISARpath', 'var') || isempty(run_ISARpath)
    run_ISARpath = defaultRunISAR;
else
    run_ISARpath = logical(run_ISARpath);
end
if ~run_RVpath && ~run_ISARpath
    warning('RV_ISAR_Path:NoPathSelected', ...
        'run_RVpath 与 run_ISARpath 均为 false，本次未执行任何链路。');
end

if ~exist('UAV_PLOT_GEOMETRY_ENABLED', 'var') || isempty(UAV_PLOT_GEOMETRY_ENABLED)
    UAV_PLOT_GEOMETRY_ENABLED = logical(simCfg.switches.uavPlotGeometryEnabled);
end
if ~exist('rvShowEchoTimeFrequency', 'var') || isempty(rvShowEchoTimeFrequency)
    rvShowEchoTimeFrequency = logical(simCfg.switches.rvShowEchoTimeFrequency);
end
if ~exist('rvTimeFrequencyPulseIndex', 'var') || isempty(rvTimeFrequencyPulseIndex)
    rvTimeFrequencyPulseIndex = simCfg.switches.rvTimeFrequencyPulseIndex;
end
if ~exist('rvLfmSpectrumMode', 'var') || isempty(rvLfmSpectrumMode)
    rvLfmSpectrumMode = simCfg.switches.rvLfmSpectrumMode;
end
if ~exist('rvCfarPlotEnabled', 'var') || isempty(rvCfarPlotEnabled)
    rvCfarPlotEnabled = logical(simCfg.switches.rvCfarPlotEnabled);
end
if ~exist('rvCfarPlotPulseIndex', 'var') || isempty(rvCfarPlotPulseIndex)
    rvCfarPlotPulseIndex = simCfg.switches.rvCfarPlotPulseIndex;
end
if ~exist('ISAR_PLOT_ENABLED', 'var') || isempty(ISAR_PLOT_ENABLED)
    ISAR_PLOT_ENABLED = logical(simCfg.switches.isarPlotEnabled);
end
if ~exist('ISAR_PLOT_DYNAMIC_RANGE_DB', 'var') || isempty(ISAR_PLOT_DYNAMIC_RANGE_DB)
    ISAR_PLOT_DYNAMIC_RANGE_DB = simCfg.switches.isarPlotDynamicRangeDb;
end
if ~exist('ISAR_QUALITY_ENABLED', 'var') || isempty(ISAR_QUALITY_ENABLED)
    ISAR_QUALITY_ENABLED = logical(simCfg.switches.isarQualityEnabled);
end
if ~exist('ISAR_QUALITY_EXPECTED_PEAK_COUNT', 'var') || isempty(ISAR_QUALITY_EXPECTED_PEAK_COUNT)
    ISAR_QUALITY_EXPECTED_PEAK_COUNT = simCfg.switches.isarQualityExpectedPeakCount;
end
if ~exist('ISAR_QUALITY_3DB_METHOD', 'var') || isempty(ISAR_QUALITY_3DB_METHOD)
    ISAR_QUALITY_3DB_METHOD = simCfg.switches.isarQuality3dBMethod;
end
if ~exist('ISAR_QUALITY_PRINT_ENABLED', 'var') || isempty(ISAR_QUALITY_PRINT_ENABLED)
    ISAR_QUALITY_PRINT_ENABLED = logical(simCfg.switches.isarQualityPrintEnabled);
end

%% 导出并校验 RV / ISAR 波形参数
cfg_RVpath = simCfg.rv;
if cfg_RVpath.Ta < cfg_RVpath.T
    error("In RVpath, Ta should be at least the same with T or larger\n");
end

cfg_ISARpath = simCfg.isar;
if exist('FULLPATH_OVERRIDES', 'var') && isstruct(FULLPATH_OVERRIDES)
    if isfield(FULLPATH_OVERRIDES, 'isarNa') && ~isempty(FULLPATH_OVERRIDES.isarNa)
        cfg_ISARpath.Na = double(FULLPATH_OVERRIDES.isarNa);
    end
    if isfield(FULLPATH_OVERRIDES, 'isarEchoMode') && ~isempty(FULLPATH_OVERRIDES.isarEchoMode)
        cfg_ISARpath.echoMode = string(FULLPATH_OVERRIDES.isarEchoMode);
    end
end
cfg_ISARpath.propagationSpeed = c;
cfg_ISARpath.B      = abs(cfg_ISARpath.fH - cfg_ISARpath.fL);
cfg_ISARpath.kr     = (cfg_ISARpath.fH - cfg_ISARpath.fL) / cfg_ISARpath.Tr;
cfg_ISARpath.fc     = (cfg_ISARpath.fH + cfg_ISARpath.fL) / 2;
cfg_ISARpath.lambda = cfg_ISARpath.propagationSpeed / cfg_ISARpath.fc;

if exist('ISAR_MOTION_OPTIONS', 'var') && ~isempty(ISAR_MOTION_OPTIONS)
    motionOptionNames = {'Mode','SpinFrequencyHz', ...
        'PrecessionFrequencyHz','PrecessionHalfAngleDeg'};
    motionConfigNames = {'bodyRotationMode','spinFrequencyHz', ...
        'precessionFrequencyHz','precessionHalfAngleDeg'};
    if ~(isstruct(ISAR_MOTION_OPTIONS) && isscalar(ISAR_MOTION_OPTIONS)) || ...
            any(~ismember(fieldnames(ISAR_MOTION_OPTIONS), motionOptionNames))
        error('RadarProject:ISARMotion:InvalidOptions', ...
            'ISAR_MOTION_OPTIONS contains an unsupported field.');
    end
    for motionOptionIndex = 1:numel(motionOptionNames)
        motionOptionName = motionOptionNames{motionOptionIndex};
        if isfield(ISAR_MOTION_OPTIONS, motionOptionName)
            cfg_ISARpath.(motionConfigNames{motionOptionIndex}) = ...
                ISAR_MOTION_OPTIONS.(motionOptionName);
        end
    end
end
cfg_ISARpath.omegaBody_B = [0, 2*pi*cfg_ISARpath.spinFrequencyHz, 0];

if exist('ISAR_POST_OPTIONS', 'var') && ~isempty(ISAR_POST_OPTIONS)
    cfg_ISARpath.postprocess = ISAR_POST_OPTIONS;
elseif isfield(cfg_ISARpath, 'postprocess') && isempty(fieldnames(cfg_ISARpath.postprocess))
    cfg_ISARpath = rmfield(cfg_ISARpath, 'postprocess');
end
if exist('ISAR_PRIOR_ACCELERATION_MPS2', 'var') && ~isempty(ISAR_PRIOR_ACCELERATION_MPS2)
    cfg_ISARpath.priorAccelerationMps2 = ISAR_PRIOR_ACCELERATION_MPS2;
elseif isfield(cfg_ISARpath, 'priorAccelerationMps2') && isempty(cfg_ISARpath.priorAccelerationMps2)
    cfg_ISARpath = rmfield(cfg_ISARpath, 'priorAccelerationMps2');
end
if cfg_ISARpath.Ta < cfg_ISARpath.T
    error("In ISARpath, Ta should be at least the same with T or larger\n");
end

%% 其余仿真参数定义
R0     = simCfg.targets.R0(:);
v_Mach = simCfg.targets.v_Mach(:);
a_g    = simCfg.targets.a_g(:);
if ~(numel(R0) == numel(v_Mach) && numel(R0) == numel(a_g))
    error('RV_ISAR_Path:TargetVectorLengthMismatch', ...
        'simCfg.targets 的 R0(%d)、v_Mach(%d)、a_g(%d) 长度必须一致。', ...
        numel(R0), numel(v_Mach), numel(a_g));
end

idxKeep = simCfg.targets.targetIndices;
if exist('FULLPATH_OVERRIDES', 'var') && isstruct(FULLPATH_OVERRIDES) && ...
        isfield(FULLPATH_OVERRIDES, 'targetIndices') && ~isempty(FULLPATH_OVERRIDES.targetIndices)
    idxKeep = FULLPATH_OVERRIDES.targetIndices;
end
if ~isempty(idxKeep)
    idxKeep = idxKeep(:);
    R0     = R0(idxKeep);
    v_Mach = v_Mach(idxKeep);
    a_g    = a_g(idxKeep);
end
n_Targets = numel(R0);

%% 目标散射体模型加载
if ~exist('targetModelTypes', 'var') || isempty(targetModelTypes)
    targetModelTypes = simCfg.scatterer.modelTypes;
end
targetModelTypes = string(targetModelTypes(:));
if isscalar(targetModelTypes)
    targetModelTypes = repmat(targetModelTypes, n_Targets, 1);
elseif numel(targetModelTypes) ~= n_Targets
    error('targetModelTypes 必须是标量或包含 n_Targets(%d) 个元素。', n_Targets);
end

if any(strcmpi(targetModelTypes, "uav"))
    uavConfig = simCfg.scatterer.uavConfigFcn();
    model = uavmodel.createModel(uavConfig);
    if UAV_PLOT_GEOMETRY_ENABLED
        uavmodel.plotGeometry(model, 0, struct( ...
            'ShowLabels', true, 'ShowFrame', true, 'ShowRotorNormals', true));
    end
end

legacyDiskSpec                  = simCfg.scatterer.legacyDiskSpec;
legacyDiskInitialAttitudeYprRad = simCfg.scatterer.legacyDiskInitialAttitudeYprRad;
rigidPointCloud_B               = simCfg.scatterer.rigidPointCloud_B;
rigidPointCloudAlpha0           = simCfg.scatterer.rigidPointCloudAlpha0;
singlePointAlpha0               = simCfg.scatterer.singlePointAlpha0;
projectedRangeOffsetM           = zeros(cfg_RVpath.Na, 1);
projectedRangePulseTimeSec      = (0:cfg_RVpath.Na-1).' .* cfg_RVpath.T;

%% 目标运动模型与雷达状态定义
losHat_R = simCfg.targets.losHat_R(:).';
losHat_R = losHat_R ./ norm(losHat_R);

radarPosition0_I      = simCfg.targets.radarPosition0_I(:).';
radarReferenceTimeSec = simCfg.targets.radarReferenceTimeSec;
radarMotionDirection_R = simCfg.targets.radarMotionDirection_R(:).';
radarSpeedMps         = simCfg.targets.radarSpeedMps;
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
radarState.velocity_I         = radarVelocity_I;
radarState.carrierFrequencyHz = cfg_RVpath.fc;
radarState.wavelengthM        = c ./ cfg_RVpath.fc;
radarState.dopplerSign        = simCfg.targets.dopplerSign;
radarState.polarization       = simCfg.targets.polarization;

targetMotionLosAngleDeg = simCfg.targets.motionLosAngleDeg;
targetMotionAzimuthDeg  = simCfg.targets.motionAzimuthDeg;
targetMotionLosAngleRad = targetMotionLosAngleDeg .* pi ./ 180;
targetMotionAzimuthRad  = targetMotionAzimuthDeg .* pi ./ 180;

losSideHat_R = [0, 1, 0];
losUpHat_R   = [0, 0, 1];
targetMotionDirectionHat_R = cos(targetMotionLosAngleRad) .* losHat_R + ...
    sin(targetMotionLosAngleRad) .* ( ...
    cos(targetMotionAzimuthRad) .* losSideHat_R + ...
    sin(targetMotionAzimuthRad) .* losUpHat_R);
targetMotionDirectionHat_R = targetMotionDirectionHat_R ./ norm(targetMotionDirectionHat_R);

targetMotionAngleXYZRad = acos(max(-1, min(1, targetMotionDirectionHat_R)));
targetMotionAngleXYZDeg = targetMotionAngleXYZRad .* 180 ./ pi;
fprintf('目标运动单位方向矢量（雷达系）：[%.6f, %.6f, %.6f]\n', ...
    targetMotionDirectionHat_R(1), targetMotionDirectionHat_R(2), ...
    targetMotionDirectionHat_R(3));
fprintf('目标运动方向与雷达系 x/y/z 轴夹角：[%.3f, %.3f, %.3f] deg\n', ...
    targetMotionAngleXYZDeg(1), targetMotionAngleXYZDeg(2), ...
    targetMotionAngleXYZDeg(3));

targetSpeed0Mps        = v_Mach(:) .* Mach;
targetAccelerationMps2 = a_g(:) .* g;
targetPosition0_I      = radarPosition0_I + R0(:) .* losHat_R;
initialLosProjection   = dot(targetMotionDirectionHat_R, losHat_R);
v_radial0              = targetSpeed0Mps .* initialLosProjection - radarVelocity_I * losHat_R.';
a_radial               = targetAccelerationMps2 .* initialLosProjection;

targetVelocity0_I     = v_radial0 .* losHat_R;
targetAcceleration0_I = a_radial .* losHat_R;

motionPulseTimeSec          = (0:cfg_ISARpath.Na-1).' .* cfg_ISARpath.Ta;
motionTrueCenterRangeM      = R0(:).' + motionPulseTimeSec .* v_radial0(:).' + ...
    0.5 .* motionPulseTimeSec.^2 .* a_radial(:).';
motionLinearPredictedRangeM = R0(:).' + motionPulseTimeSec .* v_radial0(:).';
initialLosMotionDiagnostics = struct( ...
    'model',"fixedRadar_initialLos_radialAccelerationTruth_linearPrediction", ...
    'initialLosHat',losHat_R, ...
    'projection',initialLosProjection, ...
    'radialVelocityMps',v_radial0, ...
    'radialAccelerationMps2',a_radial, ...
    'pulseTimeSec',motionPulseTimeSec, ...
    'trueCenterRangeM',motionTrueCenterRangeM, ...
    'linearPredictedRangeM',motionLinearPredictedRangeM, ...
    'accelerationResidualM',motionTrueCenterRangeM - motionLinearPredictedRangeM, ...
    'maximumTransverseVelocityMps',max(vecnorm( ...
        targetVelocity0_I - (targetVelocity0_I*losHat_R.').*losHat_R, 2, 2)), ...
    'maximumTransverseAccelerationMps2',max(vecnorm( ...
        targetAcceleration0_I - (targetAcceleration0_I*losHat_R.').*losHat_R, 2, 2)));

motion = repmat(struct('type', 'accelerated', 'referenceTimeSec', 0, ...
    'p0_I', [0, 0, 0], 'v0_I', [0, 0, 0], 'a0_I', [0, 0, 0], ...
    'attitude0YprRad', [0, 0, 0]), n_Targets, 1);
for targetIndex = 1:n_Targets
    motion(targetIndex).p0_I = targetPosition0_I(targetIndex, :);
    motion(targetIndex).v0_I = targetVelocity0_I(targetIndex, :);
    motion(targetIndex).a0_I = targetAcceleration0_I(targetIndex, :);
end

%% 将几何模板和运动状态组合成统一散射体 target 数组
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
                legacyDiskSpec, legacyMotion, struct('Name',targetName));
        case "rigidpointcloud"
            targetCells{targetIndex} = scatterermodel.createRigidPointCloudTarget( ...
                rigidPointCloud_B, rigidPointCloudAlpha0, motion(targetIndex), ...
                struct('Name',targetName));
        case "point"
            targetCells{targetIndex} = scatterermodel.createPointTarget( ...
                motion(targetIndex), singlePointAlpha0, struct('Name',targetName));
        case "projectedrange"
            targetCells{targetIndex} = scatterermodel.createProjectedRangeTarget( ...
                projectedRangeOffsetM, singlePointAlpha0, motion(targetIndex), ...
                projectedRangePulseTimeSec, struct('Name',targetName));
        otherwise
            error('不支持的 targetModelTypes(%d)="%s"。', ...
                targetIndex, targetModelTypes(targetIndex));
    end
end
targets = vertcat(targetCells{:});

%% 运动参数测量（RV_path）
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

%% ISAR 成像（ISAR_path）
if run_ISARpath
    multiuav_isarpath;
end
