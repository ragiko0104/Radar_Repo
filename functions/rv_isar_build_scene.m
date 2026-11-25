function scene = rv_isar_build_scene(simCfg)
%RV_ISAR_BUILD_SCENE 根据仿真配置构建统一的 RV + ISAR 纯计算场景结构体。
%
% 功能说明：
%   1. 将外部传入的配置结构体与 rv_isar_default_config 默认模板递归合并，并校验关键约束；
%   2. 导出物理常量、核心计算开关（runMode / rvLfmSpectrumMode）与 RV/ISAR 派生波形参数；
%   3. 解析并展开 ISAR 逐目标本振拍频起始位置 frequencyBiasHz（正=本振在目标左侧，
%      负=本振在目标右侧，默认循环 [20e6, 45e6] Hz）以及完整 ISAR 后处理/RCMC 参数
%      isarSettings（含可配置的 rcmcSearchRadiusBins 最大搜索位移半径）；
%   4. 按 targetIndices 切片目标列表，计算雷达状态、目标初始 LOS 投影与真值运动学；
%   5. 按 scatterer.modelTypes 按需实例化散射体模型（仅构造启用的模型），并组装
%      统一散射体目标数组 targets。
%   注：本函数纯计算、无任何图形弹窗副作用；几何可视化请使用 plot_scatterer_geometry。
%
% 输入参数：
%   simCfg - (可选) 1×1 struct，外部仿真配置或局部覆盖结构体（结构同
%            rv_isar_default_config 返回值）。省略或为空时使用默认配置。
%            主要字段：
%              .constants  1×1 struct，物理常量（c/Z/q/lc/h/kB/T_exp/Mach/g/win_G）
%              .switches   1×1 struct，runMode（"rv"|"isar"|"both"）与 rvLfmSpectrumMode
%              .targets    1×1 struct，目标初始距离/速度/加速度、索引切片及雷达状态
%              .scatterer  1×1 struct，散射体模型选择 modelTypes 与各模型子配置
%              .rv         1×1 struct，RV 链路六段 LFM 波形与一致性门限参数
%              .isar       1×1 struct，ISAR 链路宽带 LFM 波形、微动姿态、
%                          frequencyBiasHz（正=左侧/负=右侧）及 postprocess（含
%                          rcmcSearchRadiusBins 等全部后处理配置）
%
% 输出参数：
%   scene  - 1×1 struct，包含全链路所需的完整场景描述，字段如下：
%     .cfg            1×1 struct      合并并更新 scatterer 后的完整配置结构体
%     .constants      1×1 struct      物理常量（含派生光频 nu = c/lc，均为 1×1 double）
%     .switches       1×1 struct      归一化计算开关（runMode/runRV/runISAR/rvLfmSpectrumMode）
%     .cfgRV          1×1 struct      校验后的 RV 链路波形与一致性配置
%     .cfgISAR        1×1 struct      校验并补齐派生量（propagationSpeed/B/kr/fc/lambda/
%                                     omegaBody_B/frequencyBiasHz N×1/timeBiasSec N×1/
%                                     postprocess）的 ISAR 链路配置
%     .isarSettings   1×1 struct      解析并校验完毕的 ISAR 后处理与 RCMC 完整参数
%     .radarState     1×1 struct      雷达状态（position_I 句柄，输入 Nt×1 返回 Nt×3；
%                                     velocity_I 1×3；carrierFrequencyHz/wavelengthM/
%                                     dopplerSign 均为 1×1 double；polarization 字符向量）
%     .truth          1×1 struct      目标真值运动学与诊断信息（R0/v_radial0/a_radial 等 N×1）
%     .motion         N×1 struct      逐目标刚体平动描述数组
%     .scatterer      1×1 struct      散射体容器（含 activeModelTypes N×1、按需实例化的模型及 targets）
%     .targets        N×1 struct      统一散射体目标数组，可直接传入回波生成函数
%     .nTargets       1×1 double      有效目标个数 N

if nargin < 1 || isempty(simCfg)
    simCfg = rv_isar_default_config();
else
    simCfg = rv_isar_default_config(simCfg);
end

%% 1. 物理常量
constants = simCfg.constants;
constants.nu = constants.c / constants.lc;

%% 2. 核心计算开关解析
switches = simCfg.switches;
runModeNormalized = lower(string(switches.runMode));
switch runModeNormalized
    case "rv"
        switches.runRV   = true;
        switches.runISAR = false;
    case "isar"
        switches.runRV   = false;
        switches.runISAR = true;
    case "both"
        switches.runRV   = true;
        switches.runISAR = true;
    otherwise
        error('RV_ISAR_Path:InvalidRunMode', ...
            'simCfg.switches.runMode 必须为 "rv"、"isar" 或 "both"，当前为 "%s"。', ...
            string(switches.runMode));
end
switches.runMode = runModeNormalized;

%% 3. 目标参数切片与校验
R0     = double(simCfg.targets.R0(:));
v_Mach = double(simCfg.targets.v_Mach(:));
a_g    = double(simCfg.targets.a_g(:));
if ~(numel(R0) == numel(v_Mach) && numel(R0) == numel(a_g))
    error('RV_ISAR_Path:TargetVectorLengthMismatch', ...
        'simCfg.targets 的 R0(%d)、v_Mach(%d)、a_g(%d) 长度必须一致。', ...
        numel(R0), numel(v_Mach), numel(a_g));
end

if ~isempty(simCfg.targets.targetIndices)
    targetIndices = double(simCfg.targets.targetIndices(:));
    R0     = R0(targetIndices);
    v_Mach = v_Mach(targetIndices);
    a_g    = a_g(targetIndices);
else
    targetIndices = (1:numel(R0)).';
end
nTargets = numel(R0);

%% 4. RV / ISAR 波形参数、逐目标本振拍频偏置与后处理/RCMC 参数解析
cfgRV = simCfg.rv;
if cfgRV.Ta < cfgRV.T
    error('RV_ISAR_Path:InvalidRVTiming', ...
        'In RVpath, Ta (%.6g) must be >= T (%.6g).', cfgRV.Ta, cfgRV.T);
end

cfgISAR = simCfg.isar;
if cfgISAR.Ta < cfgISAR.T
    error('RV_ISAR_Path:InvalidISARTiming', ...
        'In ISARpath, Ta (%.6g) must be >= T (%.6g).', cfgISAR.Ta, cfgISAR.T);
end
cfgISAR.propagationSpeed = constants.c;
cfgISAR.B                = abs(cfgISAR.fH - cfgISAR.fL);
cfgISAR.kr               = (cfgISAR.fH - cfgISAR.fL) / cfgISAR.Tr;
cfgISAR.fc               = (cfgISAR.fH + cfgISAR.fL) / 2;
cfgISAR.lambda           = cfgISAR.propagationSpeed / cfgISAR.fc;
if isfield(cfgISAR, 'omegaBody_B') && ~isempty(cfgISAR.omegaBody_B)
    omegaRaw = double(cfgISAR.omegaBody_B(:).');
    if ~(isreal(omegaRaw) && numel(omegaRaw) == 3 && all(isfinite(omegaRaw)))
        error('RV_ISAR_Path:InvalidOmegaBody', ...
            'simCfg.isar.omegaBody_B 必须为包含 3 个有限实数的向量 [wx, wy, wz] (rad/s)。');
    end
    cfgISAR.omegaBody_B = omegaRaw;
else
    cfgISAR.omegaBody_B = [0, 2*pi*cfgISAR.spinFrequencyHz, 0];
end

% 解析每个目标本振相对目标的拍频起始位置（正=左侧/本振提前，负=右侧/本振滞后）
if isfield(cfgISAR, 'frequencyBiasHz') && ~isempty(cfgISAR.frequencyBiasHz)
    rawFreqBias = double(cfgISAR.frequencyBiasHz(:));
elseif isfield(cfgISAR, 'frequencyBiasPatternHz') && ~isempty(cfgISAR.frequencyBiasPatternHz)
    rawFreqBias = double(cfgISAR.frequencyBiasPatternHz(:));
else
    rawFreqBias = [20e6; 45e6];
end
if any(~isfinite(rawFreqBias)) || ~isreal(rawFreqBias)
    error('RV_ISAR_Path:InvalidFrequencyBias', ...
        'simCfg.isar.frequencyBiasHz 必须为有限实数标量或向量。');
end
if isscalar(rawFreqBias)
    freqBiasPerTarget = repmat(rawFreqBias, nTargets, 1);
elseif numel(rawFreqBias) == nTargets
    freqBiasPerTarget = rawFreqBias;
else
    freqBiasPerTarget = rawFreqBias(mod(0:nTargets-1, numel(rawFreqBias)) + 1);
end
cfgISAR.frequencyBiasHz        = freqBiasPerTarget;
cfgISAR.frequencyBiasPatternHz = freqBiasPerTarget.';
cfgISAR.timeBiasSec            = freqBiasPerTarget ./ cfgISAR.kr;

% 合并并校验 ISAR 后处理及 RCMC 搜索参数
if isfield(cfgISAR, 'postprocess') && isstruct(cfgISAR.postprocess) && ...
        isfield(cfgISAR.postprocess, 'hrrfBins') && ...
        ~isfield(cfgISAR.postprocess, 'rcmcCorrelationPoints')
    cfgISAR.postprocess.rcmcCorrelationPoints = 2 * double(cfgISAR.postprocess.hrrfBins);
end
isarSettings = isarpath_settings(cfgISAR);
validate_isar_settings(isarSettings);
cfgISAR.postprocess = extract_postprocess_fields(isarSettings);

if isfield(cfgISAR, 'priorAccelerationMps2') && isempty(cfgISAR.priorAccelerationMps2)
    cfgISAR = rmfield(cfgISAR, 'priorAccelerationMps2');
end

%% 5. 雷达状态与目标运动真值
losHat_R = double(simCfg.targets.losHat_R(:).');
losHat_R = losHat_R ./ norm(losHat_R);

radarPosition0_I       = double(simCfg.targets.radarPosition0_I(:).');
radarReferenceTimeSec  = double(simCfg.targets.radarReferenceTimeSec);
radarMotionDirection_R = double(simCfg.targets.radarMotionDirection_R(:).');
radarSpeedMps          = double(simCfg.targets.radarSpeedMps);
if norm(radarMotionDirection_R) == 0
    if radarSpeedMps ~= 0
        error('RV_ISAR_Path:InvalidRadarMotion', ...
            '非零雷达速度必须配合非零 radarMotionDirection_R。');
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
radarState.carrierFrequencyHz = cfgRV.fc;
radarState.wavelengthM        = constants.c ./ cfgRV.fc;
radarState.dopplerSign        = simCfg.targets.dopplerSign;
radarState.polarization       = simCfg.targets.polarization;

targetMotionLosAngleRad = double(simCfg.targets.motionLosAngleDeg) .* pi ./ 180;
targetMotionAzimuthRad  = double(simCfg.targets.motionAzimuthDeg)  .* pi ./ 180;
losSideHat_R = [0, 1, 0];
losUpHat_R   = [0, 0, 1];
motionDirectionHat_R = cos(targetMotionLosAngleRad) .* losHat_R + ...
    sin(targetMotionLosAngleRad) .* ( ...
    cos(targetMotionAzimuthRad) .* losSideHat_R + ...
    sin(targetMotionAzimuthRad) .* losUpHat_R);
motionDirectionHat_R = motionDirectionHat_R ./ norm(motionDirectionHat_R);
motionAngleXYZDeg    = acos(max(-1, min(1, motionDirectionHat_R))) .* 180 ./ pi;

targetSpeed0Mps        = v_Mach .* constants.Mach;
targetAccelerationMps2 = a_g .* constants.g;
targetPosition0_I      = radarPosition0_I + R0 .* losHat_R;
initialLosProjection   = dot(motionDirectionHat_R, losHat_R);
v_radial0              = targetSpeed0Mps .* initialLosProjection - radarVelocity_I * losHat_R.';
a_radial               = targetAccelerationMps2 .* initialLosProjection;

targetVelocity0_I     = v_radial0 .* losHat_R;
targetAcceleration0_I = a_radial .* losHat_R;

motionPulseTimeSec          = (0:cfgISAR.Na-1).' .* cfgISAR.Ta;
motionTrueCenterRangeM      = R0.' + motionPulseTimeSec .* v_radial0.' + ...
    0.5 .* motionPulseTimeSec.^2 .* a_radial.';
motionLinearPredictedRangeM = R0.' + motionPulseTimeSec .* v_radial0.';
motionDiagnostics = struct( ...
    'model', "fixedRadar_initialLos_radialAccelerationTruth_linearPrediction", ...
    'initialLosHat', losHat_R, ...
    'projection', initialLosProjection, ...
    'radialVelocityMps', v_radial0, ...
    'radialAccelerationMps2', a_radial, ...
    'pulseTimeSec', motionPulseTimeSec, ...
    'trueCenterRangeM', motionTrueCenterRangeM, ...
    'linearPredictedRangeM', motionLinearPredictedRangeM, ...
    'accelerationResidualM', motionTrueCenterRangeM - motionLinearPredictedRangeM, ...
    'maximumTransverseVelocityMps', max(vecnorm( ...
        targetVelocity0_I - (targetVelocity0_I * losHat_R.') .* losHat_R, 2, 2)), ...
    'maximumTransverseAccelerationMps2', max(vecnorm( ...
        targetAcceleration0_I - (targetAcceleration0_I * losHat_R.') .* losHat_R, 2, 2)));

motion = repmat(struct( ...
    'type', 'accelerated', ...
    'referenceTimeSec', 0, ...
    'p0_I', [0, 0, 0], ...
    'v0_I', [0, 0, 0], ...
    'a0_I', [0, 0, 0], ...
    'attitude0YprRad', [0, 0, 0]), nTargets, 1);
for targetIndex = 1:nTargets
    motion(targetIndex).p0_I = targetPosition0_I(targetIndex, :);
    motion(targetIndex).v0_I = targetVelocity0_I(targetIndex, :);
    motion(targetIndex).a0_I = targetAcceleration0_I(targetIndex, :);
end

%% 6. 散射体模型按需实例化与 targets 组装（统一调用 lib 的 scatterermodel.createTarget）
scatterer = simCfg.scatterer;
if isfield(scatterer, 'preset') && ~isempty(scatterer.preset) && ...
        (~isfield(scatterer, 'modelTypes') || isempty(scatterer.modelTypes))
    scatterer.modelTypes = scatterer.preset;
end

presetCells = normalize_scatterer_presets(scatterer.modelTypes, nTargets);
overrideCells = normalize_scatterer_overrides(scatterer, nTargets);

targetCells      = cell(nTargets, 1);
resolvedConfigs  = cell(nTargets, 1);
activeModelTypes = strings(nTargets, 1);
uavModelCache    = containers.Map('KeyType', 'char', 'ValueType', 'any');
uavConfigCache   = containers.Map('KeyType', 'char', 'ValueType', 'any');
firstUavModel    = [];
firstUavConfig   = [];

for targetIndex = 1:nTargets
    presetSpec   = presetCells{targetIndex};
    itemOverride = overrideCells{targetIndex};

    % projectedRange 需要结合 cfgRV 的脉冲时间轴构造兼容预设
    if is_projected_range_spec(presetSpec)
        projCfg = struct( ...
            'kind',         "projectedRange", ...
            'name',         "projected_range", ...
            'offsetM',      zeros(cfgRV.Na, 1), ...
            'pulseTimeSec', (0:cfgRV.Na-1).' .* cfgRV.T, ...
            'alpha0',       1);
        presetSpec = projCfg;
    end

    createOpts = struct('Name', sprintf('target%02d', targetIndex));
    if ~isempty(itemOverride) && isstruct(itemOverride) && ~isempty(fieldnames(itemOverride))
        createOpts.Overrides = itemOverride;
    end

    % 同款无覆盖 UAV 预设仅编译一次，复用已编译模型以节省开销
    cacheKey = make_preset_cache_key(presetSpec, itemOverride);
    if ~isempty(cacheKey) && isKey(uavModelCache, cacheKey)
        createOpts.CompiledModel = uavModelCache(cacheKey);
    end

    [tgt, compiledUav, resolvedCfg] = scatterermodel.createTarget( ...
        presetSpec, motion(targetIndex), createOpts);
    tgt.name = string(sprintf('target%02d_%s', targetIndex, tgt.kind));

    targetCells{targetIndex}      = tgt;
    resolvedConfigs{targetIndex}  = resolvedCfg;
    activeModelTypes(targetIndex) = string(tgt.kind);

    if ~isempty(compiledUav)
        if ~isempty(cacheKey) && ~isKey(uavModelCache, cacheKey)
            uavModelCache(cacheKey)  = compiledUav;
            uavConfigCache(cacheKey) = resolvedCfg;
        end
        if isempty(firstUavModel)
            firstUavModel  = compiledUav;
            firstUavConfig = resolvedCfg;
        end
    end
end

targets = vertcat(targetCells{:});
scatterer.presets          = presetCells;
scatterer.activeModelTypes = activeModelTypes;
scatterer.configs          = resolvedConfigs;
scatterer.targets          = targets;
if ~isempty(firstUavModel)
    scatterer.uav = struct('config', firstUavConfig, 'model', firstUavModel);
end
simCfg.scatterer        = scatterer;
simCfg.isar.postprocess = cfgISAR.postprocess;

%% 7. 组装统一场景输出
truth = struct( ...
    'nTargets',               nTargets, ...
    'targetIndices',          targetIndices, ...
    'R0',                     R0, ...
    'v_Mach',                 v_Mach, ...
    'a_g',                    a_g, ...
    'targetSpeed0Mps',        targetSpeed0Mps, ...
    'targetAccelerationMps2', targetAccelerationMps2, ...
    'losHat_R',               losHat_R, ...
    'motionDirectionHat_R',   motionDirectionHat_R, ...
    'motionAngleXYZDeg',      motionAngleXYZDeg, ...
    'initialLosProjection',   initialLosProjection, ...
    'v_radial0',              v_radial0, ...
    'a_radial',               a_radial, ...
    'targetPosition0_I',      targetPosition0_I, ...
    'targetVelocity0_I',      targetVelocity0_I, ...
    'targetAcceleration0_I',  targetAcceleration0_I, ...
    'motionDiagnostics',      motionDiagnostics);

scene = struct( ...
    'cfg',          simCfg, ...
    'constants',    constants, ...
    'switches',     switches, ...
    'cfgRV',        cfgRV, ...
    'cfgISAR',      cfgISAR, ...
    'isarSettings', isarSettings, ...
    'radarState',   radarState, ...
    'truth',        truth, ...
    'motion',       motion, ...
    'scatterer',    scatterer, ...
    'targets',      targets, ...
    'nTargets',     nTargets);
end

function validate_isar_settings(p)
if p.rangeFftPoints < p.timeWindowSamples
    error('RadarProject:ISARRangeFftTooShort', ...
        'rangeFftPoints 必须不小于 timeWindowSamples。');
end
if p.hrrfBins < 1 || p.hrrfBins > p.rangeFftPoints
    error('RadarProject:ISARInvalidHrrfBinCount', ...
        'hrrfBins 必须位于 [1, rangeFftPoints]。');
end
if ~ismember(upper(string(p.timeWindowMode)), ["FIXED_FIRST_PULSE", "MOVING_PULSE"])
    error('RadarProject:ISARUnsupportedTimeWindowMode', ...
        '不支持的 timeWindowMode="%s"。', string(p.timeWindowMode));
end
if ~ismember(upper(string(p.rangeMigrationMode)), ...
        ["MODEL", "XCORR", "MODEL_XCORR", "FFT2048_SINGLE", "NONE"])
    error('RadarProject:ISARUnsupportedRangeMigrationMode', ...
        '不支持的 rangeMigrationMode="%s"。', string(p.rangeMigrationMode));
end
if ~(isnumeric(p.rcmcSearchRadiusBins) && isscalar(p.rcmcSearchRadiusBins) && ...
        isfinite(p.rcmcSearchRadiusBins) && p.rcmcSearchRadiusBins >= 0 && ...
        p.rcmcSearchRadiusBins == floor(p.rcmcSearchRadiusBins))
    error('RV_ISAR_Path:InvalidRcmcSearchRadius', ...
        'rcmcSearchRadiusBins 必须为非负整数标量。');
end
end

function post = extract_postprocess_fields(p)
post = struct( ...
    'timeWindowSamples',             p.timeWindowSamples, ...
    'rangeFftPoints',                p.rangeFftPoints, ...
    'hrrfBins',                      p.hrrfBins, ...
    'timeMarginSamples',             p.timeMarginSamples, ...
    'timeWindowMode',                p.timeWindowMode, ...
    'cropCenterOffsetBinsEffective', p.cropCenterOffsetBinsEffective, ...
    'rangeMigrationMode',            p.rangeMigrationMode, ...
    'rangeAlignmentMode',            p.rangeAlignmentMode, ...
    'rcmcCorrelationPoints',         p.rcmcCorrelationPoints, ...
    'rcmcSearchRadiusBins',          p.rcmcSearchRadiusBins, ...
    'phaseCompensationMode',         p.phaseCompensationMode, ...
    'phaseThreshold',                p.phaseThreshold, ...
    'dispersionCompensationEnabled', p.dispersionCompensationEnabled, ...
    'windowType',                    p.windowType, ...
    'slowTimeZeroPaddingFactor',     p.slowTimeZeroPaddingFactor);
end

function presetCells = normalize_scatterer_presets(rawSpec, nTargets)
if isa(rawSpec, 'function_handle') || (isstruct(rawSpec) && isscalar(rawSpec))
    presetCells = repmat({rawSpec}, nTargets, 1);
elseif ischar(rawSpec) && isrow(rawSpec)
    presetCells = repmat({string(rawSpec)}, nTargets, 1);
elseif isstring(rawSpec)
    vals = rawSpec(:);
    if isscalar(vals)
        presetCells = repmat({vals}, nTargets, 1);
    elseif numel(vals) == nTargets
        presetCells = cellstr(vals);
    else
        error('RV_ISAR_Path:ScattererModelTypeLengthMismatch', ...
            'scatterer.modelTypes 必须是标量或包含 nTargets(%d) 个元素。', nTargets);
    end
elseif iscell(rawSpec)
    vals = rawSpec(:);
    if isscalar(vals)
        presetCells = repmat(vals, nTargets, 1);
    elseif numel(vals) == nTargets
        presetCells = vals;
    else
        error('RV_ISAR_Path:ScattererModelTypeLengthMismatch', ...
            'scatterer.modelTypes cell 数组必须是标量或包含 nTargets(%d) 个元素。', nTargets);
    end
else
    error('RV_ISAR_Path:InvalidScattererSpec', ...
        'scatterer.modelTypes 必须为 lib 预设函数句柄、预设名字符串、配置结构体或 N×1 cell 数组。');
end
end

function overrideCells = normalize_scatterer_overrides(scatterer, nTargets)
if ~isfield(scatterer, 'overrides') || isempty(scatterer.overrides)
    overrideCells = repmat({struct()}, nTargets, 1);
elseif isstruct(scatterer.overrides) && isscalar(scatterer.overrides)
    overrideCells = repmat({scatterer.overrides}, nTargets, 1);
elseif iscell(scatterer.overrides)
    vals = scatterer.overrides(:);
    if isscalar(vals)
        overrideCells = repmat(vals, nTargets, 1);
    elseif numel(vals) == nTargets
        overrideCells = vals;
    else
        error('RV_ISAR_Path:ScattererOverrideLengthMismatch', ...
            'scatterer.overrides cell 数组必须是标量或包含 nTargets(%d) 个元素。', nTargets);
    end
else
    error('RV_ISAR_Path:InvalidScattererOverrides', ...
        'scatterer.overrides 必须为空、1×1 struct 或 N×1 cell 数组。');
end
end

function tf = is_projected_range_spec(spec)
tf = ((ischar(spec) && isrow(spec)) || (isstring(spec) && isscalar(spec))) && ...
    strcmpi(strtrim(string(spec)), "projectedrange");
end

function key = make_preset_cache_key(presetSpec, itemOverride)
key = '';
if isstruct(itemOverride) && ~isempty(fieldnames(itemOverride))
    return;
end
if isa(presetSpec, 'function_handle')
    key = ['fh:' func2str(presetSpec)];
elseif (ischar(presetSpec) && isrow(presetSpec)) || (isstring(presetSpec) && isscalar(presetSpec))
    key = ['str:' lower(char(strtrim(string(presetSpec))))];
end
end
