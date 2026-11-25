function isarScene = isarpath_build_scene(scene, priorMotion, options)
%ISARPATH_BUILD_SCENE 构建 ISAR 成像链路（ISAR_path）纯计算场景、先验、微动与坐标轴结构体。
%
% 功能说明：
%   1. 直接读取由 rv_isar_build_scene 统一解析好的 ISAR 后处理/RCMC 参数（isarSettings）
%      及逐目标本振拍频起始位置（cfgISAR.frequencyBiasHz，正=本振在目标左侧，负=右侧）；
%   2. 导出 ISAR 采样率、高采样快时间轴、绝对慢时间轴及 ISAR 载频下的雷达状态；
%   3. 校验并整理逐目标先验运动参数（优先采用传入的 priorMotion，否则回退至 scene.truth），
%      根据符号化拍频偏置计算各目标局部本振起始时刻 referenceStartTimeSec 与单目标 prior；
%   4. 调用 isar_apply_micro_motion、isar_cross_range_geometry 与 scatterermodel.evaluateTarget
%      生成微动目标数组 targetsISAR、方位定标几何与慢时间散射点状态；
%   5. 校验加速度残差是否落在 rcmcSearchRadiusBins 范围内，并构造距离轴与方位定标轴。
%
% 输入参数：
%   scene       - 1×1 struct，由 rv_isar_build_scene 产出的统一场景结构体，至少包含：
%     .cfgISAR       1×1 struct，ISAR 波形、微动姿态与逐目标 frequencyBiasHz (N×1)
%     .isarSettings  1×1 struct，已解析并校验的后处理与 RCMC 参数（省略时自动调 isarpath_settings）
%     .targets       N×1 struct，基础散射体目标数组
%     .radarState    1×1 struct，雷达运动学与射频状态
%     .truth         1×1 struct，目标真值（至少含 R0 N×1、v_radial0 N×1，可选 motionDiagnostics）
%   priorMotion - (可选) 1×1 struct，外部传入的测距测速先验（省略或为空时使用 scene.truth）：
%     .rangeM      N×1 double，先验初始斜距 (m)
%     .velocityMps N×1 double，先验初始径向速度 (m/s)
%   options     - (可选) 1×1 struct，可选覆盖项：
%     .tauTxSec              1×1 double，发射附加时延 (s)，默认 0
%     .tauRefSec             1×1 double，本振附加时延 (s)，默认等于 tauTxSec
%     .rangeAxisMode         1×1 string，"CENTERED_OFFSET_M"（默认）或 "LOCAL_OFFSET_M"
%     .crossAxisMode         1×1 string，目前支持 "DOPPLER_M"（默认）
%     .dechirpChainMode      1×1 string，"ANALYTIC_IF"（默认）或 "REAL_RF"
%
% 输出参数：
%   isarScene   - 1×1 struct，ISAR 链路完整纯计算场景容器：
%     .cfg                       1×1 struct，ISAR 波形配置
%     .sharedParameters          1×1 struct，ISAR 后处理与 RCMC 参数（来自 scene.isarSettings）
%     .nTargets                  1×1 double，目标数 N
%     .fsHighHz, .fsScopeHz      1×1 double，高采样率与示波器/FPGA 侧降采样率 (Hz)
%     .NHigh, .NScope            1×1 double，单脉冲高采样/低采样点数
%     .fastTimeSec               1×NHigh double，高采样快时间轴 (s)
%     .pulseTimeSec              Na×1 double，绝对慢时间轴 (s)
%     .slowTimeRelativeSec       Na×1 double，相对首脉冲慢时间轴 (s)
%     .tauTxSec, .tauRefSec      1×1 double，发射与本振附加时延 (s)
%     .radarState                1×1 struct，ISAR 载频/波长下的雷达状态
%     .priorRangeM               N×1 double，生效的逐目标先验斜距 (m)
%     .priorVelocityMps          N×1 double，生效的逐目标先验径向速度 (m/s)
%     .priorSource               1×1 string，"supplied" 或 "truth"
%     .freqBiasHz                N×1 double，各目标拍频起始位置 (Hz，正=左侧，负=右侧)
%     .timeBiasSec               N×1 double，各目标本振相对目标回波的提前量 freqBiasHz/kr (s)
%     .referenceStartTimeSec     N×1 double，各目标局部本振起始时刻 (s)
%     .priors                    N×1 struct，逐目标先验结构体
%     .targetsISAR               N×1 struct，叠加微动姿态后的散射体目标数组
%     .microMotion               1×1 struct，微动姿态参数与诊断信息
%     .crossRangeGeometry        N×1 struct，各目标方位定标几何（effectiveRateRadps 等）
%     .scattererState            N×1 cell，各目标慢时间散射点状态（位置 Ns×3×Na，距离 Ns×Na）
%     .commonRangeAxisM          1×Nr double，公共距离偏移轴 (m)
%     .slowTimeDopplerAxisHz     1×Nc double，慢时间多普勒频率轴 (Hz)
%     .crossAxisM                N×1 cell，各目标方位坐标轴 (1×Nc double, m)
%     .motionDiagnostics         1×1 struct，初始 LOS 加速度残差与 RCMC 裕度诊断
%     .dechirpOptions            1×1 struct，默认去斜链参数（供 isarpath_dechirp 使用）

if nargin < 2
    priorMotion = [];
end
if nargin < 3 || isempty(options)
    options = struct();
end

cfgISAR = scene.cfgISAR;
targets = scene.targets;
nTarget = numel(targets);

fsHigh  = double(cfgISAR.fs);
fsScope = double(cfgISAR.fs_scope);
T_ISAR  = double(cfgISAR.T);
Tr_ISAR = double(cfgISAR.Tr);
Ta_ISAR = double(cfgISAR.Ta);
Na_ISAR = double(cfgISAR.Na);

tauTxSec         = double(get_opt(options, 'tauTxSec', 0));
tauRefSec        = double(get_opt(options, 'tauRefSec', tauTxSec));
rangeAxisMode    = upper(string(get_opt(options, 'rangeAxisMode', "CENTERED_OFFSET_M")));
crossAxisMode    = upper(string(get_opt(options, 'crossAxisMode', "DOPPLER_M")));
dechirpChainMode = string(get_opt(options, 'dechirpChainMode', "ANALYTIC_IF"));

%% 1. 获取由 rv_isar_build_scene 统一解析的后处理与 RCMC 参数
if isfield(scene, 'isarSettings') && isstruct(scene.isarSettings) && ~isempty(fieldnames(scene.isarSettings))
    sharedParameters = scene.isarSettings;
else
    sharedParameters = isarpath_settings(cfgISAR);
end

%% 2. 时间轴与 ISAR 载频雷达状态
NHigh  = round(fsHigh * T_ISAR);
NScope = round(fsScope * T_ISAR);
fastTimeSec         = (0:NHigh-1) ./ fsHigh;
pulseTimeSec        = double(cfgISAR.startTimeSec) + (0:Na_ISAR-1).' .* Ta_ISAR;
slowTimeRelativeSec = pulseTimeSec - pulseTimeSec(1);

radarStateISAR = scene.radarState;
radarStateISAR.carrierFrequencyHz = cfgISAR.fc;
radarStateISAR.wavelengthM        = cfgISAR.lambda;

%% 3. 先验运动参数与逐目标局部本振起始时刻
if ~isempty(priorMotion) && isstruct(priorMotion)
    priorRangeM      = double(priorMotion.rangeM(:));
    priorVelocityMps = double(priorMotion.velocityMps(:));
    priorSource      = "supplied";
else
    priorRangeM      = double(scene.truth.R0(:));
    priorVelocityMps = double(scene.truth.v_radial0(:));
    priorSource      = "truth";
end

if numel(priorRangeM) ~= nTarget || numel(priorVelocityMps) ~= nTarget
    error('RadarProject:ISARPriorSizeMismatch', ...
        '先验距离/速度数量 (%d/%d) 与目标数 (%d) 不一致。', ...
        numel(priorRangeM), numel(priorVelocityMps), nTarget);
end
if any(~isfinite(priorRangeM)) || any(~isfinite(priorVelocityMps))
    error('RadarProject:ISARPriorNotFinite', ...
        'priorRangeM 和 priorVelocityMps 中不能包含 NaN 或 Inf。');
end

if isfield(cfgISAR, 'frequencyBiasHz') && numel(cfgISAR.frequencyBiasHz) == nTarget
    freqBiasHz = double(cfgISAR.frequencyBiasHz(:));
elseif isfield(cfgISAR, 'frequencyBiasPatternHz') && ~isempty(cfgISAR.frequencyBiasPatternHz)
    pattern = double(cfgISAR.frequencyBiasPatternHz(:));
    freqBiasHz = pattern(mod(0:nTarget-1, numel(pattern)) + 1);
else
    pattern = [20e6; 45e6];
    freqBiasHz = pattern(mod(0:nTarget-1, numel(pattern)) + 1);
end
timeBiasSec           = freqBiasHz ./ cfgISAR.kr;
referenceStartTimeSec = 2 .* priorRangeM ./ cfgISAR.propagationSpeed - timeBiasSec;

if any(referenceStartTimeSec < 0) || any(referenceStartTimeSec + Tr_ISAR > T_ISAR)
    error('RadarProject:ISARReferenceOutsideRecord', ...
        '至少一个局部本振的完整 Tr 窗口超出单脉冲记录 [0, T]。');
end

phaseMode = upper(string(sharedParameters.phaseCompensationMode));
priors = repmat(struct( ...
    'epochSec', 0, 'rangeM', 0, 'velocityMps', 0, ...
    'accelerationMps2', 0, 'loStartSec', 0, 'beatOffsetHz', 0), nTarget, 1);
for targetIndex = 1:nTarget
    accPrior = 0;
    if isfield(cfgISAR, 'priorAccelerationMps2') && ~isempty(cfgISAR.priorAccelerationMps2)
        suppliedAcc = cfgISAR.priorAccelerationMps2;
        if isscalar(suppliedAcc)
            accPrior = double(suppliedAcc);
        else
            accPrior = double(suppliedAcc(targetIndex));
        end
    elseif phaseMode == "MODEL_CA"
        error('RadarProject:ISARAccelerationPriorMissing', ...
            'MODEL_CA requires cfgISAR.priorAccelerationMps2.');
    end
    priors(targetIndex) = struct( ...
        'epochSec',         pulseTimeSec(1), ...
        'rangeM',           priorRangeM(targetIndex), ...
        'velocityMps',      priorVelocityMps(targetIndex), ...
        'accelerationMps2', accPrior, ...
        'loStartSec',       referenceStartTimeSec(targetIndex), ...
        'beatOffsetHz',     freqBiasHz(targetIndex));
end

%% 4. 接入微动姿态、方位几何与慢时间散射点状态
[targetsISAR, microMotion] = isar_apply_micro_motion(targets, cfgISAR);
crossRangeGeometry = isar_cross_range_geometry( ...
    targetsISAR, radarStateISAR, cfgISAR, pulseTimeSec);

scattererState = cell(nTarget, 1);
for targetIndex = 1:nTarget
    scattererState{targetIndex} = scatterermodel.evaluateTarget( ...
        targetsISAR(targetIndex), radarStateISAR, pulseTimeSec);
end

%% 5. 初始 LOS 加速度残差与 RCMC 搜索半径检查
motionDiagnostics = struct();
if isfield(scene, 'truth') && isfield(scene.truth, 'motionDiagnostics') && ...
        ~isempty(scene.truth.motionDiagnostics)
    motionDiagnostics = scene.truth.motionDiagnostics;
    motionDiagnostics.rangeBinSpacingM = sharedParameters.rangeBinSpacingM;
    motionDiagnostics.accelerationResidualBins = ...
        motionDiagnostics.accelerationResidualM ./ sharedParameters.rangeBinSpacingM;
    motionDiagnostics.maximumAbsoluteResidualBins = max( ...
        abs(motionDiagnostics.accelerationResidualBins), [], 1).';
    if any(motionDiagnostics.maximumAbsoluteResidualBins > sharedParameters.rcmcSearchRadiusBins)
        error('RadarProject:ISARRcmcResidualOutsideSearch', ...
            'Initial-LOS acceleration residual exceeds the +/- %d bin RCMC search.', ...
            sharedParameters.rcmcSearchRadiusBins);
    end
end

%% 6. 公共距离轴与逐目标方位定标轴
switch rangeAxisMode
    case "CENTERED_OFFSET_M"
        commonRangeAxisM = ((0:sharedParameters.hrrfBins-1) - ...
            floor(sharedParameters.hrrfBins / 2)) .* sharedParameters.rangeBinSpacingM;
    case "LOCAL_OFFSET_M"
        commonRangeAxisM = (0:sharedParameters.hrrfBins-1) .* sharedParameters.rangeBinSpacingM;
    otherwise
        error('RadarProject:ISARUnsupportedRangeAxisMode', ...
            '不支持的 rangeAxisMode="%s"。', rangeAxisMode);
end

nCrossBins = Na_ISAR * sharedParameters.slowTimeZeroPaddingFactor;
switch crossAxisMode
    case "DOPPLER_M"
        slowTimeDopplerAxisHz = ((0:nCrossBins-1) - nCrossBins / 2) ./ (nCrossBins * Ta_ISAR);
    otherwise
        error('RadarProject:ISARUnsupportedCrossAxisMode', ...
            '当前 crossAxisM 只支持 crossAxisMode="DOPPLER_M"。');
end

crossAxisM = cell(nTarget, 1);
for targetIndex = 1:nTarget
    crossAxisM{targetIndex} = (cfgISAR.propagationSpeed / cfgISAR.fc) .* ...
        slowTimeDopplerAxisHz ./ (2 * crossRangeGeometry(targetIndex).effectiveRateRadps);
end

scopeLpfEdgeHz = fsScope / 2;
dechirpOptions = struct( ...
    'Mode',                       dechirpChainMode, ...
    'RemoveMirrorFrequency',      true, ...
    'RemoveSamplingAlias',        true, ...
    'RemoveTimeOverlapCrosstalk', true, ...
    'RemoveDc',                   true, ...
    'LowpassCenterHz',            scopeLpfEdgeHz / 2, ...
    'LowpassBandwidthHz',         scopeLpfEdgeHz, ...
    'FilterSuppressionDb',        120, ...
    'ReturnPerTarget',            true);

isarScene = struct( ...
    'cfg',                   cfgISAR, ...
    'sharedParameters',      sharedParameters, ...
    'nTargets',              nTarget, ...
    'fsHighHz',              fsHigh, ...
    'fsScopeHz',             fsScope, ...
    'NHigh',                 NHigh, ...
    'NScope',                NScope, ...
    'fastTimeSec',           fastTimeSec, ...
    'pulseTimeSec',          pulseTimeSec, ...
    'slowTimeRelativeSec',   slowTimeRelativeSec, ...
    'tauTxSec',              tauTxSec, ...
    'tauRefSec',             tauRefSec, ...
    'radarState',            radarStateISAR, ...
    'priorRangeM',           priorRangeM, ...
    'priorVelocityMps',      priorVelocityMps, ...
    'priorSource',           priorSource, ...
    'freqBiasHz',            freqBiasHz, ...
    'timeBiasSec',           timeBiasSec, ...
    'referenceStartTimeSec', referenceStartTimeSec, ...
    'priors',                priors, ...
    'targetsISAR',           targetsISAR, ...
    'microMotion',           microMotion, ...
    'crossRangeGeometry',    crossRangeGeometry, ...
    'scattererState',        {scattererState}, ...
    'commonRangeAxisM',      commonRangeAxisM, ...
    'slowTimeDopplerAxisHz', slowTimeDopplerAxisHz, ...
    'crossAxisM',            {crossAxisM}, ...
    'motionDiagnostics',     motionDiagnostics, ...
    'dechirpOptions',        dechirpOptions);
end

function val = get_opt(s, name, defaultVal)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    val = s.(name);
else
    val = defaultVal;
end
end
