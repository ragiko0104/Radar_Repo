function isarRes = isarpath_run(scene, priorMotion, options)
%ISARPATH_RUN ISAR 成像链路（ISAR_path）纯计算八步调度总入口。
%
% 功能说明：
%   接收由 rv_isar_build_scene 构建的统一场景结构体 scene 及可选先验运动参数，
%   按顺序执行 8 个纯计算步骤（不含绘图与质量评估）：
%     1. isarpath_build_scene           : 场景产生（时间轴、先验、微动、几何与坐标轴）
%     2. isarpath_generate_echo         : 高采样宽带复 LFM 回波产生
%     3. isarpath_dechirp               : 局部本振去斜、滤波与降采样
%     4. isarpath_window_and_dispersion : 逐目标时域截窗与多普勒色散补偿
%     5. isarpath_form_hrrp             : 距离维 FFT 与先验随动截取形成 HRRP-慢时间矩阵
%                                         （同步保留色散补偿前 hrrfBeforeDispersion）
%     6. isarpath_align_profiles        : 包络距离徙动校正（RCMC）
%     7. isarpath_phase_compensate      : 特显点/先验模型脉冲间相位聚焦补偿
%     8. isarpath_azimuth_fft           : 方位维加窗零填充 FFT 成像
%   注：指标评估请调用 functions/assess/ 或 code/assess_rv_isar_results.m，
%       图形绘制请调用 functions/plot/ 或 code/plot_rv_isar_results.m。
%
% 输入参数：
%   scene       - 1×1 struct，由 rv_isar_build_scene 产出的统一场景结构体
%   priorMotion - (可选) 1×1 struct，外部传入的先验运动参数（如 RV 链路校准结果）；
%                 省略或为空时回退至 scene.truth：
%     .rangeM      N×1 double，先验初始斜距 (m)
%     .velocityMps N×1 double，先验初始径向速度 (m/s)
%   options     - (可选) 1×1 struct，透传给各子步骤的可选控制项：
%     .sceneOptions   1×1 struct，透传给 isarpath_build_scene
%     .echoOptions    1×1 struct，透传给 isarpath_generate_echo
%     .dechirpOptions 1×1 struct，透传给 isarpath_dechirp
%
% 输出参数：
%   isarRes     - 1×1 struct，ISAR 链路纯计算输出容器：
%     .scene              1×1 struct，步骤 1 产出的 isarScene 结构体
%     .cfg                1×1 struct，ISAR 波形配置
%     .sharedParameters   1×1 struct，ISAR 后处理与 RCMC 参数
%     .targetsISAR        N×1 struct，叠加微动姿态后的散射体目标数组
%     .motion             1×1 struct，微动姿态诊断信息
%     .crossRangeGeometry N×1 struct，各目标方位定标几何
%     .scattererState     N×1 cell，各目标慢时间散射点状态（Ns×3×Na 等）
%     .echo               1×1 struct，步骤 2 产出的高采样率回波
%     .echoInfo           1×1 struct，回波元信息
%     .beat               1×1 struct，步骤 3 产出的窄带去斜拍频
%     .beatInfo           1×1 struct，去斜元信息
%     .post               1×1 struct，步骤 4–8 产出的逐目标中间量与最终成像容器：
%       .hrrfBeforeDispersion  N×1 cell，色散补偿前 HRRP-慢时间矩阵（Na×Nr complex）
%       .hrrfModelCrop         N×1 cell，色散补偿后 HRRP-慢时间矩阵（Na×Nr complex）
%       .hrrfAligned           N×1 cell，RCMC 对齐后 HRRP-慢时间矩阵（Na×Nr complex）
%       .phaseCompensatedHrrf  N×1 cell，相位补偿后 HRRP-慢时间矩阵（Na×Nr complex）
%       .image                 N×1 cell，最终归一化 ISAR 幅度像（Nc×Nr double）
%       .rangeAxisM            N×1 cell，距离坐标轴（1×Nr double, m）
%       .crossAxisM            N×1 cell，方位坐标轴（1×Nc double, m）
%       .rangeMigrationBins    N×Na double，各脉冲 RCMC 移位量 (bin)
%       .valid                 N×1 logical，各目标成像有效性标志

if nargin < 2
    priorMotion = [];
end
if nargin < 3 || isempty(options)
    options = struct();
end

sceneOpts   = get_opt(options, 'sceneOptions', struct());
echoOpts    = get_opt(options, 'echoOptions', struct());
dechirpOpts = get_opt(options, 'dechirpOptions', struct());

%% 步骤 1：ISAR 场景产生
isarScene = isarpath_build_scene(scene, priorMotion, sceneOpts);
p         = isarScene.sharedParameters;
nTarget   = isarScene.nTargets;
Na_ISAR   = isarScene.cfg.Na;

%% 步骤 2：高采样宽带复 LFM 回波产生
[isarEcho, isarEchoInfo] = isarpath_generate_echo(isarScene, echoOpts);

%% 步骤 3：局部本振去斜与降采样（dechirp）
[isarBeat, isarBeatInfo] = isarpath_dechirp(isarEcho, isarScene, dechirpOpts);

%% 预分配后处理输出容器
isarPost = struct();
isarPost.hrrfBeforeDispersion      = cell(nTarget, 1);
isarPost.hrrfModelCrop             = cell(nTarget, 1);
isarPost.hrrfAligned               = cell(nTarget, 1);
isarPost.rcmcCorrelationScore      = cell(nTarget, 1);
isarPost.phaseCompensatedHrrf      = cell(nTarget, 1);
isarPost.phaseCorrection           = cell(nTarget, 1);
isarPost.phaseReferenceDiagnostics = cell(nTarget, 1);
isarPost.image                     = cell(nTarget, 1);
isarPost.rangeAxisM                = cell(nTarget, 1);
isarPost.crossAxisM                = cell(nTarget, 1);
isarPost.segmentStartSamples       = nan(nTarget, Na_ISAR);
isarPost.cropCenterBins            = nan(nTarget, Na_ISAR);
isarPost.rangeMigrationBins        = nan(nTarget, Na_ISAR);
isarPost.phaseReferenceMode        = strings(nTarget, 1);
isarPost.valid                     = false(nTarget, 1);
isarPost.parameters = struct( ...
    'inputSource',                   "perTarget", ...
    'inputSampleRateHz',             isarScene.fsScopeHz, ...
    'timeWindowSamples',             p.timeWindowSamples, ...
    'rangeFftPoints',                p.rangeFftPoints, ...
    'hrrfBins',                      p.hrrfBins, ...
    'frequencyBinSpacingHz',         p.frequencyBinSpacingHz, ...
    'rangeBinSpacingM',              p.rangeBinSpacingM, ...
    'cropCenterOffsetBinsEffective', p.cropCenterOffsetBinsEffective, ...
    'timeWindowMode',                upper(string(p.timeWindowMode)), ...
    'rangeMigrationMode',            upper(string(p.rangeMigrationMode)), ...
    'rangeAlignmentMode',            p.rangeAlignmentMode, ...
    'rcmcCorrelationPoints',         p.rcmcCorrelationPoints, ...
    'rcmcSearchRadiusBins',          p.rcmcSearchRadiusBins, ...
    'rcmcReferencePulseOneBased',    1, ...
    'phaseCompensationMode',         upper(string(p.phaseCompensationMode)), ...
    'dispersionCompensationEnabled', p.dispersionCompensationEnabled, ...
    'windowType',                    p.windowType, ...
    'slowTimeZeroPaddingFactor',     p.slowTimeZeroPaddingFactor);
isarPost.microMotion        = isarScene.microMotion;
isarPost.crossRangeGeometry = isarScene.crossRangeGeometry;
if ~isempty(fieldnames(isarScene.motionDiagnostics))
    isarPost.motionDiagnostics = isarScene.motionDiagnostics;
end

nBeatSample = size(isarBeat.total, 2);

%% 步骤 4–8：逐目标执行时域截窗/色散补偿 -> HRRP -> RCMC -> 特显点补偿 -> 方位维 FFT
for targetIndex = 1:nTarget
    beatSignal = isarBeat.perTarget(targetIndex).signal;
    if size(beatSignal, 2) ~= nBeatSample
        error('RadarProject:ISARBeatSampleCountMismatch', ...
            '目标 %d 的逐目标拍频快时间点数与场景总拍频不一致。', targetIndex);
    end
    prior = isarScene.priors(targetIndex);

    % 步骤 4：时域截窗与多普勒色散补偿
    [segmentComp, segmentRaw, segmentStarts] = isarpath_window_and_dispersion( ...
        beatSignal, isarScene.cfg, prior, isarScene.pulseTimeSec, p);

    % 步骤 5：形成 HRRP-慢时间矩阵（始终保留色散补偿前 HRRP 供独立评估与绘图）
    [modelHrrf, cropCentres, hrrfBefore] = isarpath_form_hrrp( ...
        segmentComp, isarScene.cfg, prior, isarScene.pulseTimeSec, p, segmentRaw);
    clear segmentComp segmentRaw

    % 步骤 6：RCMC 包络距离徙动校正
    [alignedHrrf, rangeMigrationBins, rcmcScore] = isarpath_align_profiles(modelHrrf, p);

    % 步骤 7：特显点 / 先验模型脉冲间相位补偿
    [phaseCompHrrf, phaseCorr, phaseInfo] = isarpath_phase_compensate( ...
        alignedHrrf, p, isarScene.cfg, prior, isarScene.pulseTimeSec, rangeMigrationBins);

    % 步骤 8：方位维 FFT 成像
    [isarImage, isValid] = isarpath_azimuth_fft(phaseCompHrrf, p);

    isarPost.segmentStartSamples(targetIndex, :)    = segmentStarts(:).';
    isarPost.cropCenterBins(targetIndex, :)         = cropCentres(:).';
    isarPost.hrrfBeforeDispersion{targetIndex}      = hrrfBefore;
    isarPost.hrrfModelCrop{targetIndex}             = modelHrrf;
    isarPost.rangeMigrationBins(targetIndex, :)     = rangeMigrationBins(:).';
    isarPost.rcmcCorrelationScore{targetIndex}      = rcmcScore;
    isarPost.hrrfAligned{targetIndex}               = alignedHrrf;
    isarPost.phaseCompensatedHrrf{targetIndex}      = phaseCompHrrf;
    isarPost.phaseCorrection{targetIndex}           = phaseCorr;
    isarPost.phaseReferenceMode(targetIndex)        = phaseInfo.mode;
    isarPost.phaseReferenceDiagnostics{targetIndex} = phaseInfo;
    isarPost.image{targetIndex}                     = isarImage;
    isarPost.rangeAxisM{targetIndex}                = isarScene.commonRangeAxisM;
    isarPost.crossAxisM{targetIndex}                = isarScene.crossAxisM{targetIndex};
    isarPost.valid(targetIndex)                     = isValid;
end

isarRes = struct( ...
    'scene',              isarScene, ...
    'cfg',                isarScene.cfg, ...
    'sharedParameters',   p, ...
    'targetsISAR',        isarScene.targetsISAR, ...
    'motion',             isarScene.microMotion, ...
    'crossRangeGeometry', isarScene.crossRangeGeometry, ...
    'scattererState',     {isarScene.scattererState}, ...
    'echo',               isarEcho, ...
    'echoInfo',           isarEchoInfo, ...
    'beat',               isarBeat, ...
    'beatInfo',           isarBeatInfo, ...
    'post',               isarPost);
end

function val = get_opt(s, name, defaultVal)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    val = s.(name);
else
    val = defaultVal;
end
end
