function rvMeasurement = rvpath_solve_rv(rvPulseCompression, rvCfar, wf, truth, constants, options)
%RVPATH_SOLVE_RV 基于六段 LFM 寻峰结果执行一致性匹配、距速解算与偏置校准。
%
% 功能说明：
%   1. 计算固定初始 LOS 投影下的慢时间理论真值距离与径向速度；
%   2. 每脉冲选取平均峰功率最高的 3 段 LFM，估计各峰有效 SNR（支持 localCfar 与
%      robustOffTarget 两种噪声估计器）及矩形谱 LFM 时延 CRLB 标准差；
%   3. 调用 radar.Cal_rv_mat 联合求解 (1,2) 与 (1,3) 两对频带的距速候选矩阵，
%      通过解析雅可比传播时延误差形成一致性门限（snrPrecision 或 legacyResolution），
%      并结合速度上限门限与一对一峰独占约束筛选有效目标；
%   4. 按距离排序与理论真值配对，使用第 1 脉冲有效估计的均值偏置对全 CPI 做常量校准，
%      输出原始与校准后的误差统计及全部中间诊断字段。
%
% 输入参数：
%   rvPulseCompression - 1×1 struct，由 rvpath_pulse_compress 产出的脉压结果
%                        （使用 compressed{6×1}，每元素为 Na×N complex）
%   rvCfar             - 1×1 struct，由 rvpath_cfar_peaks 产出的 CFAR 与寻峰结果
%                        （使用 segmentKeys/segmentBandIndex/parameters/threshold/
%                        refinedPeaks/meanPeakPower）
%   wf                 - 1×1 struct，由 rvpath_derive_waveform 产出的派生 RV 波形结构体
%                        （使用 fs/Tr/fc/B1/B2/B3/fc_v/k_v/N/Na/slowTimeSec/rvConsistency）
%   truth              - 1×1 struct，目标初始径向真值（与 scene.truth 兼容）：
%     .R0          Ntarget×1 double，t=0 初始斜距 (m)
%     .v_radial0   Ntarget×1 double，t=0 初始径向速度 (m/s)
%     .a_radial    Ntarget×1 double，径向加速度 (m/s^2)
%   constants          - 1×1 struct，物理常量（必须包含 .c 1×1 double 与 .Mach 1×1 double）
%   options            - (可选) 1×1 struct，可选控制项：
%     .printSummary  1×1 logical，是否向命令行打印解算汇总（默认 true）
%
% 输出参数：
%   rvMeasurement - 1×1 struct，六段 LFM 距速解算与评估容器，主要字段：
%     .parameters                          1×1 struct，门限、噪声估计器与校准参数
%     .selectedSegmentIndices              Na×3 double，每脉冲选中的 3 个 LFM 段编号
%     .rangeCandidatesM                    Na×1 cell，每脉冲通过门限的候选距离行向量 (m)
%     .velocityCandidatesMps               Na×1 cell，每脉冲通过门限的候选速度行向量 (m/s)
%     .rangeM, .velocityMps                Na×Ntarget double，配对后的原始估计距离/速度
%     .truthRangeM, .truthVelocityMps      Na×Ntarget double，理论真值距离/速度
%     .matchedCount                        Na×1 double，每脉冲匹配目标数
%     .missedTargetCount                   Na×1 double，每脉冲漏检数
%     .falseEstimateCount                  Na×1 double，每脉冲额外估计数
%     .selectedPeakSnrDb                   Na×3 cell，选中三段各峰 SNR (dB)
%     .selectedPeakNoisePower              Na×3 cell，选中三段生效噪声功率
%     .selectedPeakLocalCfarNoisePower     Na×3 cell，局部 CFAR 噪声功率诊断
%     .selectedPeakRobustOffTargetNoisePower Na×3 cell，离目标稳健噪声功率诊断
%     .selectedPeakDelayStdSec             Na×3 cell，各峰时延标准差 (s)
%     .matchThirdPeakIndex                 Na×1 cell，匹配到的第 3 段峰索引矩阵
%     .matchRangeResidualM                 Na×1 cell，两对解距离残差矩阵 (m)
%     .matchVelocityResidualMps            Na×1 cell，两对解速度残差矩阵 (m/s)
%     .matchRangeThresholdM                Na×1 cell，距离一致性门限矩阵 (m)
%     .matchVelocityThresholdMps           Na×1 cell，速度一致性门限矩阵 (m/s)
%     .matchNormalizedScore                Na×1 cell，归一化失配得分矩阵
%     .preSpeedGateCandidateCount          Na×1 double，速度门限前候选数
%     .postSpeedGateCandidateCount         Na×1 double，速度门限后候选数
%     .postOneToOneCandidateCount          Na×1 double，一对一筛选后候选数
%     .consistencyRangeBiasM               Na×1 double，跨频带距离公共偏置 (m)
%     .consistencyVelocityBiasMps          Na×1 double，跨频带速度公共偏置 (m/s)
%     .calibrationSampleCount              1×1 double，第 1 脉冲参与偏置校准的目标数
%     .calibrationAvailable                1×1 logical，是否存在有效校准样本
%     .rangeMeanBiasM, .velocityMeanBiasMps 1×1 double，第 1 脉冲均值偏置
%     .calibratedRangeM                    Na×Ntarget double，校准后距离估计 (m)
%     .calibratedVelocityMps               Na×Ntarget double，校准后速度估计 (m/s)
%     .rawRangeErrorM, .rawVelocityErrorMps Na×Ntarget double，原始绝对误差
%     .calibratedRangeErrorM               Na×Ntarget double，校准后距离绝对误差 (m)
%     .calibratedVelocityErrorMps          Na×Ntarget double，校准后速度绝对误差 (m/s)
%     .rawRangeRelativeError               Na×Ntarget double，原始距离相对误差
%     .rawVelocityRelativeError            Na×Ntarget double，原始速度相对误差
%     .calibratedRangeRelativeError        Na×Ntarget double，校准后距离相对误差
%     .calibratedVelocityRelativeError     Na×Ntarget double，校准后速度相对误差
%     .maxRawRangeRelativeErrorPercent     1×1 double，原始测距最大相对误差 (%)
%     .maxRawVelocityRelativeErrorPercent  1×1 double，原始测速最大相对误差 (%)
%     .maxCalibratedRangeRelativeErrorPercent    1×1 double，校准后测距最大相对误差 (%)
%     .maxCalibratedVelocityRelativeErrorPercent 1×1 double，校准后测速最大相对误差 (%)

if nargin < 6 || isempty(options)
    options = struct();
end
printSummary = true;
if isfield(options, 'printSummary') && ~isempty(options.printSummary)
    printSummary = logical(options.printSummary);
end

c    = double(constants.c);
Mach = double(constants.Mach);
fs   = double(wf.fs);
Tr   = double(wf.Tr);
fc   = double(wf.fc);
B1   = double(wf.B1);
B2   = double(wf.B2);
B3   = double(wf.B3);
fc_v = double(wf.fc_v);
k_v  = double(wf.k_v);
N    = double(wf.N);
Na   = double(wf.Na);
slowTimeSecRV = double(wf.slowTimeSec(:));

R0        = double(truth.R0(:));
v_radial0 = double(truth.v_radial0(:));
a_radial  = double(truth.a_radial(:));
rvTargetCount = numel(R0);
if numel(v_radial0) ~= rvTargetCount || numel(a_radial) ~= rvTargetCount
    error('multiuav:rvpath:invalidTruthSize', ...
        'R0、v_radial0 和 a_radial 必须包含相同数量的目标。');
end

rvTruthRangeM = R0.' + slowTimeSecRV .* v_radial0.' + ...
    0.5 .* slowTimeSecRV.^2 .* a_radial.';
rvTruthVelocityMps = v_radial0.' + slowTimeSecRV .* a_radial.';

rvLegacyRangeConsistencyM      = c ./ min(abs([B1, B2, B3])) ./ 5;
rvLegacyVelocityConsistencyMps = (c ./ 2 ./ (Tr .* fc)) .* 2 ./ 5;
rvMaxAbsVelocityMps            = 20 .* Mach;

%% 1. 解析并校验一致性门限配置
if isfield(wf, 'rvConsistency') && ~isempty(wf.rvConsistency)
    rvConsistencyConfig = wf.rvConsistency;
else
    rvConsistencyConfig = struct();
end
rvConsistencyMode = char(string(get_field_or_default( ...
    rvConsistencyConfig, 'mode', 'snrPrecision')));
if ~any(strcmpi(rvConsistencyMode, {'snrPrecision', 'legacyResolution'}))
    error('multiuav:rvpath:invalidConsistencyMode', ...
        'cfg_RVpath.rvConsistency.mode must be snrPrecision or legacyResolution.');
end
rvConsistencySigmaMultiplier = double(get_field_or_default( ...
    rvConsistencyConfig, 'sigmaMultiplier', 2));
rvConsistencySnrFloorDb = double(get_field_or_default( ...
    rvConsistencyConfig, 'snrFloorDb', 0));
rvConsistencySnrCeilingDb = double(get_field_or_default( ...
    rvConsistencyConfig, 'snrCeilingDb', 60));
rvConsistencyDelayStdFloorBins = double(get_field_or_default( ...
    rvConsistencyConfig, 'delayStdFloorBins', 0.05));
rvConsistencyBiasCorrection = char(string(get_field_or_default( ...
    rvConsistencyConfig, 'biasCorrection', 'rankMedian')));
rvConsistencyNoiseEstimator = char(string(get_field_or_default( ...
    rvConsistencyConfig, 'noiseEstimator', 'localCfar')));
rvConsistencyMaxPeakRankDifference = double(get_field_or_default( ...
    rvConsistencyConfig, 'maxPeakRankDifference', 0));
rvConsistencyEnforceOneToOne = logical(get_field_or_default( ...
    rvConsistencyConfig, 'enforceOneToOnePeakUse', true));

if ~(isscalar(rvConsistencySigmaMultiplier) && ...
        isfinite(rvConsistencySigmaMultiplier) && ...
        rvConsistencySigmaMultiplier > 0)
    error('multiuav:rvpath:invalidConsistencySigmaMultiplier', ...
        'rvConsistency.sigmaMultiplier must be a positive finite scalar.');
end
if ~(isscalar(rvConsistencySnrFloorDb) && isfinite(rvConsistencySnrFloorDb) && ...
        isscalar(rvConsistencySnrCeilingDb) && isfinite(rvConsistencySnrCeilingDb) && ...
        rvConsistencySnrCeilingDb >= rvConsistencySnrFloorDb)
    error('multiuav:rvpath:invalidConsistencySnrRange', ...
        'rvConsistency SNR limits must be finite and ceiling >= floor.');
end
if ~(isscalar(rvConsistencyDelayStdFloorBins) && ...
        isfinite(rvConsistencyDelayStdFloorBins) && ...
        rvConsistencyDelayStdFloorBins >= 0)
    error('multiuav:rvpath:invalidDelayPrecisionFloor', ...
        'rvConsistency.delayStdFloorBins must be a finite nonnegative scalar.');
end
if ~any(strcmpi(rvConsistencyBiasCorrection, {'none', 'rankMedian'}))
    error('multiuav:rvpath:invalidConsistencyBiasCorrection', ...
        'rvConsistency.biasCorrection must be none or rankMedian.');
end
if ~any(strcmpi(rvConsistencyNoiseEstimator, {'localCfar', 'robustOffTarget'}))
    error('multiuav:rvpath:invalidConsistencyNoiseEstimator', ...
        'rvConsistency.noiseEstimator must be localCfar or robustOffTarget.');
end
if ~(isscalar(rvConsistencyMaxPeakRankDifference) && ...
        rvConsistencyMaxPeakRankDifference >= 0 && ...
        (isinf(rvConsistencyMaxPeakRankDifference) || ...
        rvConsistencyMaxPeakRankDifference == floor(rvConsistencyMaxPeakRankDifference)))
    error('multiuav:rvpath:invalidPeakRankDifference', ...
        'rvConsistency.maxPeakRankDifference must be a nonnegative integer or Inf.');
end
if ~isscalar(rvConsistencyEnforceOneToOne)
    error('multiuav:rvpath:invalidOneToOneOption', ...
        'rvConsistency.enforceOneToOnePeakUse must be a logical scalar.');
end

rvConsistencySnrFloorLinear   = 10.^(rvConsistencySnrFloorDb ./ 10);
rvConsistencySnrCeilingLinear = 10.^(rvConsistencySnrCeilingDb ./ 10);
rvConsistencyDelayStdFloorSec = rvConsistencyDelayStdFloorBins ./ fs;

%% 2. 预分配输出结构体
rvMeasurement = struct();
rvMeasurement.parameters = struct( ...
    'segmentKeys',                  rvCfar.segmentKeys, ...
    'segmentCarrierFrequencyHz',    fc_v, ...
    'segmentChirpSlopeHzPerSec',    k_v, ...
    'selectedSegmentCount',         3, ...
    'peakIndexToDelayConvention',   '(peakIndex-1)/fs', ...
    'consistencyGateMode',          rvConsistencyMode, ...
    'consistencySigmaMultiplier',   rvConsistencySigmaMultiplier, ...
    'consistencySnrFloorDb',        rvConsistencySnrFloorDb, ...
    'consistencySnrCeilingDb',      rvConsistencySnrCeilingDb, ...
    'consistencyDelayStdFloorBins', rvConsistencyDelayStdFloorBins, ...
    'consistencyNoiseEstimator',    rvConsistencyNoiseEstimator, ...
    'consistencyBiasCorrection',    rvConsistencyBiasCorrection, ...
    'maxPeakRankDifference',        rvConsistencyMaxPeakRankDifference, ...
    'enforceOneToOnePeakUse',       rvConsistencyEnforceOneToOne, ...
    'legacyRangeConsistencyM',      rvLegacyRangeConsistencyM, ...
    'legacyVelocityConsistencyMps', rvLegacyVelocityConsistencyMps, ...
    'rangeConsistencyM',            rvLegacyRangeConsistencyM, ...
    'velocityConsistencyMps',       rvLegacyVelocityConsistencyMps, ...
    'snrEstimator',                 ['peakExcessOver' rvConsistencyNoiseEstimator], ...
    'delayPrecisionModel',          'LFM rectangular-spectrum delay CRLB', ...
    'maxAbsVelocityMps',            rvMaxAbsVelocityMps, ...
    'calibrationPulseIndex',        1, ...
    'truthConvention',              'fixedInitialLosProjection');

rvMeasurement.selectedSegmentIndices                 = nan(Na, 3);
rvMeasurement.rangeCandidatesM                       = cell(Na, 1);
rvMeasurement.velocityCandidatesMps                  = cell(Na, 1);
rvMeasurement.rangeM                                 = nan(Na, rvTargetCount);
rvMeasurement.velocityMps                            = nan(Na, rvTargetCount);
rvMeasurement.truthRangeM                            = rvTruthRangeM;
rvMeasurement.truthVelocityMps                       = rvTruthVelocityMps;
rvMeasurement.matchedCount                           = zeros(Na, 1);
rvMeasurement.missedTargetCount                      = rvTargetCount .* ones(Na, 1);
rvMeasurement.falseEstimateCount                     = zeros(Na, 1);
rvMeasurement.selectedPeakSnrDb                      = cell(Na, 3);
rvMeasurement.selectedPeakNoisePower                 = cell(Na, 3);
rvMeasurement.selectedPeakLocalCfarNoisePower        = cell(Na, 3);
rvMeasurement.selectedPeakRobustOffTargetNoisePower  = cell(Na, 3);
rvMeasurement.selectedPeakDelayStdSec                = cell(Na, 3);
rvMeasurement.matchThirdPeakIndex                    = cell(Na, 1);
rvMeasurement.matchRangeResidualM                    = cell(Na, 1);
rvMeasurement.matchVelocityResidualMps               = cell(Na, 1);
rvMeasurement.matchRangeThresholdM                   = cell(Na, 1);
rvMeasurement.matchVelocityThresholdMps              = cell(Na, 1);
rvMeasurement.matchNormalizedScore                   = cell(Na, 1);
rvMeasurement.preSpeedGateCandidateCount             = zeros(Na, 1);
rvMeasurement.postSpeedGateCandidateCount            = zeros(Na, 1);
rvMeasurement.postOneToOneCandidateCount             = zeros(Na, 1);
rvMeasurement.consistencyRangeBiasM                  = zeros(Na, 1);
rvMeasurement.consistencyVelocityBiasMps             = zeros(Na, 1);

%% 3. 逐脉冲三段联立解算与一致性筛选
rvSelectedBandwidthHz = [B1, B2, B3];
for rvPulseIndex = 1:Na
    [~, rvSelectedSegments] = maxk(rvCfar.meanPeakPower(rvPulseIndex, :), 3);
    rvMeasurement.selectedSegmentIndices(rvPulseIndex, :) = rvSelectedSegments;

    rvPeak1 = rvCfar.refinedPeaks{rvSelectedSegments(1)}{rvPulseIndex};
    rvPeak2 = rvCfar.refinedPeaks{rvSelectedSegments(2)}{rvPulseIndex};
    rvPeak3 = rvCfar.refinedPeaks{rvSelectedSegments(3)}{rvPulseIndex};
    if isempty(rvPeak1) || isempty(rvPeak2) || isempty(rvPeak3)
        rvMeasurement.rangeCandidatesM{rvPulseIndex}      = zeros(1, 0);
        rvMeasurement.velocityCandidatesMps{rvPulseIndex} = zeros(1, 0);
        continue;
    end

    rvSelectedPeakLists   = {rvPeak1, rvPeak2, rvPeak3};
    rvSelectedDelayStdSec = cell(1, 3);
    for rvSelectedSlot = 1:3
        rvSelectedSegment = rvSelectedSegments(rvSelectedSlot);
        rvSelectedBand    = rvCfar.segmentBandIndex(rvSelectedSegment);
        rvSelectedPeakIndex = max(1, min(N, round(double( ...
            rvSelectedPeakLists{rvSelectedSlot}(:).'))));
        rvSelectedPowerRow = abs(rvPulseCompression.compressed{ ...
            rvSelectedSegment}(rvPulseIndex, :)).^2;
        rvSelectedPeakPower = rvSelectedPowerRow(rvSelectedPeakIndex);
        rvSelectedLocalCfarNoisePower = rvCfar.threshold{rvSelectedSegment}( ...
            rvPulseIndex, rvSelectedPeakIndex) ./ rvCfar.parameters.alpha;

        rvSelectedNoiseMask = false(1, N);
        rvSelectedFirstNoiseCell = rvCfar.parameters.trainingCellsPerSide + ...
            rvCfar.parameters.guardCellsByBand(rvSelectedBand) + 1;
        rvSelectedLastNoiseCell = N - ...
            rvCfar.parameters.trainingCellsPerSide - ...
            rvCfar.parameters.guardCellsByBand(rvSelectedBand);
        rvSelectedNoiseMask(rvSelectedFirstNoiseCell:rvSelectedLastNoiseCell) = true;
        rvSelectedExclusionHalfWidth = ...
            rvCfar.parameters.guardCellsByBand(rvSelectedBand) + ...
            rvCfar.parameters.peakSeparationBinsByBand(rvSelectedBand);
        for rvNoisePeakIndex = 1:numel(rvSelectedPeakIndex)
            rvNoiseExclusionStart = max(rvSelectedFirstNoiseCell, ...
                rvSelectedPeakIndex(rvNoisePeakIndex) - rvSelectedExclusionHalfWidth);
            rvNoiseExclusionStop = min(rvSelectedLastNoiseCell, ...
                rvSelectedPeakIndex(rvNoisePeakIndex) + rvSelectedExclusionHalfWidth);
            rvSelectedNoiseMask(rvNoiseExclusionStart:rvNoiseExclusionStop) = false;
        end
        rvSelectedNoiseSamples = rvSelectedPowerRow(rvSelectedNoiseMask);
        rvSelectedNoiseSamples = rvSelectedNoiseSamples( ...
            isfinite(rvSelectedNoiseSamples) & rvSelectedNoiseSamples >= 0);
        if isempty(rvSelectedNoiseSamples)
            rvSelectedRobustNoisePower = median(rvSelectedLocalCfarNoisePower, 'omitnan');
        else
            rvSelectedRobustNoisePower = median(rvSelectedNoiseSamples) ./ log(2);
        end
        if ~(isfinite(rvSelectedRobustNoisePower) && rvSelectedRobustNoisePower >= 0)
            rvSelectedRobustNoisePower = 0;
        end

        if strcmpi(rvConsistencyNoiseEstimator, 'localCfar')
            rvSelectedNoisePower = rvSelectedLocalCfarNoisePower;
        else
            rvSelectedNoisePower = rvSelectedRobustNoisePower .* ones(size(rvSelectedPeakPower));
        end
        rvSelectedNoiseFloor = rvSelectedPeakPower ./ max(rvConsistencySnrCeilingLinear, 1);
        rvSelectedNoisePower = max(rvSelectedNoisePower, rvSelectedNoiseFloor);
        rvSelectedPeakSnr = max(rvSelectedPeakPower ./ max(rvSelectedNoisePower, realmin) - 1, ...
            rvConsistencySnrFloorLinear);
        rvSelectedPeakSnr = min(rvSelectedPeakSnr, rvConsistencySnrCeilingLinear);

        rvSelectedBetaRmsHz = abs(rvSelectedBandwidthHz(rvSelectedBand)) ./ sqrt(12);
        rvSelectedDelayCrlbStdSec = 1 ./ (2 .* sqrt(2) .* pi .* ...
            rvSelectedBetaRmsHz .* sqrt(rvSelectedPeakSnr));
        rvSelectedDelayStdSec{rvSelectedSlot} = hypot( ...
            rvSelectedDelayCrlbStdSec, rvConsistencyDelayStdFloorSec);

        rvMeasurement.selectedPeakSnrDb{rvPulseIndex, rvSelectedSlot} = ...
            10 .* log10(rvSelectedPeakSnr);
        rvMeasurement.selectedPeakNoisePower{rvPulseIndex, rvSelectedSlot} = ...
            rvSelectedNoisePower;
        rvMeasurement.selectedPeakLocalCfarNoisePower{rvPulseIndex, rvSelectedSlot} = ...
            rvSelectedLocalCfarNoisePower;
        rvMeasurement.selectedPeakRobustOffTargetNoisePower{rvPulseIndex, rvSelectedSlot} = ...
            rvSelectedRobustNoisePower .* ones(size(rvSelectedPeakPower));
        rvMeasurement.selectedPeakDelayStdSec{rvPulseIndex, rvSelectedSlot} = ...
            rvSelectedDelayStdSec{rvSelectedSlot};
    end

    rvDelay1 = (double(rvPeak1(:).') - 1) ./ fs;
    rvDelay2 = (double(rvPeak2(:).') - 1) ./ fs;
    rvDelay3 = (double(rvPeak3(:).') - 1) ./ fs;
    [rvRangeMatrix12, rvVelocityMatrix12] = radar.Cal_rv_mat( ...
        rvDelay1, rvDelay2, fc_v(rvSelectedSegments(1)), ...
        fc_v(rvSelectedSegments(2)), k_v(rvSelectedSegments(1)), ...
        k_v(rvSelectedSegments(2)));
    [rvRangeMatrix13, rvVelocityMatrix13] = radar.Cal_rv_mat( ...
        rvDelay1, rvDelay3, fc_v(rvSelectedSegments(1)), ...
        fc_v(rvSelectedSegments(3)), k_v(rvSelectedSegments(1)), ...
        k_v(rvSelectedSegments(3)));

    rvConsistencyCommonPeakCount = min([numel(rvPeak1), numel(rvPeak2), numel(rvPeak3)]);
    if strcmpi(rvConsistencyBiasCorrection, 'rankMedian') && rvConsistencyCommonPeakCount > 0
        rvConsistencyCommonIndices = 1:rvConsistencyCommonPeakCount;
        rvConsistencyRangeBiasM = median(diag(rvRangeMatrix12( ...
            rvConsistencyCommonIndices, rvConsistencyCommonIndices)) - ...
            diag(rvRangeMatrix13(rvConsistencyCommonIndices, rvConsistencyCommonIndices)));
        rvConsistencyVelocityBiasMps = median(diag(rvVelocityMatrix12( ...
            rvConsistencyCommonIndices, rvConsistencyCommonIndices)) - ...
            diag(rvVelocityMatrix13(rvConsistencyCommonIndices, rvConsistencyCommonIndices)));
    else
        rvConsistencyRangeBiasM      = 0;
        rvConsistencyVelocityBiasMps = 0;
    end
    rvMeasurement.consistencyRangeBiasM(rvPulseIndex)      = rvConsistencyRangeBiasM;
    rvMeasurement.consistencyVelocityBiasMps(rvPulseIndex) = rvConsistencyVelocityBiasMps;

    rvQ1 = k_v(rvSelectedSegments(1)) ./ fc_v(rvSelectedSegments(1));
    rvQ2 = k_v(rvSelectedSegments(2)) ./ fc_v(rvSelectedSegments(2));
    rvQ3 = k_v(rvSelectedSegments(3)) ./ fc_v(rvSelectedSegments(3));
    rvRangeGain12_1 = c ./ 2 .* rvQ1 ./ (rvQ1 - rvQ2);
    rvRangeGain12_2 = -c ./ 2 .* rvQ2 ./ (rvQ1 - rvQ2);
    rvRangeGain13_1 = c ./ 2 .* rvQ1 ./ (rvQ1 - rvQ3);
    rvRangeGain13_3 = -c ./ 2 .* rvQ3 ./ (rvQ1 - rvQ3);
    rvVelocityGain12 = c ./ 2 ./ ( ...
        fc_v(rvSelectedSegments(1)) ./ k_v(rvSelectedSegments(1)) - ...
        fc_v(rvSelectedSegments(2)) ./ k_v(rvSelectedSegments(2)));
    rvVelocityGain13 = c ./ 2 ./ ( ...
        fc_v(rvSelectedSegments(1)) ./ k_v(rvSelectedSegments(1)) - ...
        fc_v(rvSelectedSegments(3)) ./ k_v(rvSelectedSegments(3)));
    rvRangeDifferenceJacobian = [ ...
        rvRangeGain12_1 - rvRangeGain13_1, rvRangeGain12_2, -rvRangeGain13_3];
    rvVelocityDifferenceJacobian = [ ...
        rvVelocityGain12 - rvVelocityGain13, -rvVelocityGain12, rvVelocityGain13];

    rvMatchMask                = false(size(rvRangeMatrix12));
    rvMatchThirdPeakIndex      = nan(size(rvRangeMatrix12));
    rvMatchRangeResidualM      = nan(size(rvRangeMatrix12));
    rvMatchVelocityResidualMps = nan(size(rvRangeMatrix12));
    rvMatchRangeThresholdM     = nan(size(rvRangeMatrix12));
    rvMatchVelocityThresholdMps = nan(size(rvRangeMatrix12));
    rvMatchNormalizedScore     = nan(size(rvRangeMatrix12));

    for rvCommonPeakIndex = 1:size(rvRangeMatrix12, 1)
        for rvSecondPeakIndex = 1:size(rvRangeMatrix12, 2)
            if abs(rvSecondPeakIndex - rvCommonPeakIndex) > rvConsistencyMaxPeakRankDifference
                continue;
            end
            rvRangeDifferenceM = abs(rvRangeMatrix12(rvCommonPeakIndex, rvSecondPeakIndex) - ...
                rvRangeMatrix13(rvCommonPeakIndex, :) - rvConsistencyRangeBiasM);
            rvVelocityDifferenceMps = abs(rvVelocityMatrix12(rvCommonPeakIndex, rvSecondPeakIndex) - ...
                rvVelocityMatrix13(rvCommonPeakIndex, :) - rvConsistencyVelocityBiasMps);

            if strcmpi(rvConsistencyMode, 'snrPrecision')
                rvDelayStd1 = rvSelectedDelayStdSec{1}(rvCommonPeakIndex);
                rvDelayStd2 = rvSelectedDelayStdSec{2}(rvSecondPeakIndex);
                rvDelayStd3 = rvSelectedDelayStdSec{3};
                rvRangeDifferenceStdM = sqrt( ...
                    (rvRangeDifferenceJacobian(1) .* rvDelayStd1).^2 + ...
                    (rvRangeDifferenceJacobian(2) .* rvDelayStd2).^2 + ...
                    (rvRangeDifferenceJacobian(3) .* rvDelayStd3).^2);
                rvVelocityDifferenceStdMps = sqrt( ...
                    (rvVelocityDifferenceJacobian(1) .* rvDelayStd1).^2 + ...
                    (rvVelocityDifferenceJacobian(2) .* rvDelayStd2).^2 + ...
                    (rvVelocityDifferenceJacobian(3) .* rvDelayStd3).^2);
                rvRangeGateM      = max(rvConsistencySigmaMultiplier .* rvRangeDifferenceStdM, eps);
                rvVelocityGateMps = max(rvConsistencySigmaMultiplier .* rvVelocityDifferenceStdMps, eps);
            else
                rvRangeGateM      = rvLegacyRangeConsistencyM .* ones(size(rvRangeDifferenceM));
                rvVelocityGateMps = rvLegacyVelocityConsistencyMps .* ones(size(rvVelocityDifferenceMps));
            end

            rvCompatibleMask = rvRangeDifferenceM < rvRangeGateM & ...
                rvVelocityDifferenceMps < rvVelocityGateMps;
            rvCompatibleMask = rvCompatibleMask & ...
                abs((1:numel(rvCompatibleMask)) - rvCommonPeakIndex) <= ...
                rvConsistencyMaxPeakRankDifference;

            if any(rvCompatibleMask)
                rvNormalizedMismatch = sqrt( ...
                    (rvRangeDifferenceM ./ rvRangeGateM).^2 + ...
                    (rvVelocityDifferenceMps ./ rvVelocityGateMps).^2);
                rvNormalizedMismatch(~rvCompatibleMask) = inf;
                [rvBestMismatch, rvBestThirdPeakIndex] = min(rvNormalizedMismatch);
                if isfinite(rvBestMismatch)
                    rvMatchMask(rvCommonPeakIndex, rvSecondPeakIndex) = true;
                    rvMatchThirdPeakIndex(rvCommonPeakIndex, rvSecondPeakIndex) = ...
                        rvBestThirdPeakIndex;
                    rvMatchRangeResidualM(rvCommonPeakIndex, rvSecondPeakIndex) = ...
                        rvRangeDifferenceM(rvBestThirdPeakIndex);
                    rvMatchVelocityResidualMps(rvCommonPeakIndex, rvSecondPeakIndex) = ...
                        rvVelocityDifferenceMps(rvBestThirdPeakIndex);
                    rvMatchRangeThresholdM(rvCommonPeakIndex, rvSecondPeakIndex) = ...
                        rvRangeGateM(rvBestThirdPeakIndex);
                    rvMatchVelocityThresholdMps(rvCommonPeakIndex, rvSecondPeakIndex) = ...
                        rvVelocityGateMps(rvBestThirdPeakIndex);
                    rvMatchNormalizedScore(rvCommonPeakIndex, rvSecondPeakIndex) = ...
                        rvBestMismatch;
                end
            end
        end
    end

    rvMeasurement.preSpeedGateCandidateCount(rvPulseIndex) = nnz(rvMatchMask);
    rvMatchMask = rvMatchMask & abs(rvVelocityMatrix12) <= rvMaxAbsVelocityMps;
    rvMeasurement.postSpeedGateCandidateCount(rvPulseIndex) = nnz(rvMatchMask);

    if rvConsistencyEnforceOneToOne && any(rvMatchMask(:))
        [rvOneToOneRows, rvOneToOneCols] = find(rvMatchMask);
        rvOneToOneLinear = sub2ind(size(rvMatchMask), rvOneToOneRows, rvOneToOneCols);
        [~, rvOneToOneOrder] = sort(rvMatchNormalizedScore(rvOneToOneLinear), 'ascend');
        rvOneToOneMask = false(size(rvMatchMask));
        rvUsedPeak1 = false(size(rvRangeMatrix12, 1), 1);
        rvUsedPeak2 = false(size(rvRangeMatrix12, 2), 1);
        rvUsedPeak3 = false(size(rvRangeMatrix13, 2), 1);
        for rvOneToOneCandidateIndex = rvOneToOneOrder(:).'
            rvOneToOneRow   = rvOneToOneRows(rvOneToOneCandidateIndex);
            rvOneToOneCol   = rvOneToOneCols(rvOneToOneCandidateIndex);
            rvOneToOneThird = round(rvMatchThirdPeakIndex(rvOneToOneRow, rvOneToOneCol));
            if ~rvUsedPeak1(rvOneToOneRow) && ...
                    ~rvUsedPeak2(rvOneToOneCol) && ...
                    ~rvUsedPeak3(rvOneToOneThird)
                rvOneToOneMask(rvOneToOneRow, rvOneToOneCol) = true;
                rvUsedPeak1(rvOneToOneRow)   = true;
                rvUsedPeak2(rvOneToOneCol)   = true;
                rvUsedPeak3(rvOneToOneThird) = true;
            end
        end
        rvMatchMask = rvOneToOneMask;
    end

    rvMeasurement.postOneToOneCandidateCount(rvPulseIndex) = nnz(rvMatchMask);
    rvMeasurement.matchThirdPeakIndex{rvPulseIndex}        = rvMatchThirdPeakIndex;
    rvMeasurement.matchRangeResidualM{rvPulseIndex}        = rvMatchRangeResidualM;
    rvMeasurement.matchVelocityResidualMps{rvPulseIndex}   = rvMatchVelocityResidualMps;
    rvMeasurement.matchRangeThresholdM{rvPulseIndex}       = rvMatchRangeThresholdM;
    rvMeasurement.matchVelocityThresholdMps{rvPulseIndex}  = rvMatchVelocityThresholdMps;
    rvMeasurement.matchNormalizedScore{rvPulseIndex}       = rvMatchNormalizedScore;

    rvCandidateRangeM      = rvRangeMatrix12(rvMatchMask);
    rvCandidateVelocityMps = rvVelocityMatrix12(rvMatchMask);
    [rvCandidateRangeM, rvCandidateOrder] = sort(rvCandidateRangeM(:).');
    rvCandidateVelocityMps = rvCandidateVelocityMps(:).';
    rvCandidateVelocityMps = rvCandidateVelocityMps(rvCandidateOrder);
    rvMeasurement.rangeCandidatesM{rvPulseIndex}      = rvCandidateRangeM;
    rvMeasurement.velocityCandidatesMps{rvPulseIndex} = rvCandidateVelocityMps;

    [~, rvTruthOrder] = sort(rvTruthRangeM(rvPulseIndex, :));
    rvPairCount = min(numel(rvCandidateRangeM), rvTargetCount);
    if rvPairCount > 0
        rvPairedTargetColumns = rvTruthOrder(1:rvPairCount);
        rvMeasurement.rangeM(rvPulseIndex, rvPairedTargetColumns) = ...
            rvCandidateRangeM(1:rvPairCount);
        rvMeasurement.velocityMps(rvPulseIndex, rvPairedTargetColumns) = ...
            rvCandidateVelocityMps(1:rvPairCount);
    end
    rvMeasurement.matchedCount(rvPulseIndex)       = rvPairCount;
    rvMeasurement.missedTargetCount(rvPulseIndex)  = rvTargetCount - rvPairCount;
    rvMeasurement.falseEstimateCount(rvPulseIndex) = max(numel(rvCandidateRangeM) - rvTargetCount, 0);
end

%% 4. 第 1 脉冲常量均值偏置校准与误差统计
rvCalibrationMask = isfinite(rvMeasurement.rangeM(1, :)) & ...
    isfinite(rvMeasurement.velocityMps(1, :));
rvMeasurement.calibrationSampleCount = nnz(rvCalibrationMask);
rvMeasurement.calibrationAvailable   = any(rvCalibrationMask);
if rvMeasurement.calibrationAvailable
    rvMeasurement.rangeMeanBiasM = mean( ...
        rvMeasurement.rangeM(1, rvCalibrationMask) - ...
        rvTruthRangeM(1, rvCalibrationMask));
    rvMeasurement.velocityMeanBiasMps = mean( ...
        rvMeasurement.velocityMps(1, rvCalibrationMask) - ...
        rvTruthVelocityMps(1, rvCalibrationMask));
else
    rvMeasurement.rangeMeanBiasM      = 0;
    rvMeasurement.velocityMeanBiasMps = 0;
end

rvMeasurement.calibratedRangeM      = rvMeasurement.rangeM - rvMeasurement.rangeMeanBiasM;
rvMeasurement.calibratedVelocityMps = rvMeasurement.velocityMps - rvMeasurement.velocityMeanBiasMps;

rvMeasurement.rawRangeErrorM             = rvMeasurement.rangeM - rvTruthRangeM;
rvMeasurement.rawVelocityErrorMps        = rvMeasurement.velocityMps - rvTruthVelocityMps;
rvMeasurement.calibratedRangeErrorM      = rvMeasurement.calibratedRangeM - rvTruthRangeM;
rvMeasurement.calibratedVelocityErrorMps = rvMeasurement.calibratedVelocityMps - rvTruthVelocityMps;

rvRangeReferenceM      = max(abs(rvTruthRangeM), eps);
rvVelocityReferenceMps = max(abs(rvTruthVelocityMps), eps);
rvMeasurement.rawRangeRelativeError             = abs(rvMeasurement.rawRangeErrorM) ./ rvRangeReferenceM;
rvMeasurement.rawVelocityRelativeError          = abs(rvMeasurement.rawVelocityErrorMps) ./ rvVelocityReferenceMps;
rvMeasurement.calibratedRangeRelativeError      = abs(rvMeasurement.calibratedRangeErrorM) ./ rvRangeReferenceM;
rvMeasurement.calibratedVelocityRelativeError   = abs(rvMeasurement.calibratedVelocityErrorMps) ./ rvVelocityReferenceMps;

rvMeasurement.maxRawRangeRelativeErrorPercent = max_percent_finite( ...
    rvMeasurement.rawRangeRelativeError);
rvMeasurement.maxRawVelocityRelativeErrorPercent = max_percent_finite( ...
    rvMeasurement.rawVelocityRelativeError);
rvMeasurement.maxCalibratedRangeRelativeErrorPercent = max_percent_finite( ...
    rvMeasurement.calibratedRangeRelativeError);
rvMeasurement.maxCalibratedVelocityRelativeErrorPercent = max_percent_finite( ...
    rvMeasurement.calibratedVelocityRelativeError);

if printSummary
    fprintf('\n================ 六段 LFM 距离/速度解算 ================\n');
    fprintf('R/V consistency gate: %s', rvConsistencyMode);
    if strcmpi(rvConsistencyMode, 'snrPrecision')
        fprintf([', %.2f sigma, SNR [%.1f, %.1f] dB, delay floor %.3f bin, ' ...
            'noise %s, bias %s, rank %g, one-to-one %d\n'], ...
            rvConsistencySigmaMultiplier, rvConsistencySnrFloorDb, ...
            rvConsistencySnrCeilingDb, rvConsistencyDelayStdFloorBins, ...
            rvConsistencyNoiseEstimator, rvConsistencyBiasCorrection, ...
            rvConsistencyMaxPeakRankDifference, rvConsistencyEnforceOneToOne);
    else
        fprintf(', fixed range %.6f m, velocity %.6f m/s\n', ...
            rvLegacyRangeConsistencyM, rvLegacyVelocityConsistencyMps);
    end
    fprintf('距离 mean 偏置：%.6f m；速度 mean 偏置：%.6f m/s；校准样本数：%d\n', ...
        rvMeasurement.rangeMeanBiasM, rvMeasurement.velocityMeanBiasMps, ...
        rvMeasurement.calibrationSampleCount);
    fprintf('原始测距最大相对误差：%.6f%%；校准后：%.6f%%\n', ...
        rvMeasurement.maxRawRangeRelativeErrorPercent, ...
        rvMeasurement.maxCalibratedRangeRelativeErrorPercent);
    fprintf('原始测速最大相对误差：%.6f%%；校准后：%.6f%%\n', ...
        rvMeasurement.maxRawVelocityRelativeErrorPercent, ...
        rvMeasurement.maxCalibratedVelocityRelativeErrorPercent);
    fprintf('总匹配数：%d；总漏检数：%d；总额外估计数：%d\n', ...
        sum(rvMeasurement.matchedCount), sum(rvMeasurement.missedTargetCount), ...
        sum(rvMeasurement.falseEstimateCount));
end
end

function val = get_field_or_default(s, name, defaultVal)
if isfield(s, name) && ~isempty(s.(name))
    val = s.(name);
else
    val = defaultVal;
end
end

function pct = max_percent_finite(relErr)
valid = isfinite(relErr);
if any(valid(:))
    pct = 100 .* max(relErr(valid));
else
    pct = NaN;
end
end
