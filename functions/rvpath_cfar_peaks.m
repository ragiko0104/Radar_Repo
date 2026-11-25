function rvCfar = rvpath_cfar_peaks(rvPulseCompression, wf, cfarOverrides)
%RVPATH_CFAR_PEAKS 对六段脉压功率执行 MEAN-CFAR 检测与精确子像素寻峰。
%
% 功能说明：
%   1. 根据默认（或覆盖的）CFAR 参数计算门限因子 alpha（radar.cfarth_alpha）
%      以及各频带最小峰间隔 peakSeparationBinsByBand；
%   2. 在六段脉压功率 abs(compressed{k}).^2 上依次调用 radar.tm_cfar 与
%      radar.find_Peak_exact，获得每段每脉冲的 CFAR 门限、候选峰、精确峰位及
%      top-k 峰均值二次筛选阈值。
%
% 输入参数：
%   rvPulseCompression - 1×1 struct，由 rvpath_pulse_compress 产出的脉压结果
%                        （使用 segmentKeys、segmentNames、compressed{6×1}）
%   wf                 - 1×1 struct，由 rvpath_derive_waveform 产出的派生 RV 波形结构体
%                        （使用 fs、B1、B2、B3、Na）
%   cfarOverrides      - (可选) 1×1 struct，用于覆盖默认 CFAR 参数中的部分字段
%
% 输出参数：
%   rvCfar             - 1×1 struct，六段 CFAR 检测与寻峰结果：
%     .segmentKeys       6×1 string，六段键名
%     .segmentNames      6×1 string，六段中文名
%     .segmentBandIndex  6×1 double，各段所属频带编号 [1;1;2;2;3;3]
%     .parameters        1×1 struct，生效的 CFAR 与寻峰参数（含 alpha、peakSeparationBinsByBand）
%     .threshold         6×1 cell，每元素为 Na×N double，CFAR 检测门限曲线
%     .candidatePeaks    6×1 cell，每元素为 Na×1 cell（或向量），CFAR 初检候选峰索引
%     .refinedPeaks      6×1 cell，每元素为 Na×1 cell，精确寻峰后的 1-based 峰索引
%     .meanPeakPower     Na×6 double，每脉冲每段的平均峰功率（用于选取能量最高三段）
%     .topKMeanThreshold Na×6 double，每脉冲每段的 top-k 峰均值二次门限

params = struct( ...
    'trainingCellsPerSide',       64, ...
    'guardCellsByBand',           [32, 32, 32], ...
    'trimLargestPerSide',         0, ...
    'falseAlarmProbability',      1e-4, ...
    'alphaDomain',                'P', ...
    'cfarMethod',                 'MEAN', ...
    'consecutiveDetectionLength', 0, ...
    'orderStatisticRank',         round(2 * 64 * 0.8), ...
    'hardThresholdMode',          'NO', ...
    'hardThresholdPowerByBand',   [0, 0, 0], ...
    'peakAverageMode',            'YES', ...
    'topKPeakCount',              16, ...
    'topKScale',                  0.5);

if nargin >= 3 && ~isempty(cfarOverrides)
    if ~(isstruct(cfarOverrides) && isscalar(cfarOverrides))
        error('RV_ISAR_Path:InvalidCfarOverrides', 'cfarOverrides 必须是 1×1 结构体。');
    end
    overrideNames = fieldnames(cfarOverrides);
    for k = 1:numel(overrideNames)
        params.(overrideNames{k}) = cfarOverrides.(overrideNames{k});
    end
end

params.alpha = radar.cfarth_alpha( ...
    params.trainingCellsPerSide, ...
    params.falseAlarmProbability, ...
    params.alphaDomain);
params.peakSeparationBinsByBand = round(wf.fs ./ abs([wf.B1, wf.B2, wf.B3]));

nSeg = numel(rvPulseCompression.segmentKeys);
Na   = wf.Na;

rvCfar = struct();
rvCfar.segmentKeys       = rvPulseCompression.segmentKeys;
rvCfar.segmentNames      = rvPulseCompression.segmentNames;
rvCfar.segmentBandIndex  = [1; 1; 2; 2; 3; 3];
rvCfar.parameters        = params;
rvCfar.threshold         = cell(nSeg, 1);
rvCfar.candidatePeaks    = cell(nSeg, 1);
rvCfar.refinedPeaks      = cell(nSeg, 1);
rvCfar.meanPeakPower     = zeros(Na, nSeg);
rvCfar.topKMeanThreshold = zeros(Na, nSeg);

for segIdx = 1:nSeg
    bandIdx = rvCfar.segmentBandIndex(segIdx);
    cfarPower = abs(rvPulseCompression.compressed{segIdx}).^2;

    [candidatePeaks, threshold] = radar.tm_cfar( ...
        cfarPower, params.trainingCellsPerSide, ...
        params.guardCellsByBand(bandIdx), ...
        params.trimLargestPerSide, params.alpha, ...
        params.cfarMethod, params.consecutiveDetectionLength, ...
        params.orderStatisticRank, params.hardThresholdMode, ...
        params.hardThresholdPowerByBand(bandIdx));

    [refinedPeaks, meanPeakPower, topKMeanThreshold] = radar.find_Peak_exact( ...
        cfarPower, candidatePeaks, ...
        params.peakSeparationBinsByBand(bandIdx), ...
        params.peakAverageMode, params.topKPeakCount, params.topKScale);

    rvCfar.threshold{segIdx}         = threshold;
    rvCfar.candidatePeaks{segIdx}    = candidatePeaks;
    rvCfar.refinedPeaks{segIdx}      = refinedPeaks;
    rvCfar.meanPeakPower(:, segIdx)  = meanPeakPower;
    rvCfar.topKMeanThreshold(:, segIdx) = topKMeanThreshold;
end
end
