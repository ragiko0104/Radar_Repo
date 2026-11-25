function [detail, summaryTbl, savedPath] = assess_isar_scatterers(isarResOrPost, assessCfg, saveOpts)
%ASSESS_ISAR_SCATTERERS 独立评估 ISAR 图像中强散射点检出数及二维 3dB 分辨率。
%
% 功能说明：
%   复用 lib 中的 radar.assess_isar_scatterer_resolution，对逐目标二维 ISAR 幅度像
%   提取期望数量的强亮点并统计其距离向/方位向 3dB 分辨率均值，支持按 saveOpts 导出汇总表。
%
% 输入参数：
%   isarResOrPost - 1×1 struct，isarpath_run 输出 isarRes（含 .post）或直接为 isarPost
%   assessCfg     - (可选) 1×1 struct，评估配置项：
%     .targetIndices      1×K double，待评估目标索引（默认全部目标）
%     .expectedPeakCount  1×1 或 1×K double，期望强亮点数（默认 4）
%     .resolutionMethod   char/string，3dB 宽度计算方法（默认 "linear"）
%     .resolutionConfig   1×1 struct，透传给 ResolutionConfig
%     .printSummary       1×1 logical，是否打印汇总表（默认 false）
%   saveOpts      - (可选) 1×1 struct，保存配置（enable/dir/fileName/format，默认不保存）
%
% 输出参数：
%   detail     - K×1 struct，radar.assess_isar_scatterer_resolution 返回的逐目标明细
%   summaryTbl - K×5 table，逐目标期望/实测亮点数及距离/方位 3dB 均值分辨率汇总表
%   savedPath  - 1×1 string，表格导出路径（未保存时为 ""）

if nargin < 2 || isempty(assessCfg), assessCfg = struct(); end
if nargin < 3, saveOpts = struct(); end

if isfield(isarResOrPost, 'post'), post = isarResOrPost.post; else, post = isarResOrPost; end

nTotal = numel(post.image);
if isfield(assessCfg, 'targetIndices') && ~isempty(assessCfg.targetIndices)
    idx = double(assessCfg.targetIndices(:));
else
    idx = (1:nTotal).';
end
nEval = numel(idx);

expectedPeaks = 4;
if isfield(assessCfg, 'expectedPeakCount') && ~isempty(assessCfg.expectedPeakCount)
    expectedPeaks = double(assessCfg.expectedPeakCount(:));
end
if isscalar(expectedPeaks)
    expectedPeaks = repmat(expectedPeaks, nEval, 1);
elseif numel(expectedPeaks) == nTotal
    expectedPeaks = expectedPeaks(idx);
elseif numel(expectedPeaks) ~= nEval
    error('RV_ISAR_Path:ExpectedPeakCountMismatch', ...
        'expectedPeakCount 必须为标量或与评估目标数(%d)一致。', nEval);
end

resMethod = "linear";
if isfield(assessCfg, 'resolutionMethod') && ~isempty(assessCfg.resolutionMethod)
    resMethod = string(assessCfg.resolutionMethod);
end
resCfg = struct();
if isfield(assessCfg, 'resolutionConfig') && isstruct(assessCfg.resolutionConfig)
    resCfg = assessCfg.resolutionConfig;
end
printSummary = isfield(assessCfg, 'printSummary') && logical(assessCfg.printSummary);

scatOpts = struct( ...
    'ResolutionMethod', resMethod, ...
    'ResolutionConfig', resCfg, ...
    'RangeUnit',        'm', ...
    'CrossRangeUnit',   'm');

detailCells   = cell(nEval, 1);
detectedCount = zeros(nEval, 1);
rangeResMean  = nan(nEval, 1);
crossResMean  = nan(nEval, 1);

for k = 1:nEval
    tIdx = idx(k);
    q = radar.assess_isar_scatterer_resolution( ...
        post.image{tIdx}, post.rangeAxisM{tIdx}, post.crossAxisM{tIdx}, ...
        expectedPeaks(k), scatOpts);
    q.targetIndex    = tIdx;
    detailCells{k}   = q;
    detectedCount(k) = q.detectedPeakCount;
    rangeResMean(k)  = q.summary.range3dBMean;
    crossResMean(k)  = q.summary.crossRange3dBMean;
end
detail = vertcat(detailCells{:});

summaryTbl = table(idx, expectedPeaks, detectedCount, rangeResMean, crossResMean, ...
    'VariableNames', { ...
    'Target', 'ExpectedPeakCount', 'DetectedPeakCount', ...
    'ScattererRange3dBMeanM', 'ScattererCrossRange3dBMeanM'});

if printSummary
    fprintf('\n=== ISAR 散射点二维分辨率评估 ===\n');
    disp(summaryTbl);
end

savedPath = save_output_artifact(summaryTbl, "table", saveOpts, 'assess_isar_scatterers');
end
