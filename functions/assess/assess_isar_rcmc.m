function [detail, summaryTbl, savedPath] = assess_isar_rcmc(isarResOrPost, assessCfg, saveOpts)
%ASSESS_ISAR_RCMC 独立评估 ISAR RCMC 包络对齐前后的脉冲间抖动改善效果。
%
% 功能说明：
%   复用 lib 中的 radar.assess_hrrf_alignment，对比逐目标 RCMC 对齐前的
%   hrrfModelCrop 与对齐后的 hrrfAligned，输出对齐前后距离峰值位置标准差（bins）
%   及最大偏移统计，并支持按 saveOpts 导出评估表。
%
% 输入参数：
%   isarResOrPost - 1×1 struct，isarpath_run 输出 isarRes（含 .post）或直接为 isarPost
%   assessCfg     - (可选) 1×1 struct，评估配置项：
%     .targetIndices     1×K double，待评估目标索引（默认全部目标）
%     .useSubbinPeakFit  1×1 logical，是否启用子像素峰拟合（默认按 rangeAlignmentMode 判断）
%     .printSummary      1×1 logical，是否打印汇总表（默认 false）
%   saveOpts      - (可选) 1×1 struct，保存配置（enable/dir/fileName/format，默认不保存）
%
% 输出参数：
%   detail     - K×1 struct，radar.assess_hrrf_alignment 返回的逐目标完整评估结构体
%   summaryTbl - K×4 table，逐目标对齐前/后标准差（bins）与最大 RCMC 位移汇总表
%   savedPath  - 1×1 string，表格导出路径（未保存时为 ""）

if nargin < 2 || isempty(assessCfg), assessCfg = struct(); end
if nargin < 3, saveOpts = struct(); end

if isfield(isarResOrPost, 'post'), post = isarResOrPost.post; else, post = isarResOrPost; end

nTotal = numel(post.hrrfModelCrop);
if isfield(assessCfg, 'targetIndices') && ~isempty(assessCfg.targetIndices)
    idx = double(assessCfg.targetIndices(:));
else
    idx = (1:nTotal).';
end
nEval = numel(idx);

useSubbin = false;
if isfield(post, 'parameters') && isfield(post.parameters, 'rangeAlignmentMode')
    useSubbin = (upper(string(post.parameters.rangeAlignmentMode)) == "FRACTIONAL");
end
if isfield(assessCfg, 'useSubbinPeakFit') && ~isempty(assessCfg.useSubbinPeakFit)
    useSubbin = logical(assessCfg.useSubbinPeakFit);
end
printSummary = isfield(assessCfg, 'printSummary') && logical(assessCfg.printSummary);

alignOpts = struct('UseSubbinPeakFit', useSubbin, 'AxisUnit', 'm');

detailCells   = cell(nEval, 1);
stdBeforeBins = nan(nEval, 1);
stdAfterBins  = nan(nEval, 1);
maxShiftBins  = nan(nEval, 1);

for k = 1:nEval
    tIdx = idx(k);
    q = radar.assess_hrrf_alignment( ...
        post.hrrfModelCrop{tIdx}, post.hrrfAligned{tIdx}, post.rangeAxisM{tIdx}, alignOpts);
    q.targetIndex = tIdx;
    detailCells{k} = q;
    stdBeforeBins(k) = q.summary.beforeStdBins;
    stdAfterBins(k)  = q.summary.afterStdBins;
    if isfield(post, 'rangeMigrationBins') && size(post.rangeMigrationBins, 1) >= tIdx
        maxShiftBins(k) = max(abs(post.rangeMigrationBins(tIdx, :)));
    end
end
detail = vertcat(detailCells{:});

summaryTbl = table(idx, stdBeforeBins, stdAfterBins, maxShiftBins, ...
    'VariableNames', { ...
    'Target', 'AlignmentStdBeforeBins', 'AlignmentStdAfterBins', 'MaxAbsRcmcShiftBins'});

if printSummary
    fprintf('\n=== ISAR RCMC 包络对齐质量评估 ===\n');
    disp(summaryTbl);
end

savedPath = save_output_artifact(summaryTbl, "table", saveOpts, 'assess_isar_rcmc');
end
