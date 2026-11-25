function [detail, summaryTbl, savedPath] = assess_isar_image(isarResOrPost, assessCfg, saveOpts)
%ASSESS_ISAR_IMAGE 独立评估 ISAR 二维图像的熵、归一化熵、2D PSLR 与 3dB 能量比。
%
% 功能说明：
%   复用 lib 中的 radar.assess_isar_image_quality，对逐目标二维 ISAR 幅度像计算
%   图像熵、归一化熵、二维峰值旁瓣比（PSLR）与 3dB 主瓣能量集中度，并支持按
%   saveOpts 导出汇总表。
%
% 输入参数：
%   isarResOrPost - 1×1 struct，isarpath_run 输出 isarRes（含 .post）或直接为 isarPost
%   assessCfg     - (可选) 1×1 struct，评估配置项：
%     .targetIndices     1×K double，待评估目标索引（默认全部目标）
%     .resolutionMethod  char/string，3dB 评估方法（默认 "linear"）
%     .resolutionConfig  1×1 struct，透传给 radar.assess_isar_image_quality.ResolutionConfig
%     .printSummary      1×1 logical，是否打印汇总表（默认 false）
%   saveOpts      - (可选) 1×1 struct，保存配置（enable/dir/fileName/format，默认不保存）
%
% 输出参数：
%   detail     - K×1 struct，radar.assess_isar_image_quality 返回的逐目标完整评估结构体
%   summaryTbl - K×5 table，逐目标图像熵/归一化熵/2D PSLR/3dB 能量比汇总表
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

resMethod = "linear";
if isfield(assessCfg, 'resolutionMethod') && ~isempty(assessCfg.resolutionMethod)
    resMethod = string(assessCfg.resolutionMethod);
end
resCfg = struct();
if isfield(assessCfg, 'resolutionConfig') && isstruct(assessCfg.resolutionConfig)
    resCfg = assessCfg.resolutionConfig;
end
printSummary = isfield(assessCfg, 'printSummary') && logical(assessCfg.printSummary);

imgOpts = struct( ...
    'ResolutionMethod', resMethod, ...
    'ResolutionConfig', resCfg, ...
    'RangeUnit',        'm', ...
    'CrossRangeUnit',   'm');

detailCells = cell(nEval, 1);
entropyVal  = nan(nEval, 1);
normEntVal  = nan(nEval, 1);
pslr2dVal   = nan(nEval, 1);
energy3dB   = nan(nEval, 1);

for k = 1:nEval
    tIdx = idx(k);
    q = radar.assess_isar_image_quality( ...
        post.image{tIdx}, post.rangeAxisM{tIdx}, post.crossAxisM{tIdx}, imgOpts);
    q.targetIndex  = tIdx;
    detailCells{k} = q;
    entropyVal(k)  = q.entropy;
    normEntVal(k)  = q.normalizedEntropy;
    pslr2dVal(k)   = q.pslr2dDb;
    energy3dB(k)   = q.energyRatio3dB;
end
detail = vertcat(detailCells{:});

summaryTbl = table(idx, entropyVal, normEntVal, pslr2dVal, energy3dB, ...
    'VariableNames', { ...
    'Target', 'ImageEntropy', 'ImageEntropyNormalized', ...
    'ImagePslr2dDb', 'ImageEnergyRatio3dB'});

if printSummary
    fprintf('\n=== ISAR 图像聚焦质量评估 ===\n');
    disp(summaryTbl);
end

savedPath = save_output_artifact(summaryTbl, "table", saveOpts, 'assess_isar_image');
end
