function [detail, summaryTbl, savedPath] = assess_isar_dispersion(isarResOrPost, assessCfg, saveOpts)
%ASSESS_ISAR_DISPERSION 独立评估 ISAR 脉内多普勒色散补偿前后的 HRRP 展宽与 PSLR。
%
% 功能说明：
%   复用 lib 中的 radar.assess_hrrf_quality，对逐目标色散补偿前的 HRRP
%   （hrrfBeforeDispersion）与补偿后的 HRRP（hrrfModelCrop）分别计算 3dB 主瓣宽度与
%   PSLR，量化展宽改善比例，并支持按 saveOpts 导出评估表。
%
% 输入参数：
%   isarResOrPost - 1×1 struct，isarpath_run 输出 isarRes（含 .post）或直接为 isarPost
%   assessCfg     - (可选) 1×1 struct，评估配置项：
%     .targetIndices     1×K double，待评估目标索引（默认全部目标）
%     .resolutionMethod  char/string，3dB 宽度计算方法（默认 "linear"）
%     .resolutionConfig  1×1 struct，透传给 radar.assess_hrrf_quality.ResolutionConfig
%     .printSummary      1×1 logical，是否打印汇总表（默认 false）
%   saveOpts      - (可选) 1×1 struct，保存配置（enable/dir/fileName/format，默认不保存）
%
% 输出参数：
%   detail     - K×1 struct，逐目标明细（beforeDispersion/afterDispersion/dispersionComparison）
%   summaryTbl - K×7 table，逐目标色散补偿前后 3dB 宽度、展宽缩减百分比及 PSLR 对比表
%   savedPath  - 1×1 string，表格导出路径（未保存时为 ""）

if nargin < 2 || isempty(assessCfg), assessCfg = struct(); end
if nargin < 3, saveOpts = struct(); end

post = unwrap_isar_post(isarResOrPost);
assert(isfield(post, 'hrrfBeforeDispersion') && ~isempty(post.hrrfBeforeDispersion), ...
    'isarPost 缺少 hrrfBeforeDispersion，无法对比色散补偿前后 HRRP。');

nTotal = numel(post.hrrfModelCrop);
idx = resolve_target_indices(assessCfg, nTotal);
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

hrrfOpts = struct( ...
    'ResolutionMethod', resMethod, ...
    'ResolutionConfig', resCfg, ...
    'AxisUnit',         'm');

applied = true;
if isfield(post, 'parameters') && isfield(post.parameters, 'dispersionCompensationEnabled')
    applied = logical(post.parameters.dispersionCompensationEnabled);
end

detail = repmat(struct( ...
    'targetIndex',          0, ...
    'beforeDispersion',     struct(), ...
    'afterDispersion',      struct(), ...
    'dispersionComparison', struct()), nEval, 1);

wBefore = nan(nEval, 1);
wAfter  = nan(nEval, 1);
wRedPct = nan(nEval, 1);
pBefore = nan(nEval, 1);
pAfter  = nan(nEval, 1);
pChange = nan(nEval, 1);

for k = 1:nEval
    tIdx = idx(k);
    rangeAxis = post.rangeAxisM{tIdx};
    qBefore = radar.assess_hrrf_quality(post.hrrfBeforeDispersion{tIdx}, rangeAxis, hrrfOpts);
    qAfter  = radar.assess_hrrf_quality(post.hrrfModelCrop{tIdx}, rangeAxis, hrrfOpts);

    wb = qBefore.summary.range3dBMean;
    wa = qAfter.summary.range3dBMean;
    wRed = wb - wa;
    if isfinite(wb) && isfinite(wa) && wb > 0
        pct = 100 * wRed / wb;
    elseif wb == 0 && wa == 0
        pct = 0;
    else
        pct = NaN;
    end
    pb = qBefore.summary.pslrMeanDb;
    pa = qAfter.summary.pslrMeanDb;

    comp = struct( ...
        'applied',                   applied, ...
        'meanWidthBeforeM',          wb, ...
        'meanWidthAfterM',           wa, ...
        'meanWidthReductionM',       wRed, ...
        'meanWidthReductionPercent', pct, ...
        'meanPslrBeforeDb',          pb, ...
        'meanPslrAfterDb',           pa, ...
        'meanPslrChangeDb',          pa - pb);

    detail(k).targetIndex          = tIdx;
    detail(k).beforeDispersion     = qBefore;
    detail(k).afterDispersion      = qAfter;
    detail(k).dispersionComparison = comp;

    wBefore(k) = wb;
    wAfter(k)  = wa;
    wRedPct(k) = pct;
    pBefore(k) = pb;
    pAfter(k)  = pa;
    pChange(k) = pa - pb;
end

summaryTbl = table(idx, wBefore, wAfter, wRedPct, pBefore, pAfter, pChange, ...
    'VariableNames', { ...
    'Target', 'Hrrf3dBBeforeM', 'Hrrf3dBAfterM', ...
    'DispersionWidthReductionPercent', 'HrrfPslrBeforeDb', ...
    'HrrfPslrAfterDb', 'HrrfPslrChangeDb'});

if printSummary
    fprintf('\n=== ISAR 色散补偿 HRRP 质量评估 ===\n');
    disp(summaryTbl);
end

savedPath = save_output_artifact(summaryTbl, "table", saveOpts, 'assess_isar_dispersion');
end

function post = unwrap_isar_post(x)
if isfield(x, 'post'), post = x.post; else, post = x; end
end

function idx = resolve_target_indices(cfg, nTotal)
if isfield(cfg, 'targetIndices') && ~isempty(cfg.targetIndices)
    idx = double(cfg.targetIndices(:));
else
    idx = (1:nTotal).';
end
end
