function [fig, savedPath] = plot_isar_images(isarResOrPost, plotCfg, saveOpts)
%PLOT_ISAR_IMAGES 复用 lib 中 plotting.draw_isar_img 绘制逐目标 ISAR 二维成像结果。
%
% 功能说明：
%   从 isarRes（或 isarPost）中抽取选定目标的归一化二维 ISAR 幅度像与物理距离/方位轴，
%   直接调用 matlab_research_lib 中的 plotting.draw_isar_img 完成归一化 dB 渲染，
%   并支持按 saveOpts 导出图片。
%
% 输入参数：
%   isarResOrPost - 1×1 struct，isarpath_run 输出 isarRes（含 .post/.scene）或 isarPost
%   plotCfg       - (可选) 1×1 struct，绘图配置项：
%     .targetIndices   1×K double，待绘制的目标索引（默认全部目标 1:N）
%     .dynamicRangeDb  1×1 double，dB 显示动态范围（默认 30，对应 CLimDb=[-30, 0]）
%     .colormap        char/string，色图名（默认 'turbo'）
%     .showColorbar    1×1 logical，是否显示色标（默认 true）
%     .visible         char/string，图窗可见性 'on'（默认）或 'off'
%   saveOpts      - (可选) 1×1 struct，保存配置（enable/dir/fileName/format/dpi）
%
% 输出参数：
%   fig       - 1×1 Figure 句柄
%   savedPath - 1×1 string，导出文件路径（未保存时为 ""）

if nargin < 2 || isempty(plotCfg), plotCfg = struct(); end
if nargin < 3, saveOpts = struct(); end

if isfield(isarResOrPost, 'post')
    post = isarResOrPost.post;
    priorRanges = [];
    if isfield(isarResOrPost, 'scene') && isfield(isarResOrPost.scene, 'priorRangeM')
        priorRanges = isarResOrPost.scene.priorRangeM;
    end
else
    post = isarResOrPost;
    priorRanges = [];
end

nTotal = numel(post.image);
if isfield(plotCfg, 'targetIndices') && ~isempty(plotCfg.targetIndices)
    idx = double(plotCfg.targetIndices(:).');
else
    idx = 1:nTotal;
end
dynRangeDb   = double(get_cfg(plotCfg, 'dynamicRangeDb', 30));
cmapName     = char(string(get_cfg(plotCfg, 'colormap', 'turbo')));
showColorbar = logical(get_cfg(plotCfg, 'showColorbar', true));
figVisible   = char(string(get_cfg(plotCfg, 'visible', 'on')));

imgs       = post.image(idx);
rangeAxes  = post.rangeAxisM(idx);
crossAxes  = post.crossAxisM(idx);
titles     = strings(1, numel(idx));
for k = 1:numel(idx)
    tIdx = idx(k);
    if ~isempty(priorRanges) && numel(priorRanges) >= tIdx
        titles(k) = sprintf('Target %d, R_{prior}=%.1f m, valid=%d', ...
            tIdx, priorRanges(tIdx), post.valid(tIdx));
    else
        titles(k) = sprintf('Target %d, valid=%d', tIdx, post.valid(tIdx));
    end
end

% radar.ISAR_imaging 已完成第 1 维 fftshift，此处 ApplyFftshiftDim1 设为 false
fig = plotting.draw_isar_img(imgs, rangeAxes, crossAxes, ...
    'Visible',           figVisible, ...
    'ApplyFftshiftDim1', false, ...
    'NormalizeDb',       true, ...
    'CLimDb',            [-dynRangeDb, 0], ...
    'Colormap',          cmapName, ...
    'Titles',            titles, ...
    'XLabel',            'Range offset (m)', ...
    'YLabel',            'Cross-range (m)', ...
    'ShowColorbar',      showColorbar, ...
    'ColorbarLabel',     'Normalized magnitude (dB)');

savedPath = save_output_artifact(fig, "figure", saveOpts, 'isar_images');
end

function val = get_cfg(s, name, defaultVal)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)), val = s.(name); else, val = defaultVal; end
end
