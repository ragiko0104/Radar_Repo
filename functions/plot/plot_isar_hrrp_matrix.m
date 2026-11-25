function [fig, savedPath] = plot_isar_hrrp_matrix(isarResOrPost, plotCfg, saveOpts)
%PLOT_ISAR_HRRP_MATRIX 绘制单目标各处理阶段的 HRRP-慢时间二维矩阵及 RCMC 位移曲线。
%
% 功能说明：
%   将单目标在不同处理阶段的复一维距离像-慢时间矩阵（Na×Nr）按归一化 dB 显示，
%   便于直观对比多普勒色散补偿、RCMC 包络对齐和特显点相位聚焦前后的变化；
%   可选附加 RCMC 逐脉冲移位量曲线，并支持按 saveOpts 导出图片。
%
% 输入参数：
%   isarResOrPost - 1×1 struct，isarpath_run 输出 isarRes（含 .post/.scene）或 isarPost
%   plotCfg       - (可选) 1×1 struct，绘图配置项：
%     .targetIndex     1×1 double，目标索引（默认 1）
%     .stages          1×M string，选画阶段列表，可选：
%                        "beforeDispersion" | "afterDispersion" | "afterRcmc" | "afterPhaseComp"
%                      默认 ["beforeDispersion", "afterDispersion", "afterRcmc", "afterPhaseComp"]
%     .showRcmcCurve   1×1 logical，是否附加 RCMC 逐脉冲距离移位曲线（默认 true）
%     .dynamicRangeDb  1×1 double，HRRP dB 显示动态范围（默认 40）
%     .colormap        char/string，色图名（默认 'turbo'）
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
    if isfield(isarResOrPost, 'scene') && isfield(isarResOrPost.scene, 'pulseTimeSec')
        slowTimeMs = isarResOrPost.scene.pulseTimeSec(:) * 1e3;
    else
        slowTimeMs = [];
    end
else
    post = isarResOrPost;
    slowTimeMs = [];
end

targetIdx    = double(get_cfg(plotCfg, 'targetIndex', 1));
stages       = string(get_cfg(plotCfg, 'stages', ...
    ["beforeDispersion", "afterDispersion", "afterRcmc", "afterPhaseComp"]));
showRcmc     = logical(get_cfg(plotCfg, 'showRcmcCurve', true));
dynRangeDb   = double(get_cfg(plotCfg, 'dynamicRangeDb', 40));
cmapName     = char(string(get_cfg(plotCfg, 'colormap', 'turbo')));
figVisible   = char(string(get_cfg(plotCfg, 'visible', 'on')));

rangeAxisM = post.rangeAxisM{targetIdx};
nPulses    = size(post.hrrfModelCrop{targetIdx}, 1);
if isempty(slowTimeMs) || numel(slowTimeMs) ~= nPulses
    slowTimeAxis = 1:nPulses;
    yLabelText   = 'Pulse index';
else
    slowTimeAxis = slowTimeMs;
    yLabelText   = 'Slow time (ms)';
end

stageMats   = {};
stageTitles = {};
for k = 1:numel(stages)
    switch lower(stages(k))
        case "beforedispersion"
            if isfield(post, 'hrrfBeforeDispersion') && ~isempty(post.hrrfBeforeDispersion{targetIdx})
                stageMats{end+1}   = post.hrrfBeforeDispersion{targetIdx}; %#ok<AGROW>
                stageTitles{end+1} = 'HRRP (Before Dispersion Comp)'; %#ok<AGROW>
            end
        case "afterdispersion"
            stageMats{end+1}   = post.hrrfModelCrop{targetIdx}; %#ok<AGROW>
            stageTitles{end+1} = 'HRRP (After Dispersion Comp)'; %#ok<AGROW>
        case "afterrcmc"
            stageMats{end+1}   = post.hrrfAligned{targetIdx}; %#ok<AGROW>
            stageTitles{end+1} = 'HRRP (After RCMC)'; %#ok<AGROW>
        case "afterphasecomp"
            stageMats{end+1}   = post.phaseCompensatedHrrf{targetIdx}; %#ok<AGROW>
            stageTitles{end+1} = 'HRRP (After Phase Comp)'; %#ok<AGROW>
        otherwise
            error('RV_ISAR_Path:UnsupportedHrrpStage', '不支持的 HRRP 阶段 "%s"。', stages(k));
    end
end

nPanels = numel(stageMats) + double(showRcmc && isfield(post, 'rangeMigrationBins'));
nCols   = min(2, nPanels);
nRows   = ceil(nPanels / nCols);

fig = figure('Name', sprintf('ISAR HRRP-SlowTime (Target %d)', targetIdx), ...
    'Color', 'w', 'Visible', figVisible);
layout = tiledlayout(fig, nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

for k = 1:numel(stageMats)
    ax = nexttile(layout);
    mag = abs(stageMats{k});
    peakVal = max(mag(:));
    if peakVal > 0
        imgDb = 20 * log10(max(mag / peakVal, 10^(-dynRangeDb / 20)));
    else
        imgDb = -dynRangeDb * ones(size(mag));
    end
    imagesc(ax, rangeAxisM, slowTimeAxis, imgDb);
    axis(ax, 'xy');
    clim(ax, [-dynRangeDb, 0]);
    colormap(ax, cmapName);
    cb = colorbar(ax);
    cb.Label.String = 'dB';
    xlabel(ax, 'Range offset (m)');
    ylabel(ax, yLabelText);
    title(ax, stageTitles{k});
end

if showRcmc && isfield(post, 'rangeMigrationBins')
    ax = nexttile(layout);
    plot(ax, slowTimeAxis, post.rangeMigrationBins(targetIdx, :), '-o', 'LineWidth', 1.3, 'MarkerSize', 4);
    grid(ax, 'on');
    xlabel(ax, yLabelText);
    ylabel(ax, 'RCMC shift (bins)');
    title(ax, sprintf('RCMC Shift (mode: %s)', string(post.parameters.rangeMigrationMode)));
end

title(layout, sprintf('Target %d — HRRP-SlowTime Matrices', targetIdx));

savedPath = save_output_artifact(fig, "figure", saveOpts, ...
    sprintf('isar_hrrp_matrix_target%02d', targetIdx));
end

function val = get_cfg(s, name, defaultVal)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)), val = s.(name); else, val = defaultVal; end
end
