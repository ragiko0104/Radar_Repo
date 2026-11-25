function [fig, savedPath] = plot_rv_pulse_compression(rvRes, plotCfg, saveOpts)
%PLOT_RV_PULSE_COMPRESSION 绘制 RV 链路 3×2 六段脉冲压缩功率、CFAR 与 top-k 阈值图。
%
% 功能说明：
%   在单个 3×2 布局图窗中展示六段 LFM（频带 1–3 上下扫）的脉压线性功率、
%   MEAN-CFAR 门限与 top-k 峰均值二次门限，支持按 saveOpts 导出图片。
%
% 输入参数：
%   rvRes    - 1×1 struct，rvpath_run 产出的结果结构体（使用 pulseCompression/cfar/cfg）
%   plotCfg  - (可选) 1×1 struct，绘图配置项：
%     .pulseIndex  1×1 double，待绘制的脉冲索引（1..Na，默认 1）
%     .c           1×1 double，光速 (m/s，默认 3e8)
%     .visible     char/string，图窗可见性 'on'（默认）或 'off'
%   saveOpts - (可选) 1×1 struct，保存配置（enable/dir/fileName/format/dpi）
%
% 输出参数：
%   fig       - 1×1 Figure 句柄
%   savedPath - 1×1 string，导出文件路径（未保存时为 ""）

if nargin < 2 || isempty(plotCfg), plotCfg = struct(); end
if nargin < 3, saveOpts = struct(); end

wf     = rvRes.cfg;
pc     = rvRes.pulseCompression;
rvCfar = rvRes.cfar;
Na     = wf.Na;
N      = wf.N;

pulseIdx   = double(get_cfg(plotCfg, 'pulseIndex', 1));
c          = double(get_cfg(plotCfg, 'c', 3e8));
figVisible = char(string(get_cfg(plotCfg, 'visible', 'on')));
assert(isscalar(pulseIdx) && pulseIdx >= 1 && pulseIdx <= Na && pulseIdx == floor(pulseIdx), ...
    'plotCfg.pulseIndex 必须是 1 到 %d 之间的整数。', Na);

rangeAxisM = wf.fastTimeSec .* c ./ 2;
fig = figure('Name', 'RV pulse compression and CFAR thresholds', 'Color', 'w', 'Visible', figVisible);
layout = tiledlayout(fig, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for segIdx = 1:numel(rvCfar.segmentKeys)
    ax = nexttile(layout, segIdx);
    pwr = abs(pc.compressed{segIdx}(pulseIdx, :)).^2;
    plot(ax, rangeAxisM, pwr, 'LineWidth', 1);
    hold(ax, 'on');
    plot(ax, rangeAxisM, rvCfar.threshold{segIdx}(pulseIdx, :), 'LineWidth', 1.2);
    plot(ax, rangeAxisM, rvCfar.topKMeanThreshold(pulseIdx, segIdx) .* ones(1, N), 'LineWidth', 1.2);
    grid(ax, 'on');
    xlabel(ax, '距离 (m)');
    ylabel(ax, '功率');
    title(ax, rvCfar.segmentNames(segIdx));
    if segIdx == 1
        legend(ax, {'脉冲压缩功率', 'CFAR 阈值', 'top-k 峰均值阈值'}, 'Location', 'best');
    end
end
title(layout, sprintf('第 %d 个脉冲：六段 LFM 脉压与 CFAR 结果（本振：%s）', ...
    pulseIdx, pc.localOscillatorMode));

savedPath = save_output_artifact(fig, "figure", saveOpts, ...
    sprintf('rv_pulse_compression_p%d', pulseIdx));
end

function val = get_cfg(s, name, defaultVal)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)), val = s.(name); else, val = defaultVal; end
end
