function [fig, savedPath] = plot_rv_time_frequency(rvRes, plotCfg, saveOpts)
%PLOT_RV_TIME_FREQUENCY 绘制 RV 链路去载频六段 LFM 回波时频图。
%
% 功能说明：
%   复用 lib 中的 plotting.spectrogram_uniformed，绘制指定慢时间脉冲的去载频
%   六段复合基带回波时频分布，并支持按 saveOpts 导出图片。
%
% 输入参数：
%   rvRes    - 1×1 struct，rvpath_run 产出的结果结构体（使用 echoBaseband Na×N、cfg.fs、cfg.k_v）
%   plotCfg  - (可选) 1×1 struct，绘图配置项：
%     .pulseIndex  1×1 double，待绘制的脉冲索引（1..Na，默认 1）
%     .visible     char/string，图窗可见性 'on'（默认）或 'off'
%   saveOpts - (可选) 1×1 struct，保存配置（enable/dir/fileName/format/dpi）
%
% 输出参数：
%   fig       - 1×1 Figure 句柄
%   savedPath - 1×1 string，导出文件路径（未保存时为 ""）

if nargin < 2 || isempty(plotCfg), plotCfg = struct(); end
if nargin < 3, saveOpts = struct(); end

wf = rvRes.cfg;
Na = wf.Na;
fs = wf.fs;

pulseIdx   = double(get_cfg(plotCfg, 'pulseIndex', 1));
figVisible = char(string(get_cfg(plotCfg, 'visible', 'on')));
assert(isscalar(pulseIdx) && pulseIdx >= 1 && pulseIdx <= Na && pulseIdx == floor(pulseIdx), ...
    'plotCfg.pulseIndex 必须是 1 到 %d 之间的整数。', Na);

nfft    = 2^max(1, nextpow2(fs / sqrt(max(abs(wf.k_v)))) - 1);
overlap = floor(nfft / 2);

fig = figure('Name', 'RV simulated echo time-frequency', 'Color', 'w', 'Visible', figVisible);
plotting.spectrogram_uniformed( ...
    rvRes.echoBaseband(pulseIdx, :), nfft, overlap, nfft, fs, 'yaxis');
title(sprintf('第 %d 个仿真测量脉冲时频图（已去除载频）', pulseIdx));

savedPath = save_output_artifact(fig, "figure", saveOpts, ...
    sprintf('rv_time_frequency_p%d', pulseIdx));
end

function val = get_cfg(s, name, defaultVal)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)), val = s.(name); else, val = defaultVal; end
end
