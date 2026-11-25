function Tn = identify_isar_tip_noise_study(T, fits, outDir, cfg, trials)
%IDENTIFY_ISAR_TIP_NOISE_STUDY Local parameter uncertainty vs extraction noise.
%   The released tip tracks are nearly deterministic (residuals of
%   0.002-0.13 Hz), so their parameter errors are far below what a real
%   extraction would produce. This adds white noise to the extracted peak
%   track and refits the decided structure from the clean solution, which
%   measures local sensitivity rather than global basin selection.
%
%   T       identification table from identify_isar_tip_micro_motion
%   fits    per-scene fit structures (same row order as T)
%   outDir  output directory
%   cfg     isar_tip_motion_cfg
%   trials  noise realisations per level

if nargin < 5 || isempty(trials)
    trials = 8;
end
scenes = ["R2_s2_p3p5_a7_w0", "R3_s2_p3p5_a5_w3", "R3_s2_p3p5_a7_w5"];
sigmas = [0.5, 1, 2, 5];
rows = cell(0, 1);
rng(20261021);
for s = 1:numel(scenes)
    idx = find(string(T.caseId) == scenes(s), 1);
    if isempty(idx)
        continue;
    end
    f = fits{idx};
    modelName = f.best.model;
    theta0 = f.best.theta;
    for sg = sigmas
        aEst = nan(trials, 1);
        ampEst = nan(trials, 1);
        fpEst = nan(trials, 1);
        fwEst = nan(trials, 1);
        for k = 1:trials
            yN = f.y + sg * randn(numel(f.y), 1);
            o = isar_fit_tip_motion_mode(f.t, yN, modelName, ...
                struct('cfg', cfg, 'seed', theta0));
            aEst(k) = o.theta(2);
            ampEst(k) = o.theta(3);
            fpEst(k) = o.theta(1);
            fwEst(k) = o.theta(4);
        end
        rows{end + 1, 1} = struct('caseId', scenes(s), 'model', string(modelName), ...
            'noiseRmsHz', sg, 'trials', trials, ...
            'alphaBiasDeg', mean(aEst) - T.truthAlphaDeg(idx), ...
            'alphaStdDeg', std(aEst), ...
            'ampBiasDeg', mean(ampEst) - T.truthSwingAmpDeg(idx), ...
            'ampStdDeg', std(ampEst), 'fpStdHz', std(fpEst), ...
            'fwStdHz', std(fwEst)); %#ok<AGROW>
        fprintf('noise %4.1f Hz | %-22s alpha bias %+8.4f std %7.4f | A bias %+7.4f std %7.4f\n', ...
            sg, scenes(s), mean(aEst) - T.truthAlphaDeg(idx), std(aEst), ...
            mean(ampEst) - T.truthSwingAmpDeg(idx), std(ampEst));
    end
end
if isempty(rows)
    Tn = table();
    return;
end
Tn = struct2table([rows{:}]);
writetable(Tn, fullfile(outDir, 'noise_robustness.csv'));

font = 'Microsoft YaHei';
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [40 40 2000 900]);
layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
names = unique(Tn.caseId, 'stable');
cols = lines(numel(names));
for panel = 1:2
    ax = nexttile(layout, panel);
    hold(ax, 'on');
    for k = 1:numel(names)
        sub = Tn(Tn.caseId == names(k), :);
        if panel == 1
            errorbar(ax, sub.noiseRmsHz, sub.alphaStdDeg, sub.alphaStdDeg, '-o', ...
                'Color', cols(k, :), 'LineWidth', 2, 'MarkerFaceColor', cols(k, :));
        else
            errorbar(ax, sub.noiseRmsHz, sub.ampStdDeg, sub.ampStdDeg, '-o', ...
                'Color', cols(k, :), 'LineWidth', 2, 'MarkerFaceColor', cols(k, :));
        end
    end
    set(ax, 'XScale', 'log', 'YScale', 'log');
    grid(ax, 'on');
    ax.FontSize = 22;
    ax.FontName = font;
    xlabel(ax, '叠加到提取轨迹的白噪声标准差 (Hz)', 'FontSize', 22, 'FontName', font);
    if panel == 1
        ylabel(ax, '半锥角估计标准差 (°)', 'FontSize', 22, 'FontName', font);
        title(ax, '进动半锥角 α', 'FontSize', 26, 'FontName', font, 'FontWeight', 'normal');
    else
        ylabel(ax, '摆动幅度估计标准差 (°)', 'FontSize', 22, 'FontName', font);
        title(ax, '摆动幅度 A', 'FontSize', 26, 'FontName', font, 'FontWeight', 'normal');
    end
end
legend(nexttile(layout, 1), cellstr(names), 'Location', 'northwest', ...
    'FontSize', 20, 'FontName', font);
sgtitle(layout, '参数估计精度随提取噪声的变化（局部重拟合）', 'FontSize', 28, ...
    'FontName', font, 'FontWeight', 'normal');
exportgraphics(fig, fullfile(outDir, 'tip_mode_noise_robustness.png'), 'Resolution', 120);
savefig(fig, fullfile(outDir, 'tip_mode_noise_robustness.fig'));
close(fig);
end
