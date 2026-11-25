function result = sweep_isar_micro_motion_limits(options)
%SWEEP_ISAR_MICRO_MOTION_LIMITS Sensitivity window of the frozen pipeline.
%   Tip modes use the four-model empirical F ladder. Spin detection uses
%   the F, frequency and joint-RMS gates of identify_isar_spin_from_tail_points.
%
%   For each parameter, one value at a time is swept around a baseline, with
%   independent noise realisations, and every trial is classified as a
%   success or failure by the frozen decision rules:
%
%     fp     model not static and |fp - truth| < 0.10 Hz  (half a slow-time bin)
%     alpha    |alpha - truth| < max(0.2*alpha, 0.05 deg)
%     A        swing declared AND |A - truth| < max(0.2*A, 0.2 deg)
%     fw       swing declared AND |fw - truth| < max(0.2*fw, 0.1 Hz)
%     fs       spin declared AND |fs - truth| < max(0.2*fs, 0.1 Hz)
%
%   Sweeps: fp, alpha, A, fw and fs; aperture length (the 1/T slow-time
%   resolution limit); and two coupling probes, fw near fp and fs near fp.
%
%   OPTIONS: Smoke, Trials, NoiseHz, NoiseHzForSpin, OutputDirectory,
%   GenerateFigures. SpinOnly reruns the two spin grids and writes their CSV.
%
%   Noise: the tip track carries 0.2 Hz white noise (measured residual of the
%   real tip tracks is 0.003-0.5 Hz). The spin stage carries 1 Hz white noise
%   per point by default, which stands in for the 7-8 Hz systematic residual
%   seen on the real base tracks; NoiseHzForSpin lets that be varied because
%   the spin limit scales with it.

if nargin < 1 || isempty(options)
    options = struct();
end
root = fileparts(fileparts(mfilename('fullpath')));
outDir = char(opt(options, 'OutputDirectory', ...
    fullfile(root, 'results', 'isar_sensitivity_limits_20261021')));
if ~isfolder(outDir)
    mkdir(outDir);
end
smoke = opt(options, 'Smoke', false);
nTrial = opt(options, 'Trials', 5);
sigmaTip = opt(options, 'NoiseHz', 0.20);
sigmaSpin = opt(options, 'NoiseHzForSpin', 1.0);
cfg = isar_tip_motion_cfg();
P = isar_cone_body_points();
P1 = P(1, :);
tFull = 0.0156 + (0:49)' * cfg.sampleStepSec;
rng(20261021);

base = struct('fs', 1, 'fp', 1.7, 'alpha', 3, 'A', 0, 'fw', 0.7);
if smoke
    grids = struct('fp', [0.2 1], 'alpha', [0.1 3], 'A', [0.05 3], ...
        'fw', [0.2 2], 'fs', [0.3 1], 'aperture', [5 1], ...
        'fwNearFp', [1.7 3.4], 'fsNearFp', [1.7 3.4]);
    nTrial = 1;
else
    grids = struct( ...
        'fp', [0.05 0.10 0.20 0.30 0.50 1 2 3.5 5 7], ...
        'alpha', [0.03 0.05 0.10 0.20 0.30 0.50 1 3 7 12], ...
        'A', [0.02 0.05 0.10 0.20 0.30 0.50 1 3 6 10], ...
        'fw', [0.05 0.10 0.20 0.50 0.70 1 2 3 5], ...
        'fs', [0.05 0.10 0.20 0.30 0.50 1 2 3 4 5.5], ...
        'aperture', [5 3 2 1], ...
        'fwNearFp', [1.0 1.4 1.6 1.7 1.8 2.0 2.4], ...
        'fsNearFp', [1.0 1.4 1.6 1.7 1.8 2.2 3.4]);
end

timer = tic;
if opt(options, 'SpinOnly', false)
    result = struct('spinTable', spinSweeps(base, grids, nTrial, sigmaSpin, outDir, timer), ...
        'configuration', cfg, 'trials', nTrial, 'sigmaSpin', sigmaSpin, 'baseline', base);
    return
end
rows = cell(0, 1);

% ---- tip-only sweeps ------------------------------------------------------
specs = {'fp', 'fp'; 'alpha', 'alpha'; 'A', 'A'; 'fw', 'fw'};
for s = 1:size(specs, 1)
    name = specs{s, 1};
    for v = grids.(name)
        q = base;
        q.(specs{s, 2}) = v;
        if strcmp(name, 'fw')
            q.A = 3;                     % swing must exist for fw to mean anything
        end
        for k = 1:nTrial
            r = evalTip(q, tFull, P1, cfg, sigmaTip);
            rows{end + 1, 1} = packRow(name, v, k, r, q, nTrial); %#ok<AGROW>
        end
        fprintf('%-6s %7.3g | P(success) %5.2f | %.0f s\n', name, v, ...
            successRate(rows, name, v), toc(timer));
    end
end
% ---- aperture sweep (slow-time resolution 1/T) ---------------------------
for T = grids.aperture
    nCpi = max(round(T / cfg.sampleStepSec), 5);
    tSeg = tFull(1:nCpi);
    q = base;
    q.alpha = 1;
    for k = 1:nTrial
        r = evalTip(q, tSeg, P1, cfg, sigmaTip);
        r.obsSec = T;
        rows{end + 1, 1} = packRow('aperture', T, k, r, q, nTrial); %#ok<AGROW>
    end
    fprintf('aperture %4.1f s | P(success) %5.2f | %.0f s\n', T, ...
        successRate(rows, 'aperture', T), toc(timer));
end
% ---- coupling probe: swing frequency near the precession frequency --------
for v = grids.fwNearFp
    q = base;
    q.A = 3;
    q.fw = v;
    for k = 1:nTrial
        r = evalTip(q, tFull, P1, cfg, sigmaTip);
        rows{end + 1, 1} = packRow('fwNearFp', v, k, r, q, nTrial); %#ok<AGROW>
    end
    fprintf('fw~fp %5.2f | P(success) %5.2f | %.0f s\n', v, ...
        successRate(rows, 'fwNearFp', v), toc(timer));
end
table1 = struct2table([rows{:}]);
writetable(table1, fullfile(outDir, 'sensitivity_sweeps.csv'));

table2 = spinSweeps(base, grids, nTrial, sigmaSpin, outDir, timer);
summary = struct('configuration', cfg, 'trials', nTrial, 'smoke', smoke, ...
    'sigmaTip', sigmaTip, 'sigmaSpin', sigmaSpin, 'baseline', base, ...
    'grids', grids, 'tipTable', table1, 'spinTable', table2);
save(fullfile(outDir, 'sensitivity_limits.mat'), 'summary');
if opt(options, 'GenerateFigures', true)
    render_figures(table1, table2, outDir, grids);
end
result = summary;
fprintf('DONE %s | %.0f s\n', outDir, toc(timer));
end

function table2 = spinSweeps(base, grids, nTrial, sigmaSpin, outDir, timer)
spinRows = cell(0, 1);
for v = grids.fs
    for k = 1:nTrial
        spinRows{end + 1, 1} = evalSpin(base, 'fs', v, k, ...
            20261021 + numel(spinRows), sigmaSpin, outDir); %#ok<AGROW>
    end
    fprintf('fs     %5.2f | P(success) %5.2f | %.0f s\n', v, ...
        spinSuccess(spinRows, 'fs', v), toc(timer));
end
for v = grids.fsNearFp
    for k = 1:nTrial
        spinRows{end + 1, 1} = evalSpin(base, 'fsNearFp', v, k, ...
            20261021 + numel(spinRows), sigmaSpin, outDir); %#ok<AGROW>
    end
    fprintf('fs~fp  %5.2f | P(success) %5.2f | %.0f s\n', v, ...
        spinSuccess(spinRows, 'fsNearFp', v), toc(timer));
end
table2 = struct2table([spinRows{:}]);
writetable(table2, fullfile(outDir, 'spin_sensitivity.csv'));
end

% ======================================================================
function r = evalTip(q, tTrack, P, cfg, sigma)
p = [q.fs, q.fp, q.alpha, q.A, q.fw, 0, 0, 0];
[Y, ~] = isar_simulate_cone_peaks(tTrack, P, p, cfg, sigma, [5.15, 0.14]);
% Frozen pipeline now includes the H1/H3 parsimony decision: a swing must
% beat the pure-precession fit by an F test on its three extra parameters,
% otherwise A is reported as zero. See isar_fit_tip_parsimonious.
fit = isar_fit_tip_parsimonious(tTrack, Y(:, 1), cfg, struct('ranking', 'short'));
th = fit.theta;
% The swing flag must come from the algorithm's own decision (the F ladder),
% not from an amplitude threshold applied here: scoring with the old 0.2 deg
% rule while the algorithm decides by F mixed two criteria and marked correct
% low-amplitude detections as failures.
r = struct('fp', th(1), 'alpha', th(2), 'A', th(3), 'fw', th(4), ...
    'tipRms', fit.rmsHz, 'swing', fit.parsimony.keptSwing, ...
    'obsSec', tTrack(end) - tTrack(1) + cfg.sampleStepSec);
end

function row = packRow(name, value, trial, r, q, nTrial) %#ok<INUSD>
row = struct('sweep', string(name), 'value', value, 'trial', trial, ...
    'trueFp', q.fp, 'trueAlpha', q.alpha, 'trueA', q.A, 'trueFw', q.fw, ...
    'estFp', r.fp, 'estAlpha', r.alpha, 'estA', r.A, 'estFw', r.fw, ...
    'swingDetected', r.swing, 'tipRms', r.tipRms, 'obsSec', r.obsSec);
switch name
    case 'fp'
        row.ok = abs(r.fp - q.fp) < 0.10 && r.fp > 0.05;
    case 'alpha'
        row.ok = abs(r.alpha - q.alpha) < max(0.2 * q.alpha, 0.05) && r.alpha > 0.01;
    case 'A'
        row.ok = r.swing && abs(r.A - q.A) < max(0.2 * q.A, 0.2);
    case 'fw'
        row.ok = r.swing && abs(r.fw - q.fw) < max(0.2 * q.fw, 0.1);
    case 'aperture'
        row.ok = abs(r.fp - q.fp) < 0.15 && abs(r.alpha - q.alpha) < max(0.2 * q.alpha, 0.1);
    otherwise   % coupling probes: the swing must be recovered, not merged into fp
        row.ok = r.swing && abs(r.fw - q.fw) < max(0.2 * q.fw, 0.15) ...
            && abs(r.A - q.A) < max(0.2 * q.A, 0.2);
end
end

function p = successRate(rows, name, value)
sel = cellfun(@(z) strcmp(char(z.sweep), name) && z.value == value, rows);
p = mean(cellfun(@(z) z.ok, rows(sel)));
end

function row = evalSpin(base, name, value, trial, noiseSeed, sigma, outDir)
q = base;
q.fs = value;
one = table(string(sprintf('%s_%g', name, value)), q.fs, q.fp, q.alpha, q.A, q.fw, ...
    'VariableNames', {'caseId', 'spinFrequencyHz', 'precessionFrequencyHz', ...
    'precessionHalfAngleDeg', 'swingAmplitudeDeg', 'swingFrequencyHz'});
tmp = fullfile(outDir, '_spin_tmp');
out = identify_isar_spin_from_tail_points(struct('Cases', one, ...
    'NoiseHz', sigma, 'NoiseSeed', noiseSeed, 'GenerateFigures', false, 'OutputDirectory', tmp));
T = out.table;
row = struct('sweep', string(name), 'value', value, 'trial', trial, ...
    'noiseSeed', noiseSeed, ...
    'trueFs', q.fs, 'trueFp', q.fp, 'trueAlpha', q.alpha, ...
    'estFs', T.estSpinHz(1), 'estFp', T.estFpHz(1), 'estAlpha', T.estAlphaDeg(1), ...
    'Fstat', T.Fstat(1), 'spinDetected', T.spinDetected(1), ...
    'jointRms', T.jointRmsHz(1), 'foldCount', T.foldCount(1));
row.ok = T.spinDetected(1) && abs(T.estSpinHz(1) - q.fs) < max(0.2 * q.fs, 0.1) ...
    && abs(T.estFpHz(1) - q.fp) < 0.1;
end

function p = spinSuccess(rows, name, value)
sel = cellfun(@(z) strcmp(char(z.sweep), name) && z.value == value, rows);
p = mean(cellfun(@(z) z.ok, rows(sel)));
end

function render_figures(T, S, outDir, grids)
font = 'Microsoft YaHei';
names = {'fp', 'alpha', 'A', 'fw'};
titles = {'进动频率 f_p (Hz)', '进动半锥角 α (°)', '摆动幅度 A (°)', '摆动频率 f_w (Hz)'};
xlab = {'f_p (Hz)', 'α (°)', 'A (°)', 'f_w (Hz)'};
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [40 40 2400 1300]);
layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:numel(names)
    ax = nexttile(layout, k);
    vals = grids.(names{k});
    p = arrayfun(@(v) mean(T.ok(T.sweep == names{k} & T.value == v)), vals);
    plot(ax, vals, p, '-o', 'Color', [0.00 0.45 0.35], 'LineWidth', 2.4, ...
        'MarkerFaceColor', [0.00 0.45 0.35]);
    hold(ax, 'on');
    yline(ax, 0.9, ':', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.6);
    set(ax, 'XScale', 'log');
    grid(ax, 'on');
    xlim(ax, [min(vals) * 0.9, max(vals) * 1.1]);
    ylim(ax, [-0.05 1.05]);
    ax.FontSize = 20; ax.FontName = font;
    title(ax, titles{k}, 'FontSize', 24, 'FontName', font, 'FontWeight', 'normal');
    xlabel(ax, xlab{k}, 'FontSize', 20, 'FontName', font);
    ylabel(ax, '成功率', 'FontSize', 20, 'FontName', font);
end
sgtitle(layout, '冻结算法的参数识别成功率（红线 90%）', 'FontSize', 28, ...
    'FontName', font, 'FontWeight', 'normal');
exportgraphics(fig, fullfile(outDir, 'sensitivity_tip.png'), 'Resolution', 120);
close(fig);

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [40 40 2200 900]);
layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(layout, 1);
vals = grids.fs;
p = arrayfun(@(v) mean(S.ok(S.sweep == "fs" & S.value == v)), vals);
plot(ax, vals, p, '-o', 'Color', [0.07 0.32 0.60], 'LineWidth', 2.4, ...
    'MarkerFaceColor', [0.07 0.32 0.60]);
hold(ax, 'on'); yline(ax, 0.9, ':', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.6);
grid(ax, 'on'); ax.FontSize = 20; ax.FontName = font;
xlabel(ax, 'f_s (Hz)', 'FontSize', 20, 'FontName', font);
ylabel(ax, '成功率', 'FontSize', 20, 'FontName', font);
title(ax, '自旋频率 f_s', 'FontSize', 24, 'FontName', font, 'FontWeight', 'normal');
ax = nexttile(layout, 2);
vals2 = grids.fsNearFp;
p2 = arrayfun(@(v) mean(S.ok(S.sweep == "fsNearFp" & S.value == v)), vals2);
plot(ax, vals2, p2, '-o', 'Color', [0.60 0.20 0.60], 'LineWidth', 2.4, ...
    'MarkerFaceColor', [0.60 0.20 0.60]);
hold(ax, 'on'); yline(ax, 0.9, ':', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.6);
grid(ax, 'on'); ax.FontSize = 20; ax.FontName = font;
xlabel(ax, 'f_s (Hz)，f_p = 1.7 Hz', 'FontSize', 20, 'FontName', font);
ylabel(ax, '成功率', 'FontSize', 20, 'FontName', font);
title(ax, '自旋与进动接近时的耦合', 'FontSize', 24, 'FontName', font, 'FontWeight', 'normal');
sgtitle(layout, '自旋识别窗口与 f_s–f_p 耦合', 'FontSize', 28, ...
    'FontName', font, 'FontWeight', 'normal');
exportgraphics(fig, fullfile(outDir, 'sensitivity_spin.png'), 'Resolution', 120);
close(fig);
end

function v = opt(s, name, fallback)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = fallback;
end
end
