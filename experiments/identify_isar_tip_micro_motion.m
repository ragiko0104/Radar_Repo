function result = identify_isar_tip_micro_motion(options)
%IDENTIFY_ISAR_TIP_MICRO_MOTION Blind tip micro-motion mode identification.
%   RESULT = IDENTIFY_ISAR_TIP_MICRO_MOTION(OPTIONS) consumes the extracted
%   cone-tip Doppler tracks of run_isar_tip_motion_retest (the 14-scene
%   release) and identifies the micro-motion structure without using the
%   scenario truth: it fits the candidate hypotheses of
%   isar_tip_motion_bank, reports the estimated parameters with their
%   uncertainties, and only then compares them with the theoretical prior.
%
%   The observable is the per-CPI peak of the tip range gate, i.e. 50 samples
%   of the dwell-averaged tip Doppler at 10 Hz over 5 s. Nothing in the fit
%   uses the scene label, the true frequencies or the true angles.
%
%   OPTIONS
%     ResultDirectory    directory holding all_scene_similarity.csv and the
%                        groupN/tip_extracted_tracks.csv files
%     OutputDirectory    where the tables, figures and summary are written
%     Scenes             subset of caseId values (default: all)
%     GenerateFigures    default true
%     SwingPromotionF    F3 cutoff for adding swing to precession
%                        (default: cfg.parsimonyF)
%
%   RESULT holds the identification table, the per-scene fit structures and
%   the configuration actually used.
%
%   See also ISAR_FIT_TIP_MOTION_MODE, ISAR_TIP_MOTION_BANK,
%   ISAR_TIP_MOTION_KINEMATICS, test_isar_tip_motion_mode_identification.

if nargin < 1 || isempty(options)
    options = struct();
end
root = fileparts(fileparts(mfilename('fullpath')));
resultDir = char(opt(options, 'ResultDirectory', ...
    fullfile(root, 'results', 'isar_tip_motion_retest_20261020')));
outDir = char(opt(options, 'OutputDirectory', ...
    fullfile(root, 'results', 'isar_tip_motion_identification_20261021')));
if ~isfolder(outDir)
    mkdir(outDir);
end
cfg = isar_tip_motion_cfg();
cfg.swingPromotionF = opt(options, 'SwingPromotionF', cfg.parsimonyF);
swingThresholdDeg = cfg.swingDecisionDeg;

truthFile = fullfile(resultDir, 'all_scene_similarity.csv');
if ~isfile(truthFile)
    error('RadarProject:TipMode:MissingResults', ...
        'Missing %s. Run run_isar_tip_motion_retest first.', truthFile);
end
truth = readtable(truthFile, 'TextType', 'string');
if isfield(options, 'Scenes') && ~isempty(options.Scenes)
    wanted = string(options.Scenes);
    truth = truth(ismember(string(truth.caseId), wanted), :);
end
% Stable ordering: group, then caseId.
truth = sortrows(truth, {'family', 'caseId'});
nScene = height(truth);

rows = cell(nScene, 1);
fits = cell(nScene, 1);
bank = isar_tip_motion_bank();
for s = 1:nScene
    cid = string(truth.caseId(s));
    group = "group" + extractBetween(cid, "R", "_");
    group = group(1);
    trackFile = fullfile(resultDir, char(group), 'tip_extracted_tracks.csv');
    tracks = readtable(trackFile, 'TextType', 'string');
    tracks = tracks(string(tracks.caseId) == cid, :);
    if isempty(tracks)
        error('RadarProject:TipMode:MissingTrack', 'No track rows for %s.', cid);
    end
    tracks = sortrows(tracks, 'receiveMidpointSec');
    t = tracks.receiveMidpointSec;
    y = tracks.extractedTipPeakHz;

    fprintf('[%2d/%2d] %-22s truth fp %4.1f a %2.0f A %2.0f fw %3.1f | ', ...
        s, nScene, char(cid), truth.precessionFrequencyHz(s), ...
        truth.precessionHalfAngleDeg(s), truth.swingAmplitudeDeg(s), ...
        truth.swingFrequencyHz(s));
    % One decision path for every caller: the four hypotheses are compared by
    % the nested F ladder of isar_identify_tip_mode, which replaced the
    % earlier BIC ranking plus fixed 0.2 deg amplitude rule here and a
    % separate H1/H3 test in isar_fit_tip_parsimonious.
    id = isar_identify_tip_mode(t, y, cfg);
    models = id.fits;
    best = id.selected;
    th = best.theta;
    swingDetected = id.swing;
    decidedMode = id.mode;
    truthHasSwing = truth.swingAmplitudeDeg(s) > 0 && truth.swingFrequencyHz(s) > 0;
    truthHasPrecession = truth.precessionHalfAngleDeg(s) > 0 && truth.precessionFrequencyHz(s) > 0;
    truthIndex = 1 + truthHasPrecession + 2 * truthHasSwing;

    row = struct();
    row.caseId = cid;
    row.truthMode = string(truth.mode(s));
    row.truthTipMode = string(bank(truthIndex).description);
    row.truthFpHz = truth.precessionFrequencyHz(s);
    row.truthAlphaDeg = truth.precessionHalfAngleDeg(s);
    row.truthSwingAmpDeg = truth.swingAmplitudeDeg(s);
    row.truthSwingFreqHz = truth.swingFrequencyHz(s);
    row.selectedModel = string(best.model);
    row.decidedMode = decidedMode;
    bics = cellfun(@(z) z.bic, models);
    others = bics;
    others(strcmp({bank.name}, best.model)) = [];
    row.bicMarginToRunnerUp = min(others) - best.bic;
    row.estFpHz = th(1);
    row.seFpHz = best.sd(1);
    row.estAlphaDeg = th(2);
    row.seAlphaDeg = best.sd(2);
    row.estSwingAmpDeg = th(3);
    row.seSwingAmpDeg = best.sd(3);
    row.estSwingFreqHz = th(4);
    row.seSwingFreqHz = best.sd(4);
    row.swingComponentHz = best.swingComponentHz;
    row.swingSnr = best.swingComponentHz / max(best.rmsHz, eps);
    row.residualRmsHz = best.rmsHz;
    row.bicH0 = models{1}.bic;
    row.bicH1 = models{2}.bic;
    row.bicH2 = models{3}.bic;
    row.bicH3 = models{4}.bic;
    row.swingDetected = swingDetected;
    row.modeCorrect = decidedMode == row.truthTipMode;
    row.errFpHz = th(1) - truth.precessionFrequencyHz(s);
    row.errFpPct = 100 * row.errFpHz / max(truth.precessionFrequencyHz(s), eps);
    row.errAlphaDeg = th(2) - truth.precessionHalfAngleDeg(s);
    row.errAlphaPct = 100 * row.errAlphaDeg / max(truth.precessionHalfAngleDeg(s), eps);
    row.errSwingAmpDeg = th(3) - truth.swingAmplitudeDeg(s);
    if truth.swingFrequencyHz(s) > 0
        row.errSwingFreqHz = th(4) - truth.swingFrequencyHz(s);
    else
        row.errSwingFreqHz = NaN;
    end
    rows{s} = row;
    fits{s} = struct('caseId', cid, 't', t, 'y', y, 'models', {models}, ...
        'best', best, 'truth', truth(s, :));

    fprintf('est fp %6.3f a %7.4f A %7.4f fw %6.4f | %s %s | rms %8.4f Hz\n', ...
        th(1), th(2), th(3), th(4), char(decidedMode), ...
        ternary(swingDetected, 'SWING', 'NOSWG'), best.rmsHz);
end

tableOut = struct2table([rows{:}]);
writetable(tableOut, fullfile(outDir, 'identified_parameters.csv'));
writetable(tableOut(:, {'caseId', 'truthTipMode', 'decidedMode', 'estFpHz', 'errFpHz', 'errFpPct', ...
    'estAlphaDeg', 'errAlphaDeg', 'errAlphaPct', 'estSwingAmpDeg', ...
    'errSwingAmpDeg', 'estSwingFreqHz', 'errSwingFreqHz', 'modeCorrect'}), ...
    fullfile(outDir, 'identification_errors.csv'));

summary = struct();
summary.configuration = cfg;
summary.parsimonyF = cfg.parsimonyF;
summary.swingPromotionF = cfg.swingPromotionF;
summary.swingThresholdDeg = swingThresholdDeg;
summary.resultDirectory = resultDir;
summary.outputDirectory = outDir;
summary.table = tableOut;
summary.fits = fits;
summary.nScene = nScene;
summary.modeAccuracy = mean(tableOut.modeCorrect);
summary.fpRmsErrorHz = sqrt(mean(tableOut.errFpHz .^ 2));
summary.alphaRmsErrorDeg = sqrt(mean(tableOut.errAlphaDeg .^ 2));
summary.alphaMaxAbsErrorDeg = max(abs(tableOut.errAlphaDeg));
summary.swingAmpRmsErrorDeg = sqrt(mean(tableOut.errSwingAmpDeg .^ 2));
save(fullfile(outDir, 'identification_summary.mat'), 'summary');

fprintf('\n=== identified micro-motion parameters (blind) vs theoretical prior ===\n');
fprintf('%-22s %-18s %8s %8s %8s %8s\n', 'case', 'mode', 'fp err', 'alpha err', 'A err', 'fw err');
for s = 1:nScene
    r = tableOut(s, :);
    fprintf('%-22s %-18s %+8.4f %+8.4f %+8.4f %8s\n', char(r.caseId), ...
        char(r.decidedMode + "/" + r.truthTipMode), r.errFpHz, r.errAlphaDeg, ...
        r.errSwingAmpDeg, ternary(isnan(r.errSwingFreqHz), '-', ...
        sprintf('%+.4f', r.errSwingFreqHz)));
end
fprintf('mode accuracy %d/%d | fp rms %.4f Hz | alpha rms %.4f deg (max |e| %.4f) | A rms %.4f deg\n', ...
    sum(tableOut.modeCorrect), nScene, summary.fpRmsErrorHz, ...
    summary.alphaRmsErrorDeg, summary.alphaMaxAbsErrorDeg, summary.swingAmpRmsErrorDeg);

if opt(options, 'GenerateFigures', true)
    render_figures(tableOut, fits, outDir, swingThresholdDeg);
end
if opt(options, 'NoiseStudy', true)
    noiseTable = identify_isar_tip_noise_study(tableOut, fits, outDir, cfg, ...
        opt(options, 'NoiseTrials', 8));
    summary.noiseStudy = noiseTable;
end

result = summary;
result.table = tableOut;
fprintf('DONE %s\n', outDir);
end

% ======================================================================
function render_figures(T, fits, outDir, thresholdDeg)
font = 'Microsoft YaHei';
cTruth = [0.07 0.32 0.60];
cFit = [0.00 0.45 0.35];

% ---- figure 1: track fits -------------------------------------------------
pick = ["R1_s3_p3_a1_w0", "R2_s2_p3p5_a7_w0", "R3_s2_p3p5_a5_w3", "R3_s2_p3p5_a7_w5"];
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [40 40 2200 1450]);
layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:numel(pick)
    s = find(strcmp(string(T.caseId), pick(k)), 1);
    if isempty(s)
        continue;
    end
    f = fits{s};
    ax = nexttile(layout, k);
    hold(ax, 'on');
    tf = linspace(0, 5, 2001)';
    yPred = isar_tip_motion_doppler(tf, f.best.theta, isar_tip_motion_cfg(), ...
        f.best.theta(7) + f.best.theta(8) * (tf - mean(f.t)));
    truthTheta = [T.truthFpHz(s), T.truthAlphaDeg(s), T.truthSwingAmpDeg(s), ...
        T.truthSwingFreqHz(s), 0, 0, 0, 0];
    yTruth = isar_tip_motion_doppler(tf, truthTheta, isar_tip_motion_cfg(), ...
        f.best.theta(7) + f.best.theta(8) * (tf - mean(f.t)));
    plot(ax, tf, yTruth, '--', 'Color', cTruth, 'LineWidth', 2.2);
    plot(ax, tf, yPred, '-', 'Color', cFit, 'LineWidth', 2.4);
    plot(ax, f.t, f.y, 'o', 'Color', [0.85 0.33 0.10], 'MarkerFaceColor', 'w', ...
        'MarkerSize', 7, 'LineWidth', 1.2);
    grid(ax, 'on');
    xlim(ax, [0 5]);
    ylim(ax, [min(f.y) - 5, max(f.y) + 5]);
    title(ax, sprintf('%s：%s（A=%.2f°）', char(T.caseId(s)), ...
        prettyModel(T.decidedMode(s)), T.estSwingAmpDeg(s)), ...
        'FontSize', 26, 'FontWeight', 'normal', 'FontName', font, 'Interpreter', 'none');
    ax.FontSize = 24;
    ax.FontName = font;
    ax.LineWidth = 1.2;
    if k == 1
        legend(ax, {'物理真值', '盲识别模型', '提取峰值'}, 'Location', 'northoutside', ...
            'Orientation', 'horizontal', 'FontSize', 24, 'FontName', font);
    end
end
sgtitle(layout, '锥尖慢时间轨迹：盲识别模型与物理真值', 'FontSize', 30, ...
    'FontName', font, 'FontWeight', 'normal');
set(findall(fig, 'Type', 'text'), 'FontName', font);
exportgraphics(fig, fullfile(outDir, 'tip_mode_track_fits.png'), 'Resolution', 120);
savefig(fig, fullfile(outDir, 'tip_mode_track_fits.fig'));
close(fig);

% ---- figure 2: parameter accuracy ----------------------------------------
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [40 40 2200 1450]);
layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(layout, 1);
plot(ax, T.truthFpHz, T.estFpHz, 'o', 'MarkerSize', 10, 'MarkerFaceColor', cFit, ...
    'MarkerEdgeColor', 'none');
hold(ax, 'on');
lims = [0.5 4.5];
plot(ax, lims, lims, 'k--', 'LineWidth', 1.4);
grid(ax, 'on');
xlim(ax, lims); ylim(ax, lims);
title(ax, '进动频率 f_p (Hz)', 'FontSize', 26, 'FontName', font, 'FontWeight', 'normal');
ax.FontSize = 24; ax.FontName = font;

ax = nexttile(layout, 2);
plot(ax, T.truthAlphaDeg, T.estAlphaDeg, 'o', 'MarkerSize', 10, ...
    'MarkerFaceColor', cFit, 'MarkerEdgeColor', 'none');
hold(ax, 'on');
lims = [0 8];
plot(ax, lims, lims, 'k--', 'LineWidth', 1.4);
grid(ax, 'on');
xlim(ax, lims); ylim(ax, lims);
title(ax, '进动半锥角 α (°)', 'FontSize', 26, 'FontName', font, 'FontWeight', 'normal');
ax.FontSize = 24; ax.FontName = font;

ax = nexttile(layout, 3);
semilogy(ax, T.truthSwingAmpDeg + 0.01, T.estSwingAmpDeg + 0.01, 'o', ...
    'MarkerSize', 10, 'MarkerFaceColor', cFit, 'MarkerEdgeColor', 'none');
hold(ax, 'on');
plot(ax, [0.01 10], [0.01 10], 'k--', 'LineWidth', 1.4);
grid(ax, 'on');
xlim(ax, [0.01 10]); ylim(ax, [0.01 10]);
title(ax, '摆动幅度 A (°)，log 轴', 'FontSize', 26, 'FontName', font, ...
    'FontWeight', 'normal');
ax.FontSize = 24; ax.FontName = font;

ax = nexttile(layout, 4);
m = T.truthSwingFreqHz > 0;
plot(ax, T.truthSwingFreqHz(m), T.estSwingFreqHz(m), 'o', 'MarkerSize', 10, ...
    'MarkerFaceColor', cFit, 'MarkerEdgeColor', 'none');
hold(ax, 'on');
lims = [0.5 1.5];
plot(ax, lims, lims, 'k--', 'LineWidth', 1.4);
grid(ax, 'on');
xlim(ax, lims); ylim(ax, lims);
title(ax, '摆动频率 f_w (Hz)', 'FontSize', 26, 'FontName', font, 'FontWeight', 'normal');
ax.FontSize = 24; ax.FontName = font;

sgtitle(layout, '盲识别参数与理论先验对比（虚线为理想一致线）', 'FontSize', 30, ...
    'FontName', font, 'FontWeight', 'normal');
set(findall(fig, 'Type', 'text'), 'FontName', font);
exportgraphics(fig, fullfile(outDir, 'tip_mode_parameter_accuracy.png'), 'Resolution', 120);
savefig(fig, fullfile(outDir, 'tip_mode_parameter_accuracy.fig'));
close(fig);

% ---- figure 3: swing detection margin -------------------------------------
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [40 40 2200 1100]);
ax = axes(fig);
comp = max(T.swingComponentHz, 1e-4);
bar(ax, comp, 'FaceColor', [0.55 0.72 0.85], 'EdgeColor', 'none');
hold(ax, 'on');
% The decision threshold in Doppler terms: a swing of thresholdDeg at 1 Hz
% with no precession, so only the swing term contributes.
thresholdHz = isar_tip_motion_doppler((0:0.001:5)', ...
    [1, 0, thresholdDeg, 1, 0, 0, 0, 0], isar_tip_motion_cfg());
thresholdHz = max(abs(thresholdHz));   % single-sided amplitude of a sinusoid
yline(ax, thresholdHz, '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 2);
set(ax, 'YScale', 'log');
grid(ax, 'on');
ax.YMinorGrid = 'off';
xlim(ax, [0.4 height(T) + 0.6]);
ylim(ax, [1e-4 100]);
ax.XTick = 1:height(T);
ax.XTickLabel = cellstr(T.caseId);
ax.XTickLabelRotation = 40;
ax.FontSize = 20;
ax.FontName = font;
ylabel(ax, '拟合摆动分量幅度 (Hz)', 'FontSize', 24, 'FontName', font);
title(ax, sprintf('拟合摆动分量幅度（虚线为历史 %.2f° 参考值）', thresholdDeg), ...
    'FontSize', 26, 'FontName', font, 'FontWeight', 'normal');
exportgraphics(fig, fullfile(outDir, 'tip_mode_swing_margin.png'), 'Resolution', 120);
savefig(fig, fullfile(outDir, 'tip_mode_swing_margin.fig'));
close(fig);
end

% ======================================================================
function s = prettyModel(name)
switch lower(char(name))
    case 'precession'
        s = '纯进动';
    case 'precession+swing'
        s = '进动+摆动';
    case 'swing'
        s = '纯摆动';
    otherwise
        s = '静止';
end
end

function v = ternary(cond, a, b)
if cond
    v = a;
else
    v = b;
end
end

function v = opt(s, name, fallback)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = fallback;
end
end
