function result = identify_isar_spin_from_tail_points(options)
%IDENTIFY_ISAR_SPIN_FROM_TAIL_POINTS Blind spin identification from tail points.
%   RESULT = IDENTIFY_ISAR_SPIN_FROM_TAIL_POINTS(OPTIONS) extends the tip
%   analysis to all four cone scatterers. The tip P1 lies on the spin axis and
%   cannot see fs; the three base points P2..P4 at 0.3 m radius can. The
%   estimator is staged:
%
%     1. fit the tip track alone (isar_fit_tip_motion_mode, structure H3) to
%        obtain fp, alpha, A, fw, psiP and the nuisance c0 + c1;
%     2. hold those fixed and fit only [fs, psiS] against the four tracks,
%        using the true observation operator (32-pulse dwell spectrum peak,
%        with the PRF folding that aliasing produces).
%
%   The tip stage is safe at 1 kHz PRF: its Doppler never exceeds 171 Hz in
%   these scenarios. The base points reach about 131 Hz per Hz of spin plus
%   the precession contribution, so they leave +/-500 Hz above fs ~ 2 Hz and
%   fold.
%
%   OPTIONS: OutputDirectory, Cases (table), NoiseHz, NoiseSeed (synthetic
%   stream seed, default 20261021), GenerateFigures, TrackFile.

if nargin < 1 || isempty(options)
    options = struct();
end
root = fileparts(fileparts(mfilename('fullpath')));
outDir = char(opt(options, 'OutputDirectory', ...
    fullfile(root, 'results', 'isar_spin_tail_points_20261021')));
if ~isfolder(outDir)
    mkdir(outDir);
end
cfg = isar_tip_motion_cfg();
P = isar_cone_body_points();
noiseHz = opt(options, 'NoiseHz', 0.10);
useReal = isfield(options, 'TrackFile') && ~isempty(options.TrackFile);
realTracks = table();
if useReal
    % Real per-scatterer tracks from run_isar_multi_point_tracks. The scene
    % list and truth come from the sibling case table, and the spectrum
    % padding is raised to the chain's own 64*Na so the predicted peak uses
    % the same interpolation.
    trackFile = char(options.TrackFile);
    realTracks = readtable(trackFile, 'TextType', 'string');
    cases = readtable(fullfile(fileparts(trackFile), 'multi_point_cases.csv'), ...
        'TextType', 'string');
    if isfield(options, 'Cases') && ~isempty(options.Cases)
        cases = cases(ismember(string(cases.caseId), string(options.Cases)), :);
    end
    cfg.peakFftPoints = 64 * 32;
end
if ~useReal
    cases = opt(options, 'Cases', default_cases());
    previousRng = rng;
    restoreRng = onCleanup(@() rng(previousRng)); %#ok<NASGU>
    rng(opt(options, 'NoiseSeed', 20261021), 'twister');
end
nCase = height(cases);
fprintf('四散射点自旋识别：%d 个场景，PRF %g Hz，多普勒无模糊区间 ±%g Hz\n', ...
    nCase, cfg.prfHz, cfg.prfHz / 2);

rows = cell(nCase, 1);
scans = cell(nCase, 1);
tracks = cell(nCase, 1);
for c = 1:nCase
    q = cases(c, :);
    pTrue = [q.spinFrequencyHz, q.precessionFrequencyHz, q.precessionHalfAngleDeg, ...
        q.swingAmplitudeDeg, q.swingFrequencyHz, 0, 0, 0];
    if useReal
        sel = realTracks(string(realTracks.caseId) == string(q.caseId), :);
        sel = sortrows(sel, 'receiveMidpointSec');
        tTrack = sel.receiveMidpointSec;
        Y = [sel.peakP1Hz, sel.peakP2Hz, sel.peakP3Hz, sel.peakP4Hz];
        nuisance = [5.15, 0.14];
    else
        tTrack = 0.0156 + (0:49)' * cfg.sampleStepSec;
        nuisance = [5.15, 0.14];
        [Y, ~] = isar_simulate_cone_peaks(tTrack, P, pTrue, cfg, noiseHz, nuisance);
    end

    % ---- stage 1: tip only
    tipFit = isar_fit_tip_parsimonious(tTrack, Y(:, 1), cfg, struct('ranking', 'short'));
    th = tipFit.theta;
    pBase = [0, th(1), th(2), th(3), th(4), th(5), th(6), 0];
    c0 = th(7);
    c1 = th(8);

    % ---- stage 2: [fs, psiS, alpha, psiP] against all four tracks
    % Search cap from the leading base-tone approximation. Other motion
    % terms can still cross Nyquist inside this interval.
    baseGainHzPerHz = 2 / cfg.lambdaM * 0.3 / sqrt(2) * 2 * pi;
    fsCap = cfg.prfHz / (2 * baseGainHzPerHz) - th(1);
    fsCap = max(min(fsCap, 5.0), 0);
    spin = fitSpin(tTrack, Y, P, pBase, cfg, [c0, c1], fsCap);
    fsEst = spin.fs;
    pEst = pBase;
    pEst(1) = spin.fs;
    pEst(8) = spin.psiS;
    pEst(3) = spin.alpha;
    pEst(6) = spin.psiP;
    YHat = isar_predict_cone_peaks(tTrack, P, pEst, cfg, [c0, c1]);
    % rssSpin is the profiled objective, i.e. already with per-point bias
    % removed, which is the honest residual for the real tracks.
    residRms = sqrt(spin.rssSpin / numel(Y));

    % aliasing diagnostics on the unfolded base-point Doppler
    PInst = isar_multi_point_doppler(tTrack, P, pTrue, cfg, false);
    maxBaseHz = max(abs(PInst(:, 2:4)), [], 'all');
    foldCount = maxBaseHz / (cfg.prfHz / 2);

    rows{c} = struct('caseId', q.caseId, ...
        'truthSpinHz', q.spinFrequencyHz, 'estSpinHz', fsEst, ...
        'errSpinHz', fsEst - q.spinFrequencyHz, ...
        'truthFpHz', q.precessionFrequencyHz, 'estFpHz', th(1), ...
        'errFpHz', th(1) - q.precessionFrequencyHz, ...
        'truthAlphaDeg', q.precessionHalfAngleDeg, 'estAlphaDeg', th(2), ...
        'errAlphaDeg', th(2) - q.precessionHalfAngleDeg, ...
        'jointAlphaDeg', spin.alpha, ...
        'truthSwingAmpDeg', q.swingAmplitudeDeg, 'estSwingAmpDeg', th(3), ...
        'tipRmsHz', tipFit.rmsHz, 'jointRmsHz', residRms, ...
        'rssNoSpin', spin.rssNoSpin, 'rssSpin', spin.rssSpin, ...
        'Fstat', spin.F, 'spinResidualDf', spin.residualDf, ...
        'spinDetected', spin.detected, 'fsCapHz', fsCap, ...
        'maxBaseDopplerHz', maxBaseHz, 'foldCount', foldCount, ...
        'aliased', foldCount > 1);
    scans{c} = struct('fs', spin.fsGrid, 'rss', spin.rssGrid, ...
        'truth', q.spinFrequencyHz);
    tracks{c} = struct('t', tTrack, 'y', Y, 'pred', YHat, 'caseId', q.caseId, ...
        'pTrue', pTrue, 'pEst', pEst);

    fprintf(['[%d/%d] %-22s fs %5.2f (cap %4.2f) -> %6.3f Hz (err %+7.4f) | fp %5.2f -> %5.3f | ' ...
        'a %4.1f -> %5.2f | A %4.1f -> %5.2f | max base %6.0f Hz (%.2f x) %s | F %8.1f %s\n'], ...
        c, nCase, q.caseId, q.spinFrequencyHz, fsCap, fsEst, fsEst - q.spinFrequencyHz, ...
        th(1), q.precessionFrequencyHz, th(2), q.precessionHalfAngleDeg, ...
        th(3), q.swingAmplitudeDeg, maxBaseHz, foldCount, ...
        ternary(foldCount > 1, 'ALIASED', 'clean'), spin.F, ...
        ternary(spin.detected, 'DETECTED', 'rejected'));
end

T = struct2table([rows{:}]);
writetable(T, fullfile(outDir, 'spin_identification.csv'));
summary = struct('configuration', cfg, 'noiseHz', noiseHz, 'table', T, ...
    'scans', {scans}, 'tracks', {tracks}, 'bodyPoints', P);
save(fullfile(outDir, 'spin_identification.mat'), 'summary');

fprintf('\n=== 自旋识别（1 kHz PRF，四散射点联合）===\n');
fprintf('%-22s %8s %8s %9s %7s %8s %8s %7s %s\n', 'case', 'true fs', 'est fs', 'err', ...
    'fold x', 'tipRms', 'jointRms', 'F', 'verdict');
for c = 1:nCase
    r = T(c, :);
    fprintf('%-22s %8.2f %8.3f %+9.4f %7.2f %8.3f %8.3f %7.1f %s/%s\n', char(r.caseId), ...
        r.truthSpinHz, r.estSpinHz, r.errSpinHz, r.foldCount, r.tipRmsHz, ...
        r.jointRmsHz, r.Fstat, ternary(r.aliased, 'aliased', 'clean'), ...
        ternary(r.spinDetected, 'kept', 'rejected'));
end
fprintf('自旋频率 RMS 误差 %.4f Hz | 最大 |err| %.4f Hz | 判决正确 %d/%d\n', ...
    sqrt(mean(T.errSpinHz .^ 2)), max(abs(T.errSpinHz)), ...
    sum(T.spinDetected == (T.truthSpinHz > 0)), nCase);
if opt(options, 'GenerateFigures', true)
    render_figures(T, scans, tracks, outDir, cfg);
end
result = summary;
fprintf('DONE %s\n', outDir);
end

% ======================================================================
function out = fitSpin(tTrack, Y, P, pBase, cfg, nuisance, fsCap)
% Coarse (fs, psiS) grid, then a joint refinement of [fs psiS alpha psiP]
% against all four tracks, plus the nested no-spin hypothesis fs = 0.
%
% fs is scanned from 0, not from a positive floor: a "no spin" track can
% only be explained at fs = 0, and a positive lower bound made the earlier
% version return 4.95 Hz for a genuinely non-spinning scene.
% The cap limits the frequency search; foldCount is a simulation diagnostic.
fsUpper = fsCap;
fsGrid = (0:0.05:fsUpper)';
psiGrid = (0:pi / 6:2 * pi - pi / 6)';
rssGrid = zeros(numel(fsGrid), 1);
bestRss = inf;
bestV = [0, 0, pBase(3), pBase(6)];
for i = 1:numel(fsGrid)
    rowBest = inf;
    for j = 1:numel(psiGrid)
        p = pBase;
        p(1) = fsGrid(i);
        p(8) = psiGrid(j);
        YH = isar_predict_cone_peaks(tTrack, P, p, cfg, nuisance);
        r = Y - YH;
        r = r - mean(r, 1);   % insensitive to per-point Doppler offsets
        v = sum(r(:) .^ 2);
        if v < rowBest
            rowBest = v;
        end
        if v < bestRss
            bestRss = v;
            bestV = [fsGrid(i), psiGrid(j), pBase(3), pBase(6)];
        end
    end
    rssGrid(i) = sqrt(rowBest / numel(Y));
end
% joint refinement from the best grid point and from the tip-only solution
lb = [0, -pi, max(pBase(3) - 3, 0), -pi];
ub = [fsUpper, pi, pBase(3) + 3, pi];
obj = @(v) rssAt(v, tTrack, Y, P, pBase, cfg, nuisance, fsCap);
opt2 = optimset('MaxFunEvals', 1400, 'MaxIter', 1400, 'TolFun', 1e-12, ...
    'TolX', 1e-10, 'Display', 'off');
starts = [bestV; pBase(1), pBase(8), pBase(3), pBase(6)];
fAll = zeros(size(starts, 1), 1);
vAll = zeros(size(starts));
for s = 1:size(starts, 1)
    v = fminsearch(obj, starts(s, :), opt2);
    v = min(max(v, lb), ub);
    vAll(s, :) = v;
    fAll(s) = obj(v);
end
[~, ib] = min(fAll);
v = vAll(ib, :);
rssSpin = fAll(ib);
% nested no-spin hypothesis: fs pinned to zero, everything else refitted
obj0 = @(w) rssAt([0, w(1), w(2), w(3)], tTrack, Y, P, pBase, cfg, nuisance, fsCap);
w0 = fminsearch(obj0, [v(2), v(3), v(4)], opt2);
w0 = min(max(w0, lb(2:4)), ub(2:4));
rssNoSpin = obj0(w0);
n = numel(Y);
% Four nonlinear parameters and two profiled coefficients per track.
residualDf = n - (4 + 2 * size(Y, 2));
Fstat = (rssNoSpin - rssSpin) / max(rssSpin / residualDf, realmin);
out = struct('fs', v(1), 'psiS', v(2), 'alpha', v(3), 'psiP', v(4), 'fsCap', fsCap, ...
    'rssGrid', rssGrid, 'fsGrid', fsGrid, 'rssSpin', rssSpin, ...
    'rssNoSpin', rssNoSpin, 'F', Fstat, 'residualDf', residualDf, ...
    'detected', Fstat > cfg.parsimonyF && v(1) > 0.2 ...
        && sqrt(rssSpin / n) < cfg.spinFitRmsMaxHz);
end

function f = rssAt(v, tTrack, Y, P, pBase, cfg, nuisance, fsCap) %#ok<INUSD>
% Profiled objective: the per-point Doppler offset and drift are nuisance
% terms, fitted by least squares after unwrapping each point to the branch
% nearest the model prediction. The residual RVP/range-coupling bias scales
% with a scatterer's own range offset, so the three base points do not share
% the tip's bias: the real tracks differ from the tip by 8-13 Hz.
p = pBase;
p(1) = min(max(v(1), 0), fsCap);
p(8) = min(max(v(2), -pi), pi);
p(3) = min(max(v(3), 0), 25);
p(6) = min(max(v(4), -pi), pi);
M = isar_predict_cone_peaks(tTrack, P, p, cfg, [0 0]);
A2 = [ones(numel(tTrack), 1), tTrack - mean(tTrack)];
total = 0;
for k = 1:size(M, 2)
    d = Y(:, k) - M(:, k);
    d = isar_fold_doppler(d, cfg.prfHz);
    ab = A2 \ d;
    r = d - A2 * ab;
    total = total + r' * r;
end
f = total;
end

function cases = default_cases()
% Priority case first: 1 Hz spin with 1 Hz precession stays inside the
% unambiguous Doppler interval. The rest increase the spin/precession rate
% until the base scatterers fold, then keep going.
id = ["S1_fs1_fp1_a3"; "S2_fs1_fp1_a1"; "S3_fs1_fp3_a3"; ...
    "S4_fs3_fp3_a3"; "S5_fs3p5_fp3p5_a7w5"; "S6_fs0_fp3_a3"];
fs = [1; 1; 1; 3; 3.5; 0];
fp = [1; 1; 3; 3; 3.5; 3];
alpha = [3; 1; 3; 3; 7; 3];
amp = [0; 0; 0; 0; 5; 0];
fw = [0; 0; 0; 0; 1; 0];
cases = table(id, fs, fp, alpha, amp, fw, 'VariableNames', ...
    {'caseId', 'spinFrequencyHz', 'precessionFrequencyHz', ...
    'precessionHalfAngleDeg', 'swingAmplitudeDeg', 'swingFrequencyHz'});
end

function render_figures(T, scans, tracks, outDir, cfg)
font = 'Microsoft YaHei';
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [40 40 2300 1400]);
layout = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for c = 1:height(T)
    ax = nexttile(layout, c);
    semilogy(ax, scans{c}.fs, scans{c}.rss, '-', 'Color', [0.00 0.45 0.35], ...
        'LineWidth', 2.4);
    hold(ax, 'on');
    xline(ax, scans{c}.truth, '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 2);
    xline(ax, T.estSpinHz(c), ':', 'Color', [0.07 0.32 0.60], 'LineWidth', 2);
    grid(ax, 'on');
    ax.FontSize = 18;
    ax.FontName = font;
    xlim(ax, [0 5]);
    title(ax, sprintf('%s（折叠 %.1f×）', char(T.caseId(c)), T.foldCount(c)), ...
        'FontSize', 20, 'FontName', font, 'FontWeight', 'normal', 'Interpreter', 'none');
    if c == 1
        legend(ax, {'目标函数', '真值 f_s', '估计 f_s'}, 'FontSize', 16, ...
            'FontName', font, 'Location', 'best');
    end
end
sgtitle(layout, '自旋频率扫描：橙色真值、蓝色估计，纵轴为四轨联合残差 (Hz)', ...
    'FontSize', 26, 'FontName', font, 'FontWeight', 'normal');
exportgraphics(fig, fullfile(outDir, 'spin_scan_objective.png'), 'Resolution', 120);
close(fig);

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [40 40 2200 1200]);
layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for c = 1:min(4, height(T))
    ax = nexttile(layout, c);
    hold(ax, 'on');
    cols = lines(4);
    tr = tracks{c};
    for k = 1:4
        plot(ax, tr.t, tr.y(:, k), '.', 'Color', cols(k, :), 'MarkerSize', 8);
        plot(ax, tr.t, tr.pred(:, k), '-', 'Color', cols(k, :), 'LineWidth', 1.6);
    end
    grid(ax, 'on');
    xlim(ax, [0 5]);
    ax.FontSize = 18;
    ax.FontName = font;
    title(ax, sprintf('%s：实测(点) / 预测(线)', char(T.caseId(c))), ...
        'FontSize', 20, 'FontName', font, 'FontWeight', 'normal', 'Interpreter', 'none');
    if c == 1
        legend(ax, {'P1', 'P1 模型', 'P2', 'P2 模型', 'P3', 'P3 模型', 'P4', 'P4 模型'}, ...
            'FontSize', 14, 'FontName', font, 'NumColumns', 2, 'Location', 'best');
    end
end
sgtitle(layout, '四散射点慢时间峰值轨迹与联合模型', 'FontSize', 26, ...
    'FontName', font, 'FontWeight', 'normal');
exportgraphics(fig, fullfile(outDir, 'spin_tracks_vs_model.png'), 'Resolution', 120);
close(fig);
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
