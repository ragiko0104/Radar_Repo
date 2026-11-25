function test_isar_tip_motion_mode_identification()
% Convention, spin invisibility and blind-recovery checks for the tip mode fit.
% The closed form in isar_tip_motion_kinematics is compared against an
% independent rotation-matrix evaluation of the library attitude convention,
% and the estimator is exercised on synthetic tracks whose parameters are off
% every search grid, so a pass cannot come from grid alignment.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'functions'), fullfile(root, 'experiments'));
cfg = isar_tip_motion_cfg();

% ---- 1. closed form vs explicit rotation matrices -------------------------
pairs = [2.7, 6.3, 4.1, 0.7, 0.8, -0.4;
         3.0, 1.0, 0.0, 1.0, 0.0, 0.0;
         1.4, 12.0, 2.5, 2.2, -2.0, 1.1];
t = (0:0.001:1)';
for r = 1:size(pairs, 1)
    p = pairs(r, :);
    [refRange, refRate] = reference_tip(t, p(1), p(2), p(3), p(4), p(5), p(6));
    [rangeM, rateMps] = isar_tip_motion_kinematics(t, p(1), p(2), p(3), p(4), p(5), p(6), cfg);
    assert(max(abs(rangeM - refRange)) < 1e-12, 'Range convention mismatch.');
    assert(max(abs(rateMps - refRate)) < 1e-5, 'Range-rate convention mismatch.');
end

% ---- 2. the tip cannot see the spin rate ---------------------------------
% Compare reference trajectories evaluated with different spin rates against
% each other, so the check is exact rather than limited by the finite
% difference used for the reference range rate.
[refRange0, refRate0] = reference_tip(t, 2.7, 6.3, 4.1, 0.7, 0.8, -0.4, 0);
for fs = [0.7, 3, 7.5]
    [refRange, refRate] = reference_tip(t, 2.7, 6.3, 4.1, 0.7, 0.8, -0.4, fs);
    assert(max(abs(refRange - refRange0)) < 1e-12, ...
        'Spin must not move the tip, fs=%g.', fs);
    assert(max(abs(refRate - refRate0)) < 1e-12, ...
        'Spin must not change the tip range rate, fs=%g.', fs);
    [~, rateMps] = isar_tip_motion_kinematics(t, 2.7, 6.3, 4.1, 0.7, 0.8, -0.4, cfg);
    assert(max(abs(rateMps - refRate0)) < 1e-5, ...
        'Axisymmetric tip must reduce to the spin-free closed form.');
end

% ---- 3. blind recovery on a synthetic composite track --------------------
rng(20261021);
truth = [2.7, 6.3, 4.1, 0.7, 0.8, -0.4, 5.0, 0.14];
tTrack = 0.0156 + (0:49)' * 0.1;
y = synthetic_track(tTrack, truth, 0.5, cfg);
fits = fit_bank(tTrack, y, cfg);
[~, chosen] = min([fits.rss]);
assert(strcmp(fits(chosen).model, 'H3'), ...
    'Composite track should select H3, got %s.', fits(chosen).model);
th = fits(chosen).theta;
assert(abs(th(1) - truth(1)) < 0.02, 'fp error %.4f Hz.', th(1) - truth(1));
assert(abs(th(2) - truth(2)) < 0.10, 'alpha error %.4f deg.', th(2) - truth(2));
assert(abs(th(3) - truth(3)) < 0.10, 'swing amplitude error %.4f deg.', th(3) - truth(3));
assert(abs(th(4) - truth(4)) < 0.02, 'swing frequency error %.4f Hz.', th(4) - truth(4));

% ---- 4. a precession-only track must not invent a swing ------------------
truth2 = [1.7, 4.4, 0.0, 1.0, 0.0, 0.0, 5.0, 0.14];
y2 = synthetic_track(tTrack, truth2, 0.3, cfg);
fits2 = fit_bank(tTrack, y2, cfg);
th2 = fits2(2).theta;   % H1, the physically correct structure
assert(abs(th2(1) - truth2(1)) < 0.01, 'fp error on precession-only track.');
assert(abs(th2(2) - truth2(2)) < 0.05, 'alpha error on precession-only track.');
best2 = fits2([fits2.rss] == min([fits2.rss]));
if numel(best2) ~= 1
    [~, iBest2] = min([fits2.rss]);
    best2 = fits2(iBest2);
end
assert(best2.theta(3) < cfg.swingDecisionDeg, ...
    'Precession-only track must not exceed the swing decision threshold, got %.3f deg.', ...
    best2.theta(3));
assert(best2.swingComponentHz < 1.0, 'Spurious swing component too large.');

% ---- 5. H1 phase regression, including both sides of the phase seam -------
P = isar_cone_body_points();
phaseCases = zeros(6, 5);
row = 0;
for alpha = [8 20]
    for phase = [-pi + 0.02, 0, pi - 0.02]
        rng(123);
        yPhase = isar_simulate_cone_peaks(tTrack, P(1,:), ...
            [1 0.5 alpha 0 0 phase 0 0], cfg, 0.2, [5.15 0.14]);
        id = isar_identify_tip_mode(tTrack, yPhase, cfg);
        row = row + 1;
        phaseCases(row,:) = [alpha, phase, id.Fstat, id.fits{2}.rss, ...
            id.model == "H1" && id.precession && ~id.swing && ...
            abs(id.selected.theta(1) - 0.5) < 0.01 && ...
            abs(id.selected.theta(2) - alpha) < 0.1];
    end
end
disp(array2table(phaseCases, 'VariableNames', {'alpha','phase','F3','rssH1','ok'}));
assert(all(phaseCases(:,5)), 'Pure precession phase regression failed.');
assert(phaseCases(end,4) < 1.131, 'H1 discarded a feasible blind-grid fit.');

% The overlapping 6 degree setting gives identical measurements for both spins.
rng(42);
pair = [3 0.5 6 0 0 0.4 0 0];
yFast = isar_simulate_cone_peaks(tTrack, P(1,:), pair, cfg, 0.2, [5.15 0.14]);
pair(1) = 0.8;
rng(42);
ySlow = isar_simulate_cone_peaks(tTrack, P(1,:), pair, cfg, 0.2, [5.15 0.14]);
assert(isequal(yFast, ySlow), 'Axial spin changed a tip measurement.');
idPair = isar_identify_tip_mode(tTrack, yFast, cfg);
assert(idPair.model == "H1" && abs(idPair.selected.theta(2) - 6) < 0.1);

% ---- 6. a pure-swing track must be labelled swing, not precession+swing --
% Regression guard for the F4 downgrade gate. Without it a swing at fsw = fp
% is absorbed by the H1 tone, F1 fires and the mode came out as
% "precession+swing" with alpha = 0.004 deg.
% truth layout is [fp alphaDeg ampDeg fsw psiP psiW c0 c1].
rng(20261027);
swingTruth = [1.7, 0.0, 5.0, 1.0, 0.0, 0.0, 5.0, 0.14];
y3 = synthetic_track(tTrack, swingTruth, 0.2, cfg);
id3 = isar_identify_tip_mode(tTrack, y3, cfg);
assert(id3.mode == "swing", ...
    'Pure-swing track must be labelled swing, got %s (F1 %.1f, F3 %.1f, F4 %.1f).', ...
    id3.mode, id3.ladder.Fprecession, id3.ladder.FswingGivenPrecession, ...
    id3.ladder.FprecessionGivenSwing);
assert(abs(id3.selected.theta(3) - 5.0) < 0.2, ...
    'Swing amplitude on a pure-swing track error %.3f deg.', ...
    id3.selected.theta(3) - 5.0);
assert(id3.swing && ~id3.precession && id3.model == "H2");
for n = [17 18]
    tShort = tTrack(1:n);
    rng(20261027);
    short = isar_identify_tip_mode(tShort, synthetic_track(tShort, swingTruth, 0.2, cfg), cfg);
    assert(short.ladder.downgradeEnabled == (n >= 18));
    assert(short.Fstat > cfg.parsimonyF && short.ladder.FprecessionGivenSwing <= cfg.parsimonyF);
    assert(short.swing && short.precession == (n == 17));
end
try
    isar_identify_tip_mode(tTrack(1:8), y3(1:8), cfg);
    error('Test:ExpectedError', 'Expected the zero residual df input to fail.');
catch ex
    assert(strcmp(ex.identifier, 'RadarProject:TipMode:InsufficientSamples'));
end
% Historical pure-precession extraction with the largest spurious F3 (36.08).
% Only the two columns consumed by the production classifier are retained.
track = readtable(fullfile(root, 'tests', 'fixtures', 'isar_tip_R2_s2_p2_a5_w0.csv'));
conservative = cfg;
conservative.swingPromotionF = 40;
idHistorical = isar_identify_tip_mode(track.receiveMidpointSec, track.extractedTipPeakHz, conservative);
assert(idHistorical.model == "H1" && ~idHistorical.swing);
check_experiment_scoring(tTrack, y2, y3, y);
fprintf('test_isar_tip_motion_mode_identification passed.\n');
end

function check_experiment_scoring(t, precession, swing, composite)
% Deliberately mismatched truth rows exercise full-mode scoring at the CSV entry.
folder = tempname;
mkdir(fullfile(folder, 'group1'));
cleanup = onCleanup(@() rmdir(folder, 's')); %#ok<NASGU>
caseId = ["R1_a"; "R1_b"; "R1_c"];
truth = table(caseId, ones(3,1), ["spin"; "combined"; "precession"], ...
    [0; 1.7; 2.7], [0; 3; 6.3], [0; 5; 0], [0; 1; 0], ...
    'VariableNames', {'caseId','family','mode','precessionFrequencyHz', ...
    'precessionHalfAngleDeg','swingAmplitudeDeg','swingFrequencyHz'});
writetable(truth, fullfile(folder, 'all_scene_similarity.csv'));
tracks = table(repelem(caseId, numel(t)), repmat(t,3,1), [precession; swing; composite], ...
    'VariableNames', {'caseId','receiveMidpointSec','extractedTipPeakHz'});
writetable(tracks, fullfile(folder, 'group1', 'tip_extracted_tracks.csv'));
result = identify_isar_tip_micro_motion(struct('ResultDirectory', folder, ...
    'OutputDirectory', folder, 'GenerateFigures', false, 'NoiseStudy', false));
assert(isequal(result.table.truthTipMode, ["static"; "precession+swing"; "precession"]));
assert(isequal(result.table.decidedMode, ["precession"; "swing"; "precession+swing"]));
assert(result.swingPromotionF == result.parsimonyF);
assert(~any(result.table.modeCorrect) && result.modeAccuracy == 0);
assert(isfile(fullfile(folder, 'identification_errors.csv')));
strict = identify_isar_tip_micro_motion(struct('ResultDirectory', folder, ...
    'OutputDirectory', fullfile(folder, 'strict'), ...
    'GenerateFigures', false, 'NoiseStudy', false, 'SwingPromotionF', realmax));
assert(strict.swingPromotionF == realmax && strict.parsimonyF == result.parsimonyF);
assert(isequal(strict.table.selectedModel, ["H1"; "H2"; "H1"]));
end

% ======================================================================
function fits = fit_bank(t, y, cfg)
names = {'H0', 'H1', 'H2', 'H3'};
fits = cell(numel(names), 1);
for k = 1:numel(names)
    fits{k} = isar_fit_tip_motion_mode(t, y, names{k}, struct('cfg', cfg));
end
fits = [fits{:}];
end

function y = synthetic_track(t, truth, noiseHz, cfg)
rngFree = randn(numel(t), 1);
f = isar_tip_motion_doppler(t, truth, cfg);
y = f + truth(7) + truth(8) * (t - mean(t)) + noiseHz * rngFree;
end

function [rangeM, rateMps] = reference_tip(t, fp, alphaDeg, ampDeg, fw, psiP, psiW, fs)
% Independent evaluation: assemble the attitude from explicit rotations.
%   C(t) = C0 * Rz(beta) * Ry(2*pi*fp*t + psiP) * Rz(alpha) * Ry(2*pi*fs*t)
%   beta = -A*cos(2*pi*fw*t + psiW), C0 = Rz(-135 deg), LOS = +X
% The tip sits at (0,-0.75,0) in the body frame.
if nargin < 8
    fs = 0;
end
lever = 0.75;
b = [0; -lever; 0];
phi0 = deg2rad(-135);
C0 = [cos(phi0), -sin(phi0), 0; sin(phi0), cos(phi0), 0; 0, 0, 1];
a = deg2rad(alphaDeg);
A = deg2rad(ampDeg);
rangeM = zeros(numel(t), 1);
rateMps = zeros(numel(t), 1);
h = 1e-5;
for k = 1:numel(t)
    rangeM(k) = los(t(k));
    rateMps(k) = (los(t(k) + h) - los(t(k) - h)) / (2 * h);
end
    function r = los(tau)
        beta = -A * cos(2 * pi * fw * tau + psiW);
        thp = 2 * pi * fp * tau + psiP;
        ths = 2 * pi * fs * tau;
        C = C0 * rotZ(beta) * rotY(thp) * rotZ(a) * rotY(ths);
        p = C * b;
        r = p(1);
    end
end

function R = rotZ(a)
R = [cos(a), -sin(a), 0; sin(a), cos(a), 0; 0, 0, 1];
end

function R = rotY(a)
R = [cos(a), 0, sin(a); 0, 1, 0; -sin(a), 0, cos(a)];
end
