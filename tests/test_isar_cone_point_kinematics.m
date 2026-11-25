function test_isar_cone_point_kinematics()
% Conventions and invariants of the four-point cone kinematics.
% Guards the pieces that let the tail scatterers identify the spin rate:
% the tip must be spin-invariant, the base points must not be, and the
% analytic range rate must match the derivative of the analytic range.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'functions'));
cfg = isar_tip_motion_cfg();
P = isar_cone_body_points();

% geometry template of scatterermodel.createLegacyDiskTarget
expected = [0, -0.75, 0;
            0.3 * cos(pi / 6), 0.75, 0.3 * sin(pi / 6);
            0.3 * cos(pi / 6 + 2 * pi / 3), 0.75, 0.3 * sin(pi / 6 + 2 * pi / 3);
            0.3 * cos(pi / 6 + 4 * pi / 3), 0.75, 0.3 * sin(pi / 6 + 4 * pi / 3)];
assert(max(abs(P(:) - expected(:))) < 1e-12, 'Body point template changed.');

t = (0:0.001:1)';
% ---- 1. the tip is exactly spin-invariant and matches the tip model ------
for fs = [0, 1, 3.5]
    p = [fs, 1.4, 6.3, 4.1, 0.7, 0.8, -0.4, 0.3];
    [rGen, vGen] = isar_cone_point_kinematics(t, P(1, :), p, cfg);
    [rTip, vTip] = isar_tip_motion_kinematics(t, p(2), p(3), p(4), p(5), ...
        p(6), p(7), cfg);
    assert(max(abs(rGen - rTip)) < 1e-12, 'Tip range must ignore the spin rate.');
    assert(max(abs(vGen - vTip)) < 1e-9, 'Tip rate must ignore the spin rate.');
end

% ---- 2. the base points do see the spin rate -----------------------------
p0 = [0, 1.4, 6.3, 4.1, 0.7, 0.8, -0.4, 0.3];
p1 = p0;
p1(1) = 3;
[~, v0] = isar_cone_point_kinematics(t, P(2, :), p0, cfg);
[~, v1] = isar_cone_point_kinematics(t, P(2, :), p1, cfg);
deltaHz = max(abs(2 / cfg.lambdaM * (v1 - v0)));
assert(deltaHz > 100, ...
    'A base point must react strongly to fs, got %.3f Hz.', deltaHz);

% ---- 3. analytic rate equals the derivative of the analytic range -------
h = 1e-5;
for k = 1:4
    [~, v] = isar_cone_point_kinematics(t, P(k, :), p1, cfg);
    rp = isar_cone_point_kinematics(t + h, P(k, :), p1, cfg);
    rm = isar_cone_point_kinematics(t - h, P(k, :), p1, cfg);
    numeric = (rp - rm) / (2 * h);
    assert(max(abs(v - numeric)) < 1e-5, ...
        'Analytic rate mismatch for point %d: %.3e m/s.', k, ...
        max(abs(v - numeric)));
end

% ---- 4. dwell operator folds into the PRF interval ----------------------
p2 = [3.5, 3.5, 7, 5, 1, 0, 0, 0];
f = isar_multi_point_doppler(t, P, p2, cfg, true);
assert(all(abs(f(:)) <= cfg.prfHz / 2 + 1e-9), ...
    'Folded Doppler must stay inside +/-prfHz/2.');
fInst = isar_multi_point_doppler(t, P, p2, cfg, false);
assert(max(abs(fInst(:, 2))) > cfg.prfHz / 2, ...
    'This case is expected to exceed the unambiguous interval.');
% An independent sinusoidal midpoint sum checks tip folding and bias order.
theta = [4.8 25 0 0 0 0 8 2];
frame = 0.0156 + (0:49)' * cfg.sampleStepSec;
offset = ((1:8) - 4.5) * (0.032 / 8);
w = 2*pi*theta(1);
gain = 2/cfg.lambdaM * cfg.tipLeverM/sqrt(2) * sind(theta(2)) * w;
unfolded = gain * mean(sin(w*(frame + offset)), 2);
bias = theta(7) + theta(8)*(frame-mean(frame));
expected = mod(unfolded + bias + cfg.prfHz/2, cfg.prfHz) - cfg.prfHz/2;
assert(any(abs(unfolded + bias) > cfg.prfHz/2));
assert(max(abs(isar_tip_motion_doppler(frame, theta, cfg, bias) - expected)) < 1e-9);
fit = isar_fit_tip_motion_mode(frame, expected, 'H1', struct('cfg', cfg, 'seed', theta));
assert(fit.rmsHz < 1e-8, 'The fit residual must use the same folded observation.');
assert(isequal(isar_fold_doppler([-500 500], 1000), [-500 -500]));
fprintf('test_isar_cone_point_kinematics passed.\n');
end
