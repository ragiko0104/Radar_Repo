function out = isar_range_fit_noise(options)
%ISAR_RANGE_FIT_NOISE Acceleration from range measurements only.
%   The tracker supplies R(t); no image peak and no labelled scatterer
%   channel is used. A batch least-squares fit of R = r0 + v0 t + 0.5 a t^2
%   gives a_hat.
%
%   TWO TRACKER MODELS, and only one of them is the nominal case:
%
%     track="centre"  the tracker reports the BODY range. Micro-motion is a
%       relative motion of the scatterers about that centre, so it never
%       enters R; the fit sees measurement noise only (about 0.02 m/s^2 at
%       2 s). This is the nominal architecture, and the one used by
%       ISAR_MULTI_FRAME_10TARGETS_V2.
%
%     track="tip"     the tracker locks onto the single strongest scatterer,
%       so R carries that scatterer's micro-motion (+-0.36 m at alpha 20 deg)
%       and the quadratic fit absorbs part of it. The error rises to
%       0.1-0.8 m/s^2. This is a SENSITIVITY UPPER BOUND for a point-locked
%       tracker, NOT the nominal figure of merit.
%
%   Both are far below the 10.55 m/s^2 tolerance, but the distinction matters
%   when quoting a number: the nominal is noise-limited only.

if nargin < 1
    options = struct();
end
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'functions'));
outDir = char(opt(options, 'OutputDirectory', ...
    fullfile(root, 'results', 'isar_range_fit_noise_20261123')));
if ~isfolder(outDir)
    mkdir(outDir);
end
cfg = isar_tip_motion_cfg();
P = isar_cone_body_points();
lam = 3e8 / 14.8e9;
Ta = 1e-3;
nPulse = 32;
T = (nPulse - 1) * Ta;
aCrit = lam / (2 * T^2);
aTrue = -15 * 10;   % radial projection of the golden -15 g target
v0 = 10.2 * 340;    % radial projection of the 10.2 Mach target
R0 = 10050;
sigR = 0.02;
nSeed = 10;
fprintf('sigma_R %.3f m, a_crit %.2f m/s2, a_true %.1f\n', sigR, aCrit, aTrue);
fprintf('%4s %5s %5s %10s %10s %10s %9s\n', 'Tobs', 'alpha', 'fp', ...
    'max|da|', 'std da', 'k=max/a_crit', 'mean da');
rows = struct('Tobs', {}, 'alpha', {}, 'fp', {}, 'track', {}, ...
    'maxDa', {}, 'stdDa', {}, 'meanDa', {});
for track = ["tip", "centre"]
    for Tobs = [2 4 5]
        t = (0:0.1:Tobs)';
        for alpha = [10 20]
            for fp = [0.3 0.5 1 2]
                p = [0, fp, alpha, 0, 0.7, 0, 0, 0];
                [Rm, ~] = isar_cone_point_kinematics(t, P(1, :), p, cfg);
                if track == "tip"
                    rMicro = Rm;
                else
                    rMicro = zeros(size(t));
                end
                Rtrue = R0 + v0 * t + 0.5 * aTrue * t .^ 2 + rMicro;
                X = [ones(numel(t), 1), t, 0.5 * t .^ 2];
                da = zeros(nSeed, 1);
                for s = 1:nSeed
                    rng(3000 + 100 * Tobs + 10 * alpha + s);
                    th = X \ (Rtrue + sigR * randn(size(Rtrue)));
                    da(s) = th(3) - aTrue;
                end
                rows(end + 1) = struct('Tobs', Tobs, 'alpha', alpha, 'fp', fp, ...
                    'track', track, 'maxDa', max(abs(da)), 'stdDa', std(da), ...
                    'meanDa', mean(da)); %#ok<AGROW>
                if track == "tip"
                    fprintf('%4.0f %5.0f %5.1f %10.3f %10.3f %10.3f %9.3f\n', ...
                        Tobs, alpha, fp, max(abs(da)), std(da), ...
                        max(abs(da)) / aCrit, mean(da));
                end
            end
        end
    end
end
tbl = struct2table(rows);
writetable(tbl, fullfile(outDir, 'range_fit_noise.csv'));
save(fullfile(outDir, 'range_fit_noise.mat'), 'tbl', 'aCrit', 'aTrue', 'sigR');
fprintf('\n--- body-centre control (no micro-motion) ---\n');
for Tobs = [2 4 5]
    m = tbl.Tobs == Tobs & tbl.track == "centre";
    fprintf('T=%.0f s: max|da| %.3f m/s2, k %.3f\n', Tobs, ...
        max(tbl.maxDa(m)), max(tbl.maxDa(m)) / aCrit);
end
fprintf('\nworst over all cases: %.3f m/s2 (k = %.3f)\n', max(tbl.maxDa), ...
    max(tbl.maxDa) / aCrit);
fprintf('DONE %s\n', outDir);
out = struct('table', tbl);
end

function v = opt(s, name, fallback)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = fallback;
end
end
