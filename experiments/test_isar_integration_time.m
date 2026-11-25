function out = test_isar_integration_time(outputDirectory)
%TEST_ISAR_INTEGRATION_TIME Effect of 5 s -> 2 s integration on the current modes.
%   Part 1: the six current scenes, identified from 5 s (50 frames) and from
%   2 s (21 frames) tracks, five trials each.
%   Part 2: the precession-frequency lower bound at 2 s, since 1/T doubles
%   from 0.2 Hz to 0.5 Hz.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'functions'));
if nargin < 1
    outputDirectory = fullfile(root, 'results', 'isar_integration_time_20261027');
end
if ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end
cfg = isar_tip_motion_cfg();
P = isar_cone_body_points();
P1 = P(1, :);
sigma = 0.2;
nuisance = [5.15 0.14];
nTrial = 5;
scenes = struct( ...
    'tag', {'S1_a1', 'S1_a3', 'S2_a3', 'S2_a5', 'S3_w3', 'S3_w5'}, ...
    'fs', {3, 3, 3, 3, 3, 3}, 'fp', {3, 3, 3.5, 3.5, 3.5, 3.5}, ...
    'alpha', {1, 3, 3, 5, 5, 5}, 'A', {0, 0, 0, 0, 3, 5}, ...
    'fw', {0, 0, 0, 0, 1, 1});
durations = [5, 2];
rows = struct('scene', {}, 'T', {}, 'trial', {}, 'mode', {}, 'fp', {}, ...
    'alpha', {}, 'A', {}, 'fw', {}, 'rmsHz', {});
for s = 1:numel(scenes)
    sc = scenes(s);
    p = [sc.fs, sc.fp, sc.alpha, sc.A, sc.fw, 0, 0, 0];
    for T = durations
        n = round(T / 0.1);
        t = 0.0156 + (0:n - 1)' * 0.1;
        for k = 1:nTrial
            rng(9000 + 10 * s + k);
            Y = isar_simulate_cone_peaks(t, P1, p, cfg, sigma, nuisance);
            id = isar_identify_tip_mode(t, Y(:, 1), cfg);
            th = id.selected.theta;
            rows(end + 1) = struct('scene', sc.tag, 'T', T, 'trial', k, ...
                'mode', id.mode, 'fp', th(1), 'alpha', th(2), 'A', th(3), ...
                'fw', th(4), 'rmsHz', id.selected.rmsHz); %#ok<AGROW>
        end
    end
end
tbl = struct2table(rows);
writetable(tbl, fullfile(outputDirectory, 'integration_time_scenes.csv'));
fprintf('=== part 1: six current scenes ===\n');
fprintf('%-8s %3s %10s %10s %10s %10s %10s %6s\n', 'scene', 'T', 'fp err', ...
    'alpha err', 'A err', 'fw err', 'rms', 'ok');
for s = 1:numel(scenes)
    sc = scenes(s);
    for T = durations
        m = string(tbl.scene) == sc.tag & tbl.T == T;
        sub = tbl(m, :);
        ok = mean(sub.mode == string(mode_of(sc)));
        fprintf('%-8s %3.0f %10.4f %10.4f %10.4f %10.4f %10.3f %6.2f\n', ...
            sc.tag, T, mean(abs(sub.fp - sc.fp)), mean(abs(sub.alpha - sc.alpha)), ...
            mean(abs(sub.A - sc.A)), mean(abs(sub.fw - sc.fw)), ...
            mean(sub.rmsHz), ok);
    end
end

% ------------------------------------------------------------------ part 2
fprintf('\n=== part 2: fp lower bound at T = 2 s (alpha 3 deg, no swing) ===\n');
fprintf('%8s %8s %10s %10s %10s %6s\n', 'fp', 'T', 'est mean', 'rms err', 'resid', 'P');
for fp = [0.3 0.5 0.75 1.0 1.5 2.0]
    for T = [5 2]
        n = round(T / 0.1);
        t = 0.0156 + (0:n - 1)' * 0.1;
        est = zeros(nTrial, 1);
        res = zeros(nTrial, 1);
        nOk = 0;
        for k = 1:nTrial
            rng(9500 + 100 * round(fp * 100) + 10 * T + k);
            Y = isar_simulate_cone_peaks(t, P1, [0, fp, 3, 0, 0.7, 0, 0, 0], cfg, sigma, nuisance);
            id = isar_identify_tip_mode(t, Y(:, 1), cfg);
            est(k) = id.selected.theta(1);
            res(k) = id.selected.rmsHz;
            nOk = nOk + (abs(est(k) - fp) < 0.1 * fp);
        end
        fprintf('%8.2f %8.0f %10.3f %10.4f %10.3f %6.2f\n', fp, T, mean(est), ...
            sqrt(mean((est - fp) .^ 2)), mean(res), nOk / nTrial);
    end
end
out = struct('scenes', tbl, 'sigma', sigma);
save(fullfile(outputDirectory, 'integration_time.mat'), 'out');
fprintf('DONE %s\n', outputDirectory);
end

function m = mode_of(sc)
if sc.A > 0
    m = 'precession+swing';
else
    m = 'precession';
end
end
