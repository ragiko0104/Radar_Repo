function out = test_isar_integration_alpha_grid(outputDirectory)
%TEST_ISAR_INTEGRATION_ALPHA_GRID 2 s / 4 s sufficiency for alpha 3-20 deg.
%   The binding quantity at short integration is the tip modulation amplitude
%   328.7*sin(alpha)*fp against the track noise, not alpha alone, so the grid
%   sweeps alpha and fp together for T = 2 s and T = 4 s (T = 5 s is the
%   reference that already passed the sensitivity sweep).

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'functions'));
if nargin < 1
    outputDirectory = fullfile(root, 'results', 'isar_integration_alpha_grid_20261116');
end
if ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end
cfg = isar_tip_motion_cfg();
P = isar_cone_body_points();
P1 = P(1, :);
sigma = 0.2;
nuisance = [5.15 0.14];
nTrial = 3;
alphas = [3 5 8 12 20];
fps = [0.3 0.5 1.0 2.0 3.5];
durations = [2 4 5];
rows = struct('T', {}, 'alpha', {}, 'fp', {}, 'trial', {}, 'mode', {}, ...
    'estFp', {}, 'estAlpha', {}, 'fpOk', {}, 'alphaOk', {}, 'modeOk', {});
for T = durations
    n = round(T / 0.1);
    t = 0.0156 + (0:n - 1)' * 0.1;
    for a = alphas
        for fp = fps
            for k = 1:nTrial
                rng(11000 + 1000 * T + 100 * a + 10 * round(fp * 10) + k);
                Y = isar_simulate_cone_peaks(t, P1, [0, fp, a, 0, 0.7, 0, 0, 0], ...
                    cfg, sigma, nuisance);
                id = isar_identify_tip_mode(t, Y(:, 1), cfg);
                th = id.selected.theta;
                rows(end + 1) = struct('T', T, 'alpha', a, 'fp', fp, 'trial', k, ...
                    'mode', id.mode, 'estFp', th(1), 'estAlpha', th(2), ...
                    'fpOk', double(abs(th(1) - fp) < 0.1 * fp), ...
                    'alphaOk', double(abs(th(2) - a) < 0.1 * a), ...
                    'modeOk', double(id.mode == "precession")); %#ok<AGROW>
            end
        end
    end
end
tbl = struct2table(rows);
writetable(tbl, fullfile(outputDirectory, 'alpha_grid.csv'));
fprintf('=== pass rate over 3 trials, alpha 3-20 deg, fp 0.3-3.5 Hz ===\n');
fprintf('%3s %7s %8s %8s %8s %8s\n', 'T', 'alpha', 'fp', 'P(fp)', 'P(alpha)', 'P(mode)');
for T = durations
    for a = alphas
        line = '';
        for fp = fps
            m = tbl.T == T & tbl.alpha == a & tbl.fp == fp;
            sub = tbl(m, :);
            line = [line sprintf('  %.2f/%.2f/%.2f', mean(sub.fpOk), ...
                mean(sub.alphaOk), mean(sub.modeOk))]; %#ok<AGROW>
        end
        fprintf('%3.0f %7.0f %8s %s\n', T, a, '', line);
    end
end
fprintf('\nfp columns: %s\n', sprintf('%.1f  ', fps));
fprintf('\n=== worst cases (P(fp) < 1 or P(alpha) < 1) ===\n');
for T = durations
    for a = alphas
        for fp = fps
            m = tbl.T == T & tbl.alpha == a & tbl.fp == fp;
            sub = tbl(m, :);
            if mean(sub.fpOk) < 1 || mean(sub.alphaOk) < 1
                fprintf('T %3.0f s | alpha %4.0f deg | fp %4.1f Hz | fp ok %.2f | alpha ok %.2f | mode ok %.2f | dAlpha %6.2f\n', ...
                    T, a, fp, mean(sub.fpOk), mean(sub.alphaOk), mean(sub.modeOk), ...
                    mean(abs(sub.estAlpha - a)));
            end
        end
    end
end
out = struct('table', tbl, 'alphas', alphas, 'fps', fps, 'durations', durations);
save(fullfile(outputDirectory, 'alpha_grid.mat'), 'out');
fprintf('DONE %s\n', outputDirectory);
end
