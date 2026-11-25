function out = test_isar_micro_motion_resolution(outputDirectory)
%TEST_ISAR_MICRO_MOTION_RESOLUTION Upper bounds to 25 deg and parameter resolution.
%   Part 1: push the alpha and swing-amplitude sweeps up to 25 deg, the search
%   bound, to check whether the previous "12/10 deg" ceilings were real or just
%   the edge of the test grid.
%
%   Part 2: measure resolution the way it is normally defined - perturb one
%   parameter by a small step delta, refit, and ask whether the two values can
%   be told apart. The criterion is
%
%       |est - (p0 + delta)| < delta/2
%
%   i.e. the estimate must land closer to the perturbed value than to the
%   original one, with margin. The smallest delta that succeeds in every trial
%   is reported as the resolution.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'functions'));
if nargin < 1
    outputDirectory = fullfile(root, 'results', 'isar_resolution_20261027');
end
if ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end
cfg = isar_tip_motion_cfg();
P = isar_cone_body_points();
P = P(1, :);
t = 0.0156 + (0:49)' * 0.1;
sigma = 0.2;
nuisance = [5.15 0.14];
nTrial = 5;
fprintf('=== part 1: upper bound up to %.1f deg ===\n', cfg.maxAngleDeg);

% ---------------------------------------------------------------- part 1
upper = struct('param', {}, 'value', {}, 'estMean', {}, 'relErr', {}, 'ok', {});
q = 0;
for alphaDeg = [12 15 18 21 25]
    for k = 1:nTrial
        p = [0, 1.7, alphaDeg, 0, 0.7];
        rng(2000 + k);
        Y = isar_simulate_cone_peaks(t, P, [p, 0, 0, 0], cfg, sigma, nuisance);
        id = isar_identify_tip_mode(t, Y(:, 1), cfg);
        q = q + 1;
        upper(q).param = 'alpha';
        upper(q).value = alphaDeg;
        upper(q).estMean = id.selected.theta(2);
        upper(q).relErr = abs(id.selected.theta(2) - alphaDeg) / alphaDeg;
        upper(q).ok = double(upper(q).relErr < 0.10 && id.mode == "precession");
    end
end
for ampDeg = [10 14 18 21 25]
    for k = 1:nTrial
        p = [0, 1.7, 3, ampDeg, 0.7];
        rng(3000 + k);
        Y = isar_simulate_cone_peaks(t, P, [p, 0, 0, 0], cfg, sigma, nuisance);
        id = isar_identify_tip_mode(t, Y(:, 1), cfg);
        q = q + 1;
        upper(q).param = 'A';
        upper(q).value = ampDeg;
        upper(q).estMean = id.selected.theta(3);
        upper(q).relErr = abs(id.selected.theta(3) - ampDeg) / ampDeg;
        upper(q).ok = double(upper(q).relErr < 0.10 && id.swing);
    end
end
params = {'alpha', 'A'};
fprintf('%-6s %6s %10s %10s %8s %6s\n', 'param', 'truth', 'est mean', 'rel err', 'P', 'ok');
for ii = 1:numel(params)
    nm = params{ii};
    for v = unique([upper(strcmp({upper.param}, nm)).value])
        m = strcmp({upper.param}, nm) & [upper.value] == v;
        e = mean([upper(m).estMean]);
        rel = mean([upper(m).relErr]);
        ok = mean([upper(m).ok]);
        fprintf('%-6s %6.1f %10.3f %10.4f %8.2f %6d\n', nm, v, e, rel, ok, ok >= 0.9);
    end
end

% ---------------------------------------------------------------- part 2
fprintf('\n=== part 2: resolution by delta-separation ===\n');
base = [0, 1.7, 3.0, 3.0, 0.7];          % fs fp alpha A fw
% p is [fs fp alpha A fw ...] but the fit theta is [fp alpha A fw psiP psiW c0 c1],
% so the two index sets differ.
ladders = struct( ...
    'name', {'fp', 'alpha', 'A', 'fw'}, ...
    'idxP', {2, 3, 4, 5}, ...
    'idxT', {1, 2, 3, 4}, ...
    'centre', {1.7, 3.0, 3.0, 0.7}, ...
    'delta', {[1e-5 3e-5 1e-4 3e-4 1e-3 3e-3 1e-2], ...
              [1e-4 3e-4 1e-3 3e-3 1e-2 3e-2 1e-1], ...
              [1e-3 3e-3 1e-2 3e-2 1e-1 3e-1 1.0], ...
              [1e-5 3e-5 1e-4 3e-4 1e-3 3e-3 1e-2]});
resolution = struct('name', {}, 'delta', {}, 'pass', {});
for qi = 1:numel(ladders)
    L = ladders(qi);
    fprintf('%-6s centre %6.3f\n', L.name, L.centre);
    hit = NaN;
    for d = L.delta
        nPass = 0;
        for k = 1:nTrial
            p = base;
            p(L.idxP) = L.centre + d;
            rng(4000 + k);
            Y = isar_simulate_cone_peaks(t, P, [p, 0, 0, 0], cfg, sigma, nuisance);
            fit = isar_fit_tip_motion_mode(t, Y(:, 1), 'H3', ...
                struct('cfg', cfg, 'ranking', 'short'));
            if abs(fit.theta(L.idxT) - (L.centre + d)) < d / 2
                nPass = nPass + 1;
            end
        end
        fprintf('   delta %9.1e -> P %4.2f\n', d, nPass / nTrial);
        if nPass == nTrial && isnan(hit)
            hit = d;
        end
    end
    resolution(end + 1) = struct('name', L.name, 'delta', hit, 'pass', 1); %#ok<AGROW>
end
fprintf('\n=== resolution (smallest delta passing every trial) ===\n');
for r = resolution
    fprintf('%-6s %.1e\n', r.name, r.delta);
end
out = struct('upper', upper, 'resolution', resolution, 'sigma', sigma, 'base', base);
save(fullfile(outputDirectory, 'resolution_results.mat'), 'out');
end
