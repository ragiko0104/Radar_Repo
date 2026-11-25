function out = sweep_isar_bandwidth(outputDirectory)
%SWEEP_ISAR_BANDWIDTH End-to-end effect of range resolution on the frozen chain.
%   The estimator itself has no bandwidth parameter: bandwidth enters through
%   the range gate of step 2. Four bandwidths are tested, 6 GHz being the
%   reference (c/2B = 2.5 cm). For each bandwidth the four-scatterer frame is
%   synthesised with a Hann range point spread function at that resolution,
%   the track is extracted with the production gate, and the frozen
%   identification is run on it.
%
%   Noise scales with bandwidth: a matched filter keeps the peak and passes
%   noise proportional to B, so the per-cell noise is sigma_ref*sqrt(B/B_ref).

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'functions'));
if nargin < 1
    outputDirectory = fullfile(root, 'results', 'isar_bandwidth_20261027');
end
if ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end
cfg = isar_tip_motion_cfg();
P = isar_cone_body_points();
tTrack = 0.0156 + (0:49)' * 0.1;
Ta = 1e-3;
Na = 32;
nf = 2048;
nFrame = 50;
gain = [1, .3, .3, .3];
nuisance = [5.15 0.14];
scenes = struct('tag', {'precession', 'precession_swing', 'spin3_precession', ...
    'spin3_precession_swing'}, ...
    'p', {[0, 1.7, 3, 0, 0.7], [0, 1.7, 3, 3, 0.7], ...
    [3, 3, 3, 0, 0.7], [3, 3, 3, 3, 0.7]});
bandwidths = [0.1, 1, 3, 6] * 1e9;
sd = 0.2;
rows = struct('scene', {}, 'bandwidthGHz', {}, 'rangeResM', {}, ...
    'gateWidthM', {}, 'pointsInGate', {}, 'mode', {}, 'fp', {}, 'alpha', {}, ...
    'A', {}, 'fw', {}, 'rmsHz', {});
for si = 1:numel(scenes)
    sc = scenes(si);
    p = [sc.p, 0, 0, 0];
    for B = bandwidths
        dr = 3e8 / (2 * B);
        sigma = 0.03 * sqrt(0.025 / dr);
        % range grid: cone plus the gate width on each side
        halfSpan = 1.1 + 4 * dr;
        rGrid = (-halfSpan:dr:halfSpan)';
        nR = numel(rGrid);
        fd = (-nf / 2:nf / 2 - 1)' / (Ta * nf);
        wAz = .5 - .5 * cos(2 * pi * (0:Na - 1)' / (Na - 1));
        yTrack = zeros(nFrame, 1);
        nInGate = zeros(nFrame, 1);
        gateM = 0;
        for fi = 1:nFrame
            tk = (fi - 1) * 0.1 + (0:Na - 1)' * Ta;
            tMid = (fi - 1) * 0.1 + (Na - 1) / 2 * Ta;
            H = zeros(Na, nR);
            rMid = zeros(1, 4);
            for k = 1:4
                [R, ~] = isar_cone_point_kinematics(tk, P(k, :), p, cfg);
                ph = exp(1j * 4 * pi * R / cfg.lambdaM);
                [Rm, ~] = isar_cone_point_kinematics(tMid, P(k, :), p, cfg);
                rMid(k) = Rm;
                bin = (Rm - rGrid(1)) / dr;
                x = repmat(bin, 1, nR) - repmat((0:nR - 1), Na, 1);
                psf = .5 * sincLocal(x) + .25 * sincLocal(x - 1) + .25 * sincLocal(x + 1);
                H = H + gain(k) * (ph * ones(1, nR)) .* psf;
            end
            rng(7000 + 100 * si + fi);
            H = H + sigma * (randn(Na, nR) + 1j * randn(Na, nR)) / sqrt(2);
            prof = sum(abs(H) .^ 2, 1).';
            [~, cBin] = max(prof);
            gate = cBin + (-3:3);
            gate = gate(gate >= 1 & gate <= nR);
            if fi == 1
                gateM = rGrid(gate(end)) - rGrid(gate(1)) + dr;
                lo = rGrid(gate(1)) - dr / 2;
                hi = rGrid(gate(end)) + dr / 2;
                nInGate(1) = sum(rMid >= lo & rMid <= hi);
            end
            Xg = fftshift(fft(H(:, gate) .* wAz, nf, 1), 1) / Na;
            Pg = sum(abs(Xg) .^ 2, 2);
            [~, ig] = max(Pg);
            off = 0;
            if ig > 1 && ig < nf
                y = log(max(Pg(ig + (-1:1)), realmin));
                den = y(1) - 2 * y(2) + y(3);
                if den < 0
                    off = .5 * (y(1) - y(3)) / den;
                end
            end
            yTrack(fi) = fd(ig) + off * (fd(2) - fd(1));
        end
        rng(8000 + si);
        yMeas = yTrack + sd * randn(nFrame, 1);
        id = isar_identify_tip_mode(tTrack, yMeas, cfg);
        th = id.selected.theta;
        rows(end + 1) = struct('scene', sc.tag, 'bandwidthGHz', B / 1e9, ...
            'rangeResM', dr, 'gateWidthM', gateM, 'pointsInGate', nInGate(1), ...
            'mode', id.mode, 'fp', th(1), 'alpha', th(2), 'A', th(3), ...
            'fw', th(4), 'rmsHz', id.selected.rmsHz); %#ok<AGROW>
        fprintf('%-16s B %5.1f GHz  dr %6.3f m  gate %6.2f m  inGate %d | %-16s fp %.3f a %6.2f A %5.2f fw %.3f rms %.3f\n', ...
            sc.tag, B / 1e9, dr, gateM, nInGate(1), id.mode, th(1), th(2), th(3), th(4), ...
            id.selected.rmsHz);
    end
end
tbl = struct2table(rows);
writetable(tbl, fullfile(outputDirectory, 'bandwidth_sweep.csv'));
out = struct('table', tbl, 'scenes', scenes, 'bandwidths', bandwidths);
save(fullfile(outputDirectory, 'bandwidth_sweep.mat'), 'out');
fprintf('DONE %s\n', outputDirectory);
end

function y = sincLocal(x)
y = ones(size(x));
nz = x ~= 0;
y(nz) = sin(pi * x(nz)) ./ (pi * x(nz));
end
