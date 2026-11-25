function out = measure_tip_2d_resolution(options)
%MEASURE_TIP_2D_RESOLUTION Range and Doppler width of the tip point response.
%   Three scenes at the golden waveform parameters, all processed through the
%   fitted-acceleration chain (record -> fit a -> MODEL_CA -> azimuth FFT):
%
%     ideal      no micro-motion (point target control)
%     precession spin 1 Hz + precession 1 Hz / 10 deg
%     prec+nut   spin 1 Hz + precession 1 Hz / 10 deg + nutation 1 Hz / 5 deg
%
%   Measure the two -3 dB cuts on the isolated tip channel of the middle
%   frame. The ideal control has one reflecting point, fixed in the body,
%   with translation only. All scenes use the same tracker record.
%   Set IdealOnly=true to run the control on its own.

if nargin < 1
    options = struct();
end
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'functions')); ensure_md_zt_paths();
defaultDir = fullfile(root, 'results', 'tip_2d_resolution_20261124');
idealOnly = opt(options, 'IdealOnly', false);
if idealOnly
    defaultDir = fullfile(defaultDir, 'ideal_only');
end
outDir = char(opt(options, 'OutputDirectory', defaultDir));
if ~isfolder(outDir)
    mkdir(outDir);
end
fontName = 'Microsoft YaHei';
Tobs = 2;
cfg = struct('fs', 1e9, 'fs_scope', 200e6, 'T', 1e-3, 'Tr', 100e-6, ...
    'fL', 11.8e9, 'fH', 17.8e9, 'kr', 6e13, 'propagationSpeed', 3e8, ...
    'fc', 14.8e9, 'lambda', 3e8 / 14.8e9, 'Ta', 1e-3, 'Na', 32, ...
    'startTimeSec', 0, 'nativeWorkers', 8, ...
    'generatePriorWaveform', false, 'recordPolicy', "TIGHT_EMPTY_TAIL");
sensor = struct('position_I', [0 0 0], 'velocity_I', [0 0 0], ...
    'carrierFrequencyHz', 14.8e9, 'wavelengthM', cfg.lambda, 'dopplerSign', -1);
range0 = 10050;
vel0 = 3468;
aTrue = -150;
bias = 20e6;
spec = struct('Ntip', 1, 'radiusM', .3, 'phase0Rad', pi / 6, ...
    'axialSeparationM', 1.5, 'alpha0', [1; .3; .3; .3]);
tLocal = (0:cfg.Na - 1).' * cfg.Ta;
nPad = 8;
nf = cfg.Na * nPad;
fd = (-nf / 2:nf / 2 - 1)' / (cfg.Na * cfg.Ta * nPad);
tF = (0:0.1:Tobs)';
nF = numel(tF);
kMid = round(nF / 2);
scenes = struct( ...
    'tag', {'ideal', 'precession', 'prec_nutation'}, ...
    'mode', {"spin", "coneSpin", "coneSwingSpin"}, ...
    'fs', {0, 1, 1}, 'fp', {0, 1, 1}, 'alpha', {0, 10, 10}, ...
    'A', {0, 0, 5}, 'fw', {0, 0, 1});
if idealOnly
    scenes = scenes(1);
end
% Nominal range resolution and rectangular-aperture Doppler FWHM.
rangeLimitM = cfg.propagationSpeed / (2 * (cfg.fH - cfg.fL));
dopLimitHz = 0.885893 / (cfg.Na * cfg.Ta);
fprintf('waveform limits: range c/2B = %.3f cm, Doppler 0.886/T = %.2f Hz\n', ...
    rangeLimitM * 100, dopLimitHz);
rec = isar_cone_measurement_record(tF, ...
    range0 + vel0 * tF + 0.5 * aTrue * tF .^ 2, vel0 + aTrue * tF, ...
    struct('Mode', "measured", 'RangeSigmaM', 0.02, ...
    'VelocitySigmaMps', 0.3, 'Seed', 22));
aHat = isar_ca_record_acceleration(rec);
rows = struct('scene', {}, 'rangeFwhmCm', {}, 'dopFwhmHz', {}, ...
    'rangeRatio', {}, 'dopRatio', {}, 'peakFd', {}, 'drCm', {}, ...
    'recordSeed', {}, 'spinHz', {}, 'precessionHz', {}, 'halfAngleDeg', {}, ...
    'swingHz', {}, 'swingDeg', {});
for s = 1:numel(scenes)
    q = scenes(s);
    sceneSpec = spec;
    if strcmp(q.tag, 'ideal')
        sceneSpec.alpha0 = [1; 0; 0; 0];
    end
    img = [];
    for k = 1:nF
        tk = tF(k);
        prior = isar_cone_processing_prior(rec, k, aHat, bias, cfg);
        % The echo/dechirp window uses the range+velocity prior; the phase
        % compensation below keeps aHat. Zeroing the struct in place would
        % silently downgrade MODEL_CA to a CV correction.
        echoPrior = prior;
        echoPrior.accelerationMps2 = 0;
        ck = cfg;
        ck.T = ceil((prior.loStartSec + ck.Tr + 16e-6) * ck.fs / 32000) ...
            * 32000 / ck.fs;
        pk = isarpath_settings(ck);
        pk.rangeAlignmentMode = "INTEGER";
        pk.rangeMigrationMode = "FFT2048_SINGLE";
        motion = struct('type', 'accelerated', 'p0_I', [range0 0 0], ...
            'v0_I', vel0 * [1 0 0], 'a0_I', aTrue * [1 0 0], ...
            'attitude0YprRad', [-3 * pi / 4 0 0]);
        base = scatterermodel.createLegacyDiskTarget(sceneSpec, motion);
        work = ck;
        work.bodyRotationMode = q.mode;
        work.spinFrequencyHz = q.fs;
        work.precessionFrequencyHz = q.fp;
        work.precessionHalfAngleDeg = q.alpha;
        work.swingFrequencyHz = q.fw;
        work.swingAngleMinDeg = -q.A;
        work.swingAngleMaxDeg = q.A;
        [target, applied] = isar_apply_micro_motion(base, work);
        if strcmp(q.tag, 'ideal')
            assert(all([applied.spinFrequencyHz, applied.precessionFrequencyHz, ...
                applied.precessionHalfAngleDeg, applied.swingFrequencyHz, ...
                applied.swingAngleMinDeg, applied.swingAngleMaxDeg] == 0), ...
                'The ideal control must have translation only.');
        end
        [H, ~, ~] = isarpath_synthesize_profiles(target, sensor, ck, echoPrior, ...
            tLocal + tk, pk);
        if strcmp(q.tag, 'ideal')
            point = H(:, :, 2);
            assert(any(point(:) ~= 0) && ...
                norm(H(:, :, 1) - point, 'fro') <= 1e-10 * norm(point, 'fro'), ...
                'The ideal control must have exactly one reflecting channel.');
        end
        [~, sh] = isarpath_align_profiles(H(:, :, 1), pk);
        H = isarpath_align_profiles(H, pk, sh);
        ref = isarpath_prior_phase(ck, prior, tLocal + tk, pk, sh, "MODEL_CA");
        % Measure the tip's OWN point response. The isolated tip channel is a
        % labelled diagnostic tap, used here only for resolution metrology;
        % nothing in the processing or identification path reads it.
        tipCube = H(:, :, 2) .* ref.phaseCorrection;
        img = radar.ISAR_imaging(tipCube, nPad, char(pk.windowType));
        if k == kMid
            midImg = img;
            dr = pk.rangeBinSpacingM;
            peakFd = fd;
        end
    end
    % ---- 2-D peak and cuts on the mid-aperture image
    A = abs(midImg);
    [~, lin] = max(A(:));
    [r0, c0] = ind2sub(size(A), lin);
    rangeCut = A(r0, :);            % through the peak Doppler
    dopCut = A(:, c0);              % through the peak range bin
    rangeFwhm = fwhm_axis(rangeCut, c0, dr);
    dopFwhm = fwhm_axis(dopCut, r0, abs(fd(2) - fd(1)));
    rows(end + 1) = struct('scene', q.tag, 'rangeFwhmCm', rangeFwhm * 100, ...
        'dopFwhmHz', dopFwhm, 'rangeRatio', rangeFwhm / rangeLimitM, ...
        'dopRatio', dopFwhm / dopLimitHz, 'peakFd', fd(r0), 'drCm', dr * 100, ...
        'recordSeed', rec.seed, 'spinHz', applied.spinFrequencyHz, ...
        'precessionHz', applied.precessionFrequencyHz, ...
        'halfAngleDeg', applied.precessionHalfAngleDeg, ...
        'swingHz', applied.swingFrequencyHz, 'swingDeg', q.A); %#ok<AGROW>
    fprintf(['%-14s range FWHM %5.2f cm (%.2f x limit)   Doppler FWHM %6.2f Hz ' ...
        '(%.2f x limit)   peak %.1f Hz\n'], q.tag, rangeFwhm * 100, ...
        rangeFwhm / rangeLimitM, dopFwhm, dopFwhm / dopLimitHz, fd(r0));
    if s == 3
        % ---- figure for the precession + nutation case
        fig = figure('Visible', 'off', 'Color', 'w', 'Position', [60 60 1500 620]);
        subplot(1, 2, 1);
        plot((1:numel(rangeCut)) - c0, 20 * log10(rangeCut / max(rangeCut) + eps), ...
            '-', 'LineWidth', 2.0, 'Color', [0.10 0.36 0.65]);
        grid on; xlim([-15 15]); ylim([-40 2]);
        xlabel('距离单元偏移', 'FontName', fontName, 'FontSize', 18);
        ylabel('归一化幅度 / dB', 'FontName', fontName, 'FontSize', 18);
        title(sprintf('距离向剖面（−3 dB 宽 %.2f cm）', rangeFwhm * 100), ...
            'FontName', fontName, 'FontSize', 18);
        set(gca, 'FontName', fontName, 'FontSize', 16);
        subplot(1, 2, 2);
        plot(fd, 20 * log10(dopCut / max(dopCut) + eps), '-', 'LineWidth', 2.0, ...
            'Color', [0.10 0.52 0.42]);
        grid on; xlim([-200 200]); ylim([-40 2]);
        xlabel('多普勒 / Hz', 'FontName', fontName, 'FontSize', 18);
        ylabel('归一化幅度 / dB', 'FontName', fontName, 'FontSize', 18);
        title(sprintf('方位向剖面（−3 dB 宽 %.1f Hz）', dopFwhm), ...
            'FontName', fontName, 'FontSize', 18);
        set(gca, 'FontName', fontName, 'FontSize', 16);
        exportgraphics(fig, fullfile(outDir, 'tip_2d_cuts_prec_nutation.png'), ...
            'Resolution', 180);
        close(fig);
    end
end
tbl = struct2table(rows);
writetable(tbl, fullfile(outDir, 'tip_2d_resolution.csv'));
save(fullfile(outDir, 'tip_2d_resolution.mat'), 'tbl', 'rangeLimitM', ...
    'dopLimitHz');
out = struct('table', tbl, 'rangeLimitM', rangeLimitM, ...
    'dopLimitHz', dopLimitHz, 'outputDirectory', string(outDir));
fprintf('DONE %s\n', outDir);
end

% ----------------------------------------------------------------------
function w = fwhm_axis(cut, idx, step)
% -3 dB width in physical units around the peak, both sides interpolated.
cut = cut(:);
peak = max(cut);
half = peak / sqrt(2);
w = NaN;
lo = find(cut(1:idx) < half, 1, 'last');
hi = find(cut(idx:end) < half, 1, 'first') + idx - 1;
if isempty(lo) || isempty(hi) || lo + 1 > idx || hi - 1 < idx
    return;
end
xLo = lo + (half - cut(lo)) / (cut(lo + 1) - cut(lo));
xHi = hi - 1 + (cut(hi - 1) - half) / (cut(hi - 1) - cut(hi));
w = (xHi - xLo) * step;
end

function v = opt(s, name, fallback)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = fallback;
end
end
