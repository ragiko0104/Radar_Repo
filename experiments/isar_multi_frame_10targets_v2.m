function out = isar_multi_frame_10targets_v2(options)
%ISAR_MULTI_FRAME_10TARGETS_V2 Ten-target imaging with a range-fitted acceleration.
%   OutputDirectory is a parent directory; every invocation gets its own run.
%   Radial (line-of-sight) kinematics: the golden configuration stores the
%   target speed as 10.2 Mach and the acceleration as -15..-10.1 g along a
%   body velocity vector 45 deg off the LOS, so the quantities that enter the
%   echo and the processing are the radial projections v_r = speed/sqrt(2) and
%   a_r = a/sqrt(2) (cfg.radial_velocity_mps, cfg.radial_acceleration_mps2).
%   This script uses those radial values, not the full vector magnitudes.
%
%   Acceleration comes from range measurements only: R(t) + 2 cm noise, a
%   batch fit of r0 + v0 t + 0.5 a t^2. No image peak and no labelled
%   scatterer channel take part in the estimate. The images are formed from
%   the total field H(:,:,1) and the tip mainlobe is measured on that field.

if nargin < 1
    options = struct();
end
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'functions')); ensure_md_zt_paths();
outParent = char(opt(options, 'OutputDirectory', ...
    fullfile(root, 'results', 'isar_multi_frame_10targets_v2_20261123')));
Tobs = opt(options, 'DurationSec', 2);
recMode = lower(string(opt(options, 'RecordMode', 'measured')));
sigR = 0.02;
sigV = 0.3;

cfg = struct('fs', 1e9, 'fs_scope', 200e6, 'T', 1e-3, 'Tr', 100e-6, ...
    'fL', 11.8e9, 'fH', 17.8e9, 'kr', 6e13, 'propagationSpeed', 3e8, ...
    'fc', 14.8e9, 'lambda', 3e8 / 14.8e9, 'Ta', 1e-3, 'Na', 32, ...
    'startTimeSec', 0, 'nativeWorkers', 8, ...
    'generatePriorWaveform', false, 'recordPolicy', "TIGHT_EMPTY_TAIL");
sensor = struct('position_I', [0 0 0], 'velocity_I', [0 0 0], ...
    'carrierFrequencyHz', 14.8e9, 'wavelengthM', cfg.lambda, 'dopplerSign', -1);
% Per-target parameters come from the approved golden config record, not from
% a copied constant table: the local-oscillator offset alternates 20/45 MHz,
% and beatOffsetHz and loStartSec must both use the SAME recorded value.
gcfg = isar_golden_target_config(opt(options, 'GoldenConfigPath', []));
runInfo = isar_ca_run_info(outParent,mfilename,"centre",recMode,sigR,sigV);
runInfo.goldenConfigPath = string(gcfg.sourcePath);
outDir = char(runInfo.outputDirectory);
ranges = gcfg.rangeM;
biasHz = gcfg.beatOffsetHz;
velR = gcfg.velocityMps;
accR = gcfg.accelerationMps2;
spec = struct('Ntip', 1, 'radiusM', .3, 'phase0Rad', pi / 6, ...
    'axialSeparationM', 1.5, 'alpha0', [1; .3; .3; .3]);
tLocal = (0:cfg.Na - 1).' * cfg.Ta;
nPad = 8;
fd = (-cfg.Na * nPad / 2:cfg.Na * nPad / 2 - 1)' / (cfg.Na * nPad * cfg.Ta);
tF = (0:0.1:Tobs)';
nF = numel(tF);
nTarget = numel(ranges);
img = cell(nTarget, 1);
stats = struct('target', {}, 'rangeM', {}, 'vRadial', {}, 'aTrue', {}, ...
    'aHat', {}, 'aErr', {}, 'fwhmCa', {}, 'fwhmCv', {}, 'fwhmProm', {}, ...
    'peakHz', {}, 'goldenBiasHz', {}, 'usedBiasHz', {}, 'loStartSec', {}, ...
    'auditFrameIndex',{},'auditEpochSec',{},'auditRangeM',{},'recordSeed',{});
loRows = cell(nF,nTarget); % MATLAB one-based frame index, Golden zero-based target ID
fprintf(['ten targets, radial v %.1f m/s, radial a %.2f..%.2f m/s^2, ' ...
    '%d frames over %.0f s\n'], velR(1), min(accR), max(accR), nF, Tobs);
fprintf('%6s %8s %10s %10s %10s %12s\n', 'id0', 'range', 'goldenBias', ...
    'usedBias', 'loStart', 'a_r');
for it = 1:nTarget
    range0 = ranges(it);
    aTrue = accR(it);
    bias = biasHz(it);
    timer = tic;
    % ---- acceleration from range measurements only
    Rtrue = range0 + velR(it) * tF + 0.5 * aTrue * tF .^ 2;
    rec = isar_cone_measurement_record(tF, Rtrue, ...
        velR(it) + aTrue * tF, struct('Mode', recMode, ...
        'RangeSigmaM', sigR, 'VelocitySigmaMps', sigV, 'Seed', 7300 + it));
    aHat = isar_ca_record_acceleration(rec);
    fwhm = nan(nF, 3);
    peakHz = nan(nF, 3);
    for k = 1:nF
        tk = tF(k);
        prior = isar_cone_processing_prior(rec, k, aHat, bias, cfg);
        prior.accelerationMps2 = 0;    % echo/dechirp prior: CV for the window
        loRows{k,it} = isar_ca_lo_audit(gcfg.targetId(it),k,biasHz(it),prior);
        ck = cfg;
        ck.T = ceil((prior.loStartSec + ck.Tr + 16e-6) * ck.fs / 32000) ...
            * 32000 / ck.fs;
        pk = isarpath_settings(ck);
        pk.rangeAlignmentMode = "INTEGER";
        pk.rangeMigrationMode = "FFT2048_SINGLE";
        motion = struct('type', 'accelerated', 'p0_I', [range0 0 0], ...
            'v0_I', velR(it) * [1 0 0], 'a0_I', aTrue * [1 0 0], ...
            'attitude0YprRad', [-3 * pi / 4 0 0]);
        base = scatterermodel.createLegacyDiskTarget(spec, motion);
        work = ck;
        work.bodyRotationMode = "coneSpin";
        work.spinFrequencyHz = 1;
        work.precessionFrequencyHz = 1;
        work.precessionHalfAngleDeg = 20;
        work.swingFrequencyHz = 1;
        work.swingAngleMinDeg = 0;
        work.swingAngleMaxDeg = 0;
        target = isar_apply_micro_motion(base, work);
        [H, ~, ~] = isarpath_synthesize_profiles(target, sensor, ck, prior, ...
            tLocal + tk, pk);
        [~, sh] = isarpath_align_profiles(H(:, :, 1), pk);
        H = isarpath_align_profiles(H, pk, sh);
        schemes = ["VARIANCE_UNIT", "MODEL_CV", "MODEL_CA"];
        for si = 1:numel(schemes)
            pw = pk;
            pw.phaseCompensationMode = schemes(si);
            acc = 0;
            if schemes(si) == "MODEL_CA"
                acc = aHat;
            end
            pr = prior;
            pr.accelerationMps2 = acc;
            ref = isarpath_prior_phase(ck, pr, tLocal + tk, pw, sh, ...
                schemes(si));
            [~, ~, corr] = isarpath_phase_image(H, pw, ref);
            total = H(:, :, 1) .* corr;            % total field only
            [fwhm(k, si), peakHz(k, si)] = tip_width( ...
                radar.ISAR_imaging(total, nPad, char(pw.windowType)), fd);
            if si == 3 && k == round(nF / 2)
                img{it} = radar.ISAR_imaging(total, nPad, char(pw.windowType));
            end
        end
    end
    a = loRows{1,it}; % summary audit fields all describe frame 1
    stats(end + 1) = struct('target', a.targetId, 'rangeM', range0, ...
        'vRadial', velR(it), 'aTrue', aTrue, 'aHat', aHat, 'aErr', aHat - aTrue, ...
        'goldenBiasHz', a.goldenBiasHz, 'usedBiasHz', a.usedBiasHz, ...
        'loStartSec', a.loStartSec,'auditFrameIndex',a.frameIndex, ...
        'auditEpochSec',a.epochSec,'auditRangeM',a.rangeM,'recordSeed',rec.seed, ...
        'fwhmCa', median(fwhm(:, 3), 'omitnan'), ...
        'fwhmCv', median(fwhm(:, 2), 'omitnan'), ...
        'fwhmProm', median(fwhm(:, 1), 'omitnan'), ...
        'peakHz', median(peakHz(:, 3), 'omitnan')); %#ok<AGROW>
    fprintf(['target %2d  r=%6.0f m  bias %2.0f MHz (golden %2.0f)  a_r=%7.2f  a_hat=%7.2f  err=%5.2f  ' ...
        'FWHM prom/cv/ca = %6.1f / %7.1f / %6.1f Hz   (%.0f s)\n'], ...
        a.targetId, range0, a.usedBiasHz / 1e6, a.goldenBiasHz / 1e6, aTrue, aHat, aHat - aTrue, stats(end).fwhmProm, ...
        stats(end).fwhmCv, stats(end).fwhmCa, toc(timer));
end
tbl = isar_ca_tag_summary(struct2table(stats),runInfo);
loAudit = struct2table(vertcat(loRows{:}));
writetable(loAudit,fullfile(outDir,'lo_audit_per_frame.csv'));
writetable(tbl, fullfile(outDir, 'multi_frame_10targets_v2.csv'));
save(fullfile(outDir, 'multi_frame_10targets_v2.mat'), 'img', 'tbl', 'fd', ...
    'ranges', 'accR', 'velR', 'Tobs','runInfo','loAudit','-v7.3');
out = struct('table', tbl, 'outputDirectory', string(outDir),'runInfo',runInfo);
fprintf('DONE %s\n', outDir);
end

% ----------------------------------------------------------------------
function [fwhm, peakHz] = tip_width(img, fd)
amp = abs(img);
[~, lin] = max(amp(:));
[row, col] = ind2sub(size(amp), lin);
peak = amp(row, col);
n = size(amp, 1);
m0 = floor(n / 2) + 1;
line = circshift(amp(:, col), m0 - row);
off = ((1:n).' - m0) * abs(fd(2) - fd(1));
fwhm = NaN;
half = peak / sqrt(2);
below = line < half;
up = find(below(1:m0 - 1), 1, 'last');
dn = find(below(m0 + 1:end), 1, 'first') + m0;
if ~isempty(up) && ~isempty(dn) && up + 1 <= m0 && dn - 1 >= m0
    fLo = off(up) + (half - line(up)) * (off(up + 1) - off(up)) / ...
        (line(up + 1) - line(up));
    fHi = off(dn - 1) + (half - line(dn - 1)) * (off(dn) - off(dn - 1)) / ...
        (line(dn) - line(dn - 1));
    fwhm = fHi - fLo;
end
peakHz = fd(row);
end

function v = opt(s, name, fallback)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = fallback;
end
end
