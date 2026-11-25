function out = isar_chain_identify_ch1(options)
%ISAR_CHAIN_IDENTIFY_CH1 Identification entry fed by the total-field image.
%   Two corrections to the earlier entry point:
%
%     * the acceleration comes from the range measurements alone
%       (R_true + 2 cm noise, batch quadratic fit); the image is never used
%       to estimate it, so no labelled scatterer channel is involved;
%     * the slow-time Doppler track is taken from the total field H(:,:,1)
%       (the four-scatterer sum that a real receiver produces), gated to
%       +/-3 bins around that field's own range peak, with parabolic
%       interpolation of the peak. H(:,:,2) is used only to report how much
%       of the gate energy belongs to the tip, as a diagnostic.
%
%   OutputDirectory is the parent of a new per-run directory.
%   Options: DurationSec, RangeModel (centre/tip), RecordMode (measured/exact).

if nargin < 1
    options = struct();
end
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'functions')); ensure_md_zt_paths();
outParent = char(opt(options, 'OutputDirectory', ...
    fullfile(root, 'results', 'isar_chain_identify_ch1_20261123')));
Tobs = opt(options, 'DurationSec', 5);
sigR = 0.02;

cfg = struct('fs', 1e9, 'fs_scope', 200e6, 'T', 1e-3, 'Tr', 100e-6, ...
    'fL', 11.8e9, 'fH', 17.8e9, 'kr', 6e13, 'propagationSpeed', 3e8, ...
    'fc', 14.8e9, 'lambda', 3e8 / 14.8e9, 'Ta', 1e-3, 'Na', 32, ...
    'startTimeSec', 0, 'nativeWorkers', 8, ...
    'generatePriorWaveform', false, 'recordPolicy', "TIGHT_EMPTY_TAIL");
sensor = struct('position_I', [0 0 0], 'velocity_I', [0 0 0], ...
    'carrierFrequencyHz', 14.8e9, 'wavelengthM', cfg.lambda, 'dopplerSign', -1);
range0 = 10050;
vel0 = 10.2 * 340;      % radial projection onto the LOS (3468 m/s)
aTrue = -15 * 10;       % radial projection onto the LOS (-150 m/s^2)
spec = struct('Ntip', 1, 'radiusM', .3, 'phase0Rad', pi / 6, ...
    'axialSeparationM', 1.5, 'alpha0', [1; .3; .3; .3]);
tLocal = (0:cfg.Na - 1).' * cfg.Ta;
nPad = 8;
nf = cfg.Na * nPad;
fd = (-nf / 2:nf / 2 - 1)' / (cfg.Na * cfg.Ta * nPad);
scenes = struct( ...
    'tag', {'prec_fp1_a10', 'prec_fp3_a5', 'prec_swing_w3', 'prec_swing_w5'}, ...
    'mode', {"coneSpin", "coneSpin", "coneSwingSpin", "coneSwingSpin"}, ...
    'fp', {1, 3, 1, 2}, 'alpha', {10, 5, 10, 20}, ...
    'fw', {0, 0, 1, 1}, 'A', {0, 0, 3, 5});
idCfg = isar_tip_motion_cfg();
tF = (0:0.1:Tobs)';
nF = numel(tF);
% Processing-side data sources, declared once and reported with the results.
rangeModel = lower(string(opt(options, 'RangeModel', "centre")));
recMode = lower(string(opt(options, 'RecordMode', "measured")));
if ~ismember(rangeModel, ["centre", "tip"])
    error('ISAR:Chain:RangeModel', 'RangeModel must be "centre" or "tip".');
end
if ~ismember(recMode, ["measured", "exact"])
    error('ISAR:Chain:RecordMode', 'RecordMode must be "measured" or "exact".');
end
sigV = 0.3;
beatOffsetHz = opt(options, 'BeatOffsetHz', 20e6);
runInfo = isar_ca_run_info(outParent,mfilename,rangeModel,recMode,sigR,sigV);
outDir = char(runInfo.outputDirectory);
runInfo.beatOffsetHz = beatOffsetHz;
fprintf(['processing sources: range=%s(%s), velocity=%s(%s), ' ...
    'acceleration=range-fit, beatOffset=%.0f MHz\n'], ...
    runInfo.rangeModel, recMode, runInfo.velocityModel, recMode, beatOffsetHz / 1e6);
rows = struct('scene', {}, 'aHat', {}, 'mode', {}, 'fp', {}, 'alpha', {}, ...
    'A', {}, 'fw', {}, 'rmsHz', {}, 'tipToCoherentRatio', {}, 'fpTrue', {}, 'recordSeed', {});
for s = 1:numel(scenes)
    q = scenes(s);
    timer = tic;
    % ---- measurement record: the only place scene truth is consumed.
    % RangeModel "centre" (default) measures the body-centre range, which is
    % what a tracker reports; micro-motion is a relative motion of the
    % scatterers about that centre and never enters it. "tip" is an explicit,
    % separately labelled option that models a point-locked tracker, and it
    % passes ONLY the tip's 3-vector.
    switch rangeModel
        case "centre"
            rMicro = 0;
        case "tip"
            bodyPoints = isar_cone_body_points();
            rMicro = isar_cone_point_kinematics(tF, bodyPoints(1, :), ...
                [0, q.fp, q.alpha, q.A, max(q.fw, 1e-6), 0, 0, 0], ...
                idCfg);
        otherwise
            error('ISAR:Chain:RangeModel', 'Unknown RangeModel.');
    end
    Rtrue = range0 + vel0 * tF + 0.5 * aTrue * tF .^ 2 + rMicro;
    Vtrue = vel0 + aTrue * tF;
    rec = isar_cone_measurement_record(tF, Rtrue, Vtrue, ...
        struct('Mode', recMode, 'RangeSigmaM', sigR, ...
        'VelocitySigmaMps', sigV, 'Seed', 9100 + s));
    aHat = isar_ca_record_acceleration(rec);
    % ---- frames
    yTrack = zeros(nF, 1);
    tipRatio = nan(nF, 1);
    for k = 1:nF
        tk = tF(k);
        % Processing consumes the record only; no truth argument exists here.
        prior = isar_cone_processing_prior(rec, k, aHat, beatOffsetHz, cfg);
        ck = cfg;
        ck.T = ceil((prior.loStartSec + ck.Tr + 16e-6) * ck.fs / 32000) ...
            * 32000 / ck.fs;
        pk = isarpath_settings(ck);
        pk.rangeAlignmentMode = "INTEGER";
        pk.rangeMigrationMode = "FFT2048_SINGLE";
        motion = struct('type', 'accelerated', 'p0_I', [range0 0 0], ...
            'v0_I', vel0 * [1 0 0], 'a0_I', aTrue * [1 0 0], ...
            'attitude0YprRad', [-3 * pi / 4 0 0]);
        base = scatterermodel.createLegacyDiskTarget(spec, motion);
        work = ck;
        work.bodyRotationMode = q.mode;
        work.spinFrequencyHz = 1;
        work.precessionFrequencyHz = q.fp;
        work.precessionHalfAngleDeg = q.alpha;
        work.swingFrequencyHz = max(q.fw, 1e-6);
        work.swingAngleMinDeg = -q.A;
        work.swingAngleMaxDeg = q.A;
        target = isar_apply_micro_motion(base, work);
        [H, ~, ~] = isarpath_synthesize_profiles(target, sensor, ck, prior, ...
            tLocal + tk, pk);
        [~, sh] = isarpath_align_profiles(H(:, :, 1), pk);
        H = isarpath_align_profiles(H, pk, sh);
        ref = isarpath_prior_phase(ck, prior, tLocal + tk, pk, sh, "MODEL_CA");
        corr = ref.phaseCorrection;
        total = H(:, :, 1) .* corr;                 % total field only
        prof = sum(abs(total) .^ 2, 1);
        [~, b] = max(prof);
        gate = b + (-3:3);
        gate = gate(gate >= 1 & gate <= size(total, 2));
        spec_ = fftshift(fft(total(:, gate), nf, 1), 1) / cfg.Na;
        P = sum(abs(spec_) .^ 2, 2);
        [~, i0] = max(P);
        offs = 0;
        if i0 > 1 && i0 < nf
            y = log(max(P(i0 + (-1:1)), realmin));
            den = y(1) - 2 * y(2) + y(3);
            if den < 0
                offs = 0.5 * (y(1) - y(3)) / den;
            end
        end
        yTrack(k) = fd(i0) + offs * (fd(2) - fd(1));
        % Diagnostic only: tip-channel energy over coherent total-field
        % energy. It never feeds the gate, the track or the mode decision.
        tipRatio(k) = isar_tip_to_coherent_ratio(H, gate);
    end
    id = isar_identify_tip_mode(tF, yTrack, idCfg);
    th = id.selected.theta;
    rows(end + 1) = struct('scene', q.tag, 'aHat', aHat, 'mode', id.mode, ...
        'fp', th(1), 'alpha', th(2), 'A', th(3), 'fw', th(4), ...
        'rmsHz', id.selected.rmsHz, ...
        'tipToCoherentRatio', median(tipRatio, 'omitnan'), ...
        'fpTrue', q.fp, 'recordSeed', rec.seed); %#ok<AGROW>
    fprintf(['%-16s aHat %7.2f  %-16s fp %5.2f (true %4.1f)  alpha %5.2f ' ...
        '(true %4.1f)  A %5.2f  fw %5.2f  rms %5.2f  tip/coherent %.3f  [%4.0f s]\n'], ...
        q.tag, aHat, id.mode, th(1), q.fp, th(2), q.alpha, th(3), th(4), ...
        id.selected.rmsHz, median(tipRatio, 'omitnan'), toc(timer));
end
tbl = isar_ca_tag_summary(struct2table(rows),runInfo);
writetable(tbl, fullfile(outDir, 'chain_identify_ch1.csv'));
save(fullfile(outDir, 'chain_identify_ch1.mat'), 'tbl', 'scenes', 'Tobs','runInfo');
fprintf('DONE %s\n', outDir);
out = struct('table', tbl,'outputDirectory',runInfo.outputDirectory,'runInfo',runInfo);
end

function v = opt(s, name, fallback)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = fallback;
end
end
