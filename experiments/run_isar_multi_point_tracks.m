function result = run_isar_multi_point_tracks(options)
%RUN_ISAR_MULTI_POINT_TRACKS Per-scatterer slow-time tracks from full-cone RF.
%   Same native ISARpath frontend and CA prior correction as
%   run_isar_tip_slowtime_tracks, but the four scatterers are extracted
%   separately: P1 (tip) plus P2..P4 (base). The tip alone cannot see the
%   spin rate, so the base points are the observable that carries fs.
%
%   Each scatterer is isolated by the library echo generator's per-scatterer
%   channels, and its range gate is chosen data-driven: the maximum of that
%   channel's own range profile plus/minus three bins. No scenario truth is
%   used for gating, and the gate is re-selected every CPI.
%
%   OPTIONS: Cases (table), OutputDirectory, NativeWorkers, GenerateFigures,
%   DurationSec, FrameIndices.
%
%   Writes <out>/multi_point_tracks.csv and <out>/multi_point_cases.csv.
%
%   See also IDENTIFY_ISAR_SPIN_FROM_TAIL_POINTS, ISAR_CONE_BODY_POINTS.

if nargin < 1
    options = struct();
end
[baseCfg, sensor, base, cases] = isar_grid_setup();
if isfield(options, 'Cases') && ~isempty(options.Cases)
    cases = options.Cases;
else
    cases = cases(cases.spinFrequencyHz == 1, :);
end
root = fileparts(fileparts(mfilename('fullpath')));
out = char(opt(options, 'OutputDirectory', ...
    fullfile(root, 'results', 'isar_multi_point_tracks_20261021')));
if ~isfolder(out)
    mkdir(out);
end
duration = opt(options, 'DurationSec', 5);
frames = opt(options, 'FrameIndices', 1:round(duration / .1));
cfg = struct('fs', 1e9, 'fs_scope', 200e6, 'T', 1e-3, 'Tr', 100e-6, ...
    'fL', 11.8e9, 'fH', 17.8e9, 'kr', 6e13, 'propagationSpeed', 3e8, ...
    'fc', 14.8e9, 'lambda', 3e8 / 14.8e9, 'Ta', .001, 'Na', 32, ...
    'startTimeSec', 0, 'nativeWorkers', opt(options, 'NativeWorkers', 8), ...
    'generatePriorWaveform', false, 'recordPolicy', "TIGHT_EMPTY_TAIL");
p = isarpath_settings(cfg);
p.rangeAlignmentMode = "FRACTIONAL";
p.phaseCompensationMode = "MODEL_CA";
nf = 64 * cfg.Na;
fd = (-nf / 2:nf / 2 - 1)' / (cfg.Ta * nf);
rx = (-512:511) * p.rangeBinSpacingM ...
    + p.cropCenterOffsetBinsEffective * p.rangeBinSpacingM;
gateHalfWidthBins = 3;
% Library channel order: 1 coherent sum of four, 2 isolated tip, 3:5 tails.
channels = [2, 3, 4, 5];

diary(fullfile(out, 'run_log.txt'));
diaryGuard = onCleanup(@() diary('off')); %#ok<NASGU>
timer = tic;
tracks = table();
for ci = 1:height(cases)
    q = cases(ci, :);
    work = baseCfg;
    work.bodyRotationMode = q.mode;
    work.spinFrequencyHz = q.spinFrequencyHz;
    work.precessionFrequencyHz = q.precessionFrequencyHz;
    work.precessionHalfAngleDeg = q.precessionHalfAngleDeg;
    work.swingFrequencyHz = q.swingFrequencyHz;
    work.swingAngleMinDeg = -q.swingAmplitudeDeg;
    work.swingAngleMaxDeg = q.swingAmplitudeDeg;
    target = isar_apply_micro_motion(base, work);
    local = table();
    for fi = 1:numel(frames)
        frame = frames(fi);
        epoch = (frame - 1) * .1;
        t = epoch + (0:cfg.Na - 1)' * cfg.Ta;
        prior = struct('epochSec', epoch, ...
            'rangeM', 10050 + 3400 * epoch + 50 * epoch^2, ...
            'velocityMps', 3400 + 100 * epoch, 'accelerationMps2', 100, ...
            'beatOffsetHz', 20e6);
        prior.loStartSec = 2 * prior.rangeM / 3e8 - 20e6 / cfg.kr;
        cfg.T = ceil((prior.loStartSec + cfg.Tr + 16e-6) * cfg.fs / 32000) ...
            * 32000 / cfg.fs;
        [H, ~, audit] = isarpath_synthesize_profiles(target, sensor, cfg, prior, t, p);
        assert(audit.maximumRawBeatHz < 100e6 && audit.coherentSumError < 1e-11);
        [~, shifts] = isarpath_align_profiles(H(:, :, 1), p);
        H = isarpath_align_profiles(H, p, shifts);
        ref = isarpath_prior_phase(cfg, prior, t, p, shifts, "MODEL_CA");
        H = H .* ref.phaseCorrection;
        peakHz = zeros(1, 4);
        gateM = zeros(1, 4);
        mainInGate = zeros(1, 4);
        for k = 1:4
            profile = sum(abs(H(:, :, channels(k))) .^ 2, 1);
            [~, centreBin] = max(profile);
            gate = centreBin + (-gateHalfWidthBins:gateHalfWidthBins);
            gate = gate(gate >= 1 & gate <= size(H, 2));
            spectrum = fftshift(fft(H(:, gate, channels(k)), nf, 1), 1) / cfg.Na;
            peakHz(k) = peak_frequency(sum(abs(spectrum) .^ 2, 2), fd);
            gateM(k) = rx(centreBin);
            mainInGate(k) = any(gate == centreBin);
        end
        receiveMid = mean(t) + 2 * prior.rangeM / 3e8 + cfg.Tr / 2;
        local = [local; table(q.caseId, frame, epoch, receiveMid, ...
            peakHz(1), peakHz(2), peakHz(3), peakHz(4), ...
            gateM(1), gateM(2), gateM(3), gateM(4), ...
            all(mainInGate), audit.maximumRawBeatHz, audit.coherentSumError, ...
            'VariableNames', {'caseId', 'frameIndex', 'epochSec', ...
            'receiveMidpointSec', 'peakP1Hz', 'peakP2Hz', 'peakP3Hz', 'peakP4Hz', ...
            'gateP1M', 'gateP2M', 'gateP3M', 'gateP4M', 'allPointsGated', ...
            'maximumRawBeatHz', 'coherentSumError'})]; %#ok<AGROW>
        fprintf('%s | %d/%d | P1 %8.3f P2 %8.3f P3 %8.3f P4 %8.3f Hz | %.1f s\n', ...
            q.caseId, fi, numel(frames), peakHz, toc(timer));
    end
    tracks = [tracks; local]; %#ok<AGROW>
end
writetable(tracks, fullfile(out, 'multi_point_tracks.csv'));
writetable(cases, fullfile(out, 'multi_point_cases.csv'));
result = struct('cases', cases, 'tracks', tracks, 'outputDirectory', string(out), ...
    'configuration', cfg, 'processing', p, 'gateHalfWidthBins', gateHalfWidthBins, ...
    'channelOrder', channels);
save(fullfile(out, 'multi_point_tracks.mat'), 'result');
fprintf('DONE %s | %.1f s\n', out, toc(timer));
end

function hz = peak_frequency(power, f)
[~, i] = max(power);
offset = 0;
if i > 1 && i < numel(power)
    y = log(max(power(i + (-1:1)), realmin));
    den = y(1) - 2 * y(2) + y(3);
    if den < 0
        offset = .5 * (y(1) - y(3)) / den;
    end
end
hz = f(i) + offset * (f(2) - f(1));
end

function value = opt(s, name, fallback)
if isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = fallback;
end
end
