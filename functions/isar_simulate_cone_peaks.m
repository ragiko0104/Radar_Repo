function [fPeak, fInstant] = isar_simulate_cone_peaks(tTrack, bodyPoints, p, cfg, noiseHz, nuisance)
%ISAR_SIMULATE_CONE_PEAKS Simulate the per-CPI peak Doppler of each scatterer.
%   [FPEAK,FINSTANT] = ISAR_SIMULATE_CONE_PEAKS(TTRACK,BODYPOINTS,P,CFG,
%   NOISEHZ,NUISANCE) reproduces what run_isar_tip_slowtime_tracks measures:
%   for every CPI it samples the scatterer phase over the 32 pulses of the
%   dwell, takes the zero-padded slow-time spectrum, and picks the peak with
%   the same log-parabolic interpolation the experiment uses.
%
%   This is the faithful forward model: the intra-CPI spectrum is periodic in
%   the PRF, so when a base scatterer exceeds +/-prfHz/2 its peak folds
%   exactly as it would in hardware. FPEAK is that folded peak; FINSTANT is
%   the dwell-averaged instantaneous Doppler for comparison.
%
%   NUISANCE is [c0 c1]: a constant Doppler offset and a drift, standing in
%   for the residual RVP/range-coupling term the CA prior leaves behind
%   (about +5.15 Hz and +0.14 Hz/s in the released tip data).

if nargin < 5 || isempty(noiseHz)
    noiseHz = 0.05;
end
if nargin < 6 || isempty(nuisance)
    nuisance = [5.15, 0.14];
end
tTrack = double(tTrack(:));
nPoint = size(bodyPoints, 1);
na = round(cfg.dwellSec / 1e-3);
ta = cfg.dwellSec / na;
nf = cfg.peakFftPoints;
fd = (-nf / 2:nf / 2 - 1)' / (ta * nf);
tRef = mean(tTrack);
fPeak = zeros(numel(tTrack), nPoint);
fInstant = isar_multi_point_doppler(tTrack, bodyPoints, p, cfg, true);
for k = 1:nPoint
    for m = 1:numel(tTrack)
        pulses = tTrack(m) - cfg.dwellSec / 2 + (0:na - 1)' * ta;
        [rangeM, ~] = isar_cone_point_kinematics(pulses, bodyPoints(k, :), p, cfg);
        phase = 4 * pi / cfg.lambdaM * rangeM;
        spectrum = abs(fft(exp(1j * phase), nf)) .^ 2;
        spectrum = fftshift(spectrum);
        hz = peakFrequency(spectrum, fd);
        y = hz + nuisance(1) + nuisance(2) * (tTrack(m) - tRef) + noiseHz * randn();
        fPeak(m, k) = isar_fold_doppler(y, cfg.prfHz);
    end
end
end

function hz = peakFrequency(power, f)
[~, i] = max(power);
offset = 0;
if i > 1 && i < numel(power)
    y = log(max(power(i + (-1:1)), realmin));
    den = y(1) - 2 * y(2) + y(3);
    if den < 0
        offset = 0.5 * (y(1) - y(3)) / den;
    end
end
hz = f(i) + offset * (f(2) - f(1));
end
