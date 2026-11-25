function f = isar_predict_cone_peaks(tTrack, bodyPoints, p, cfg, nuisance)
%ISAR_PREDICT_CONE_PEAKS Predicted per-CPI peak Doppler of every scatterer.
%   F = ISAR_PREDICT_CONE_PEAKS(TTRACK,BODYPOINTS,P,CFG,NUISANCE) applies the
%   same observation operator as the measurement chain: for each CPI the
%   scatterer phase is sampled over the 32 pulses of the dwell, the
%   zero-padded slow-time spectrum is formed, and its peak is located with
%   log-parabolic interpolation. The intra-CPI spectrum is periodic in the
%   PRF, so a base scatterer beyond +/-prfHz/2 folds exactly as in hardware.
%
%   P is [fs fp alphaDeg ampDeg fw psiP psiW psiS]; NUISANCE is [c0 c1], the
%   constant Doppler offset and drift left by the CA prior.
%
%   Using the real operator (rather than the mean instantaneous Doppler
%   alone) matters for the base scatterers: at 280 Hz of Doppler swing the
%   instantaneous value changes by ~56 Hz inside one 32 ms dwell, and the
%   spectrum peak then differs from the dwell mean by about 1.3 Hz.

if nargin < 5 || isempty(nuisance)
    nuisance = [0, 0];
end
tTrack = double(tTrack(:));
na = round(cfg.dwellSec / 1e-3);
ta = cfg.dwellSec / na;
nf = cfg.peakFftPoints;
fd = (-nf / 2:nf / 2 - 1)' / (ta * nf);
tMat = tTrack - cfg.dwellSec / 2 + (0:na - 1) * ta;
nPoint = size(bodyPoints, 1);
if isvector(nuisance)
    nuisance = repmat(reshape(double(nuisance), 1, []), nPoint, 1);
end
nCpi = numel(tTrack);
tRef = mean(tTrack);
f = zeros(nCpi, nPoint);
for k = 1:nPoint
    [rangeM, ~] = isar_cone_point_kinematics(tMat(:), bodyPoints(k, :), p, cfg);
    phase = 4 * pi / cfg.lambdaM * reshape(rangeM, nCpi, na);
    spectrum = abs(fftshift(fft(exp(1j * phase), nf, 2), 2)) .^ 2;
    for m = 1:nCpi
        hz = peakFrequency(spectrum(m, :), fd);
        y = hz + nuisance(k, 1) + nuisance(k, 2) * (tTrack(m) - tRef);
        f(m, k) = isar_fold_doppler(y, cfg.prfHz);
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
