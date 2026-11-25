function cfg = isar_tip_motion_cfg(overrides)
%ISAR_TIP_MOTION_CFG Radar and geometry constants for tip mode identification.
%   Centralises the numbers shared by the kinematics, the measurement
%   operator, the identification experiment and its test, so the observable
%   is defined in exactly one place.
%
%   tipLeverM          axial tip-to-centre distance of the 1.5 m cone
%   lambdaM            wavelength at the 14.8 GHz carrier
%   dwellSec           coherent processing interval (32 pulses at 1 kHz)
%   dwellSubSamples    sub-samples used to average over the dwell
%   prfHz              pulse repetition frequency of the coherent dwell; the
%                      intra-CPI Doppler is unambiguous only within
%                      [-prfHz/2, prfHz/2). Base scatterers can exceed it,
%                      which is exactly the aliasing question studied by
%                      identify_isar_spin_from_tail_points
%   sampleStepSec      slow-time sampling of the extracted track
%   minFrequencyHz     lower search bound for fp and fsw
%   maxFrequencyHz     upper search bound; the track is sampled at 10 Hz, so
%                      any motion frequency above its 5 Hz Nyquist limit is
%                      indistinguishable from an aliased one and is excluded
%   minAngleDeg        lower bound for the precession half angle
%   maxAngleDeg        upper bound for precession half angle and swing amplitude
%   swingDecisionDeg   historical amplitude reference used in plots
%   swingPromotionF    optional F3 cutoff for H1 to H3, defaults to parsimonyF

if nargin < 1 || isempty(overrides)
    overrides = struct();
end
cfg = struct( ...
    'tipLeverM', 0.75, ...
    'lambdaM', 3e8 / 14.8e9, ...
    'dwellSec', 32e-3, ...
    'dwellSubSamples', 8, ...
    'sampleStepSec', 0.1, ...
    'minFrequencyHz', 0.30, ...
    'maxFrequencyHz', 4.90, ...
    'minAngleDeg', 0.0, ...
    'maxAngleDeg', 25.0, ...
    'prfHz', 1000, ...
    'peakFftPoints', 512, ...
    'parsimonyF', 10, ...
    'spinFitRmsMaxHz', 10, ...
    'swingDecisionDeg', 0.20);
names = fieldnames(overrides);
for k = 1:numel(names)
    cfg.(names{k}) = overrides.(names{k});
end
end
