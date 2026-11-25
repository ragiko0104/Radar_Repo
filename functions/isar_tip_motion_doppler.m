function f = isar_tip_motion_doppler(t, theta, cfg, biasHz)
%ISAR_TIP_MOTION_DOPPLER Dwell-averaged tip Doppler for a parameter vector.
%   F = ISAR_TIP_MOTION_DOPPLER(T,THETA,CFG) returns the slow-time Doppler
%   frequency (Hz) that one CPI would report for the instantaneous tip
%   motion. THETA is [fp alphaDeg ampDeg fsw psiP psiW] as produced by
%   isar_fit_tip_motion_mode. Optional biasHz is added before PRF folding.
%
%   A CPI integrates 32 pulses over CFG.dwellSec, and the extracted peak is
%   a dwell average of the instantaneous Doppler, not its value at a single
%   instant. Ignoring that operator leaves a residual tone at fp with about
%   1.4 % of the modulation amplitude, which biases the precession half angle
%   low by the same amount and lets richer hypotheses absorb the artefact.
%   The average is evaluated numerically so that the harmonics of fp and fsw
%   are attenuated correctly.

t = double(t(:));
if nargin < 4
    biasHz = 0;
elseif ~isscalar(biasHz)
    if ~isvector(biasHz) || numel(biasHz) ~= numel(t)
        error('RadarProject:TipMode:InvalidBias', 'Supply a scalar bias or one bias per time sample.');
    end
    biasHz = biasHz(:);
end
sub = cfg.dwellSubSamples;
offset = ((0:sub - 1) + 0.5) / sub * cfg.dwellSec - cfg.dwellSec / 2;
acc = zeros(numel(t), 1);
for k = 1:sub
    [~, rate] = isar_tip_motion_kinematics(t + offset(k), theta(1), theta(2), ...
        theta(3), theta(4), theta(5), theta(6), cfg);
    acc = acc + rate;
end
f = isar_fold_doppler(2 / cfg.lambdaM * acc / sub + biasHz, cfg.prfHz);
end
