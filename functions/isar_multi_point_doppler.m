function f = isar_multi_point_doppler(t, bodyPoints, p, cfg, fold)
%ISAR_MULTI_POINT_DOPPLER Dwell-averaged, PRF-folded Doppler of each point.
%   F = ISAR_MULTI_POINT_DOPPLER(T,BODYPOINTS,P,CFG,FOLD) returns an
%   numel(T)-by-numel(points) matrix. P is [fs fp alphaDeg ampDeg fw psiP
%   psiW psiS] (see isar_cone_point_kinematics).
%
%   Two effects separate this model from the instantaneous Doppler:
%     * the CPI is a 32 ms / 32 pulse coherent dwell, so one CPI reports the
%       dwell-averaged rate, not its value at a single instant;
%     * the intra-CPI Doppler is only unambiguous within [-prfHz/2,prfHz/2).
%       A base scatterer at 0.3 m radius reaches about 131 Hz per Hz of spin,
%       so at fs = 3.5 Hz it can leave that interval and fold. FOLD (default
%       true) applies the modulo, which is what the estimator must reproduce.
%
%   Range migration inside one dwell is not modelled: at the rates here it is
%   a few centimetres, far below the range-gate width.

if nargin < 5 || isempty(fold)
    fold = true;
end
t = double(t(:));
nPoint = size(bodyPoints, 1);
sub = cfg.dwellSubSamples;
offset = ((0:sub - 1) + 0.5) / sub * cfg.dwellSec - cfg.dwellSec / 2;
f = zeros(numel(t), nPoint);
for k = 1:nPoint
    acc = zeros(numel(t), 1);
    for s = 1:sub
        [~, rate] = isar_cone_point_kinematics(t + offset(s), bodyPoints(k, :), p, cfg);
        acc = acc + rate;
    end
    f(:, k) = 2 / cfg.lambdaM * acc / sub;
end
if fold
    f = isar_fold_doppler(f, cfg.prfHz);
end
end
