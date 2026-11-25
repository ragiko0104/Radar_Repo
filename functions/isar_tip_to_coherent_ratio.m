function [ratio, info] = isar_tip_to_coherent_ratio(H, gate)
%ISAR_TIP_TO_COHERENT_RATIO Tip-channel energy over coherent total energy.
%   [R,INFO] = ISAR_TIP_TO_COHERENT_RATIO(H,GATE) reports
%
%       ratio = sum_{gate} |H(:,gate,2)|^2 / sum_{gate} |H(:,gate,1)|^2
%
%   i.e. the energy of the simulator's isolated tip channel divided by the
%   energy of the COHERENT TOTAL FIELD. The denominator is NOT the sum of the
%   per-scatterer energies: the scatterers interfere, so the ratio is not
%   bounded by one and is not a "fraction of the tip". Values above one are
%   expected when the total field cancels inside the gate; do not clamp.
%
%   When the denominator is exactly zero the ratio is NaN and INFO.valid is
%   false, with INFO.reason set to 'zero-coherent-energy'.
%
%   This is a diagnostic only. It must not feed gate selection, track
%   extraction or the mode decision; nothing in those paths calls it.
%
%   INFO: valid, tipEnergy, coherentEnergy, reason

if ndims(H) < 3 || size(H, 3) < 2
    error('ISAR:Chain:RatioChannels', ...
        'H needs at least two channels (1 coherent total, 2 tip).');
end
gate = gate(:);
if any(gate < 1) || any(gate > size(H, 2))
    error('ISAR:Chain:RatioGate', 'Gate indices outside the range axis.');
end
tipEnergy = sum(sum(abs(double(H(:, gate, 2))) .^ 2, 1));
coherentEnergy = sum(sum(abs(double(H(:, gate, 1))) .^ 2, 1));
info = struct('valid', true, 'tipEnergy', tipEnergy, ...
    'coherentEnergy', coherentEnergy, 'reason', 'ok');
if coherentEnergy == 0
    ratio = NaN;
    info.valid = false;
    info.reason = 'zero-coherent-energy';
else
    ratio = tipEnergy / coherentEnergy;
end
end
