function folded = isar_fold_doppler(hz, prfHz)
% Principal Doppler interval [-PRF/2, PRF/2), including both Nyquist ties.
folded = hz - prfHz .* floor(hz ./ prfHz + 0.5);
end
