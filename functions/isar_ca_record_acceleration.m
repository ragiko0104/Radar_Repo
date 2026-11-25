function aHat = isar_ca_record_acceleration(rec)
% Fit the same frozen range record that supplies the processing prior.
t = rec.t(:);
r = rec.rangeM(:);
if numel(r) ~= numel(t) || any(~isfinite([t;r])) || numel(unique(t)) < 3
    error('ISAR:Chain:FitRecord', 'Need finite ranges at at least three distinct epochs.');
end
X = [ones(size(t)), t, 0.5*t.^2];
coeff = X \ r;
aHat = coeff(3);
end
