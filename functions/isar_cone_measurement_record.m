function rec = isar_cone_measurement_record(t, rangeM, velocityMps, options)
%ISAR_CONE_MEASUREMENT_RECORD Explicit measurement record for the chain.
%   REC = ISAR_CONE_MEASUREMENT_RECORD(T,R,V,OPTIONS) freezes the range and
%   velocity that a tracker would report, so that downstream processing
%   consumes a record instead of re-deriving kinematics from scene truth.
%
%   The scene truth may be used HERE, to generate the observations. Everything
%   after this point must go through REC: isar_cone_processing_prior does not
%   accept a truth argument, which is what keeps the processing side from
%   reading R or V back out of the simulation.
%
%   OPTIONS
%     Mode            "measured" (default) or "exact"
%     RangeSigmaM     range measurement noise, default 0.02 m
%     VelocitySigmaMps  velocity measurement noise, default 0.3 m/s
%     Seed            RNG seed for the noise draw, default 0
%
%   REC fields: t, rangeM, velocityMps, mode, rangeSigmaM,
%   velocitySigmaMps, seed.

if nargin < 4 || isempty(options)
    options = struct();
end
t = double(t(:));
rangeM = double(rangeM(:));
velocityMps = double(velocityMps(:));
if numel(t) ~= numel(rangeM) || numel(t) ~= numel(velocityMps)
    error('ISAR:Chain:RecordLength', ...
        't, rangeM and velocityMps must have the same number of samples.');
end
mode = lower(string(opt(options, 'Mode', "measured")));
if ~isscalar(mode) || ~ismember(mode, ["measured", "exact"])
    error('ISAR:Chain:RecordMode', ...
        'Mode must be "measured" or "exact", received %s.', mat2str(mode));
end
sigR = double(opt(options, 'RangeSigmaM', 0.02));
sigV = double(opt(options, 'VelocitySigmaMps', 0.3));
seed = double(opt(options, 'Seed', 0));
rec = struct('t', t, 'rangeM', rangeM, 'velocityMps', velocityMps, ...
    'mode', mode, 'rangeSigmaM', sigR, 'velocitySigmaMps', sigV, ...
    'seed', seed);
if mode == "measured"
    saved = rng();
    cleanup = onCleanup(@() rng(saved)); %#ok<NASGU>
    rng(seed);
    rec.rangeM = rangeM + sigR * randn(size(rangeM));
    rec.velocityMps = velocityMps + sigV * randn(size(velocityMps));
end
end

function v = opt(s, name, fallback)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = fallback;
end
end
