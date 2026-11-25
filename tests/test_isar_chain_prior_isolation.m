function test_isar_chain_prior_isolation()
%TEST_ISAR_CHAIN_PRIOR_ISOLATION Processing priors must not read scene truth.
%   Truth may generate the echo and the measurement record; after that the
%   processing side may only consume the record. This test fixes the record
%   and asks whether changing the scoring truth can alter anything that is
%   handed to the compensator.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'functions'));
cfg = struct('propagationSpeed', 3e8, 'kr', 6e13);
t = (0:0.1:2)';
range0 = 10050;
v0 = 3468;
aTrue = -150;
truth = struct('rangeM', range0 + v0 * t + 0.5 * aTrue * t .^ 2, ...
    'velocityMps', v0 + aTrue * t, 'accelerationMps2', aTrue);

% ---- generation side (truth is allowed here)
rec = isar_cone_measurement_record(t, truth.rangeM, truth.velocityMps, ...
    struct('Mode', "measured", 'RangeSigmaM', 0.02, ...
    'VelocitySigmaMps', 0.3, 'Seed', 17));
X = [ones(numel(t), 1), t, 0.5 * t .^ 2];
aHat = (X \ rec.rangeM);
aHat = aHat(3);
bias = 20e6;

% ---- processing side: no truth argument exists
priorBefore = isar_cone_processing_prior(rec, 7, aHat, bias, cfg);

% change only the scoring truth
truth.accelerationMps2 = -999;
truth.rangeM = truth.rangeM + 1234;
truth.velocityMps = truth.velocityMps * 2;

priorAfter = isar_cone_processing_prior(rec, 7, aHat, bias, cfg);
assert(isequal(priorBefore, priorAfter), ...
    'Scoring truth must not reach the processing prior.');

% ---- the prior is the record, not the truth
assert(priorBefore.rangeM == rec.rangeM(7), 'Range must come from the record.');
assert(priorBefore.velocityMps == rec.velocityMps(7), ...
    'Velocity must come from the record.');
assert(priorBefore.accelerationMps2 == aHat, ...
    'Acceleration must be the fitted value.');
noisy = abs(priorBefore.rangeM - (range0 + v0 * t(7) + 0.5 * (-150) * t(7) ^ 2)) > 0;
assert(noisy, 'Measured record must not equal truth bit-for-bit.');

% ---- exact mode is a separately named control
recExact = isar_cone_measurement_record(t, truth.rangeM, truth.velocityMps, ...
    struct('Mode', "exact"));
assert(recExact.rangeM(1) == truth.rangeM(1), 'Exact mode returns truth.');

% ---- the prior rejects a truth-shaped struct masquerading as a record
threw = false;
try
    isar_cone_processing_prior(truth, 1, aHat, bias, cfg);
catch err
    threw = strcmp(err.identifier, 'ISAR:Chain:BadRecord');
end
assert(threw, 'A raw truth struct must be rejected as a record.');

fprintf(['test_isar_chain_prior_isolation passed (aHat %.3f, record mode %s).\n'], ...
    aHat, string(rec.mode));
end
