function test_isar_golden_lo_offsets()
%TEST_ISAR_GOLDEN_LO_OFFSETS Per-target LO offsets must come from the golden.
%   The ten targets alternate 20 MHz / 45 MHz. Both beatOffsetHz and
%   loStartSec must use the same recorded value for the same target, and the
%   zero-based golden target id must map to the one-based MATLAB index.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'functions'));
gcfg = isar_golden_target_config();
cfg = struct('propagationSpeed', gcfg.propagationSpeed, ...
    'kr', gcfg.chirpRateHzPerS);
n = numel(gcfg.rangeM);
assert(n == 10, 'Expected ten golden targets, found %d.', n);
assert(isequal(gcfg.targetId, 0:9), 'Golden target ids must be zero-based 0..9.');
assert(max(abs(gcfg.beatOffsetHz - repmat([20e6, 45e6], 1, 5))) == 0, ...
    'Recorded bias must alternate 20/45 MHz.');
assert(any(gcfg.beatOffsetHz == 45e6), 'The 45 MHz branch must be present.');

t = (0:0.1:2)';
rec = isar_cone_measurement_record(t, ...
    gcfg.rangeM(1) + gcfg.velocityMps(1) * t, ...
    repmat(gcfg.velocityMps(1), size(t)), struct('Mode', "exact"));
fprintf('%6s %10s %10s %14s %16s\n', 'id0', 'rangeM', 'goldenBias', ...
    'usedBias', 'loStartSec');
for it = 1:n
    bias = gcfg.beatOffsetHz(it);
    prior = isar_cone_processing_prior(rec, 1, gcfg.accelerationMps2(it), ...
        bias, cfg);
    assert(prior.beatOffsetHz == gcfg.beatOffsetHz(it), ...
        'Target %d used a beat offset that is not the recorded one.', it - 1);
    expectedLo = 2 * rec.rangeM(1) / gcfg.propagationSpeed ...
        - gcfg.beatOffsetHz(it) / gcfg.chirpRateHzPerS;
    assert(abs(prior.loStartSec - expectedLo) < 1e-15, ...
        'Target %d LO start does not use the same recorded bias.', it - 1);
    fprintf('%6d %10.0f %10.0f %14.0f %16.9f\n', it - 1, rec.rangeM(1), ...
        gcfg.beatOffsetHz(it), prior.beatOffsetHz, prior.loStartSec);
end

% mixing a 20 MHz target with a 45 MHz offset must be detectable
lo20 = 2 * gcfg.rangeM(1) / gcfg.propagationSpeed - 20e6 / gcfg.chirpRateHzPerS;
lo45 = 2 * gcfg.rangeM(1) / gcfg.propagationSpeed - 45e6 / gcfg.chirpRateHzPerS;
delta = abs(lo45 - lo20);
assert(delta > 1e-7, 'Wrong-bias mixing must produce a visible LO shift.');
fprintf(['LO shift between 20 and 45 MHz: %.6f ms; the per-target bias is ' ...
    'therefore not interchangeable.\n'], delta * 1e3);
fprintf('test_isar_golden_lo_offsets passed.\n');
end
