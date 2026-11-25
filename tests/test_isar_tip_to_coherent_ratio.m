function test_isar_tip_to_coherent_ratio()
%TEST_ISAR_TIP_TO_COHERENT_RATIO Behaviour of the tip/coherent-energy ratio.
%   The denominator is the coherent total field, so the ratio is not bounded
%   by one and must not be clamped; a zero denominator yields NaN plus an
%   explicit invalid marker. The quantity is diagnostic only.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'functions'));
nPulse = 8;
nRange = 12;
gate = 5:7;

% ---- 1. ordinary case: tip below the coherent total
H = complex(zeros(nPulse, nRange, 5));
H(:, gate, 1) = 2 + 0.5i;            % coherent total
H(:, gate, 2) = 1 + 0.2i;            % isolated tip
ratio = isar_tip_to_coherent_ratio(H, gate);
tipE = numel(H(:, gate, 2)) * abs(1 + 0.2i) ^ 2;
cohE = numel(H(:, gate, 1)) * abs(2 + 0.5i) ^ 2;
assert(abs(ratio - tipE / cohE) < 1e-12, 'Ordinary ratio is wrong.');
assert(ratio < 1, 'Tip below the total must give a ratio below one.');

% ---- 2. cancellation: the coherent field nearly cancels but the tip does not
H2 = complex(zeros(nPulse, nRange, 5));
H2(:, gate, 1) = 1 - 0.99;           % near cancellation of two scatterers
H2(:, gate, 2) = 1;
ratio2 = isar_tip_to_coherent_ratio(H2, gate);
assert(ratio2 > 1, 'A cancelling total field must give a ratio above one.');
assert(abs(ratio2 - 1 / 0.01 ^ 2) < 1e-6, 'Cancellation ratio is wrong.');

% ---- 3. zero denominator: explicit invalid marker, not a clamp
H3 = complex(zeros(nPulse, nRange, 5));
H3(:, gate, 1) = 0;
H3(:, gate, 2) = 1 + 1i;
[r3, info3] = isar_tip_to_coherent_ratio(H3, gate);
assert(isnan(r3), 'Zero coherent energy must give NaN.');
assert(~info3.valid && strcmp(info3.reason, 'zero-coherent-energy'), ...
    'Zero denominator must be flagged.');

% ---- 4. shape guards
threw = false;
try
    isar_tip_to_coherent_ratio(H(:, :, 1), gate);
catch err
    threw = strcmp(err.identifier, 'ISAR:Chain:RatioChannels');
end
assert(threw, 'A single-channel cube must be rejected.');
threw = false;
try
    isar_tip_to_coherent_ratio(H, [1 2 99]);
catch err
    threw = strcmp(err.identifier, 'ISAR:Chain:RatioGate');
end
assert(threw, 'A gate outside the range axis must be rejected.');

% ---- 5. the diagnostic must not feed gating, tracks or the decision
txt = fileread(fullfile(root, 'experiments', 'isar_chain_identify_ch1.m'));
hits = numel(regexp(txt, 'isar_tip_to_coherent_ratio', 'start'));
assert(hits == 1, ['The ratio must appear once, as a reporting call; found ' ...
    '%d occurrences.'], hits);
lines = splitlines(string(txt));
callers = lines(contains(lines, 'isar_tip_to_coherent_ratio'));
assert(isscalar(callers) && startsWith(strtrim(callers), 'tipRatio('), ...
    'The ratio must only be assigned to the tipRatio diagnostic.');
assignsGate = ~cellfun(@isempty, regexp(cellstr(lines), '^\s*(gate|yTrack|i0|b)\s*=', 'once'));
assert(~any(assignsGate & contains(lines, 'isar_tip_to_coherent_ratio')), ...
    'The ratio must not feed the gate, the track or the bin index.');
fprintf('test_isar_tip_to_coherent_ratio passed.\n');
end
