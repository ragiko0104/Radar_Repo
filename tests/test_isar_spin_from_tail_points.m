function test_isar_spin_from_tail_points()
% End-to-end spin identification, including a folded synthetic scene.
% One scene spins at 1 Hz with 1 Hz precession (inside the 1 kHz PRF
% unambiguous interval), the other does not spin at all and must be rejected
% by the nested no-spin hypothesis rather than assigned an arbitrary rate.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'functions'), fullfile(root, 'experiments'));
cases = table(["T1_fs1_fp1_a3"; "T2_fs0_fp3_a3"], [1; 0], [1; 3], [3; 3], ...
    [0; 0], [0; 0], 'VariableNames', {'caseId', 'spinFrequencyHz', ...
    'precessionFrequencyHz', 'precessionHalfAngleDeg', ...
    'swingAmplitudeDeg', 'swingFrequencyHz'});
tmp = [tempname, '_spin'];
cases(end+1,:) = {"T3_folded", 2, 1.7, 12, 0, 0};
cleanup = onCleanup(@() rmdir(tmp, 's')); %#ok<NASGU>
result = identify_isar_spin_from_tail_points(struct( ...
    'Cases', cases, 'OutputDirectory', tmp, 'GenerateFigures', false));
T = result.table;
assert(T.spinDetected(1), 'A 1 Hz spin must be detected.');
assert(abs(T.errSpinHz(1)) < 0.02, ...
    'Spin estimate error %.4f Hz is too large.', T.errSpinHz(1));
assert(~T.spinDetected(2), 'A non-spinning scene must be rejected.');
assert(T.estSpinHz(2) < 0.2, ...
    'Non-spinning scene returned %.3f Hz.', T.estSpinHz(2));
assert(abs(T.errFpHz(1)) < 0.01, 'Precession frequency error too large.');
assert(abs(T.errAlphaDeg(1)) < 0.2, 'Half angle error too large.');
assert(all(T.spinResidualDf == 188));
expectedF = (T.rssNoSpin - T.rssSpin) ./ (T.rssSpin / 188);
assert(all(abs(T.Fstat - expectedF) <= 1e-10 * max(1, abs(expectedF))));
assert(T.aliased(3) && T.spinDetected(3), ...
    'The empirical spin gate can accept folded synthetic tracks.');
fprintf('test_isar_spin_from_tail_points passed.\n');
end
