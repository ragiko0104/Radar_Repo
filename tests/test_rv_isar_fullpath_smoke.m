%TEST_RV_ISAR_FULLPATH_SMOKE Unattended end-to-end verification of RV + ISAR chain.
% RV always runs in fastTimeExact (required for intra-pulse Doppler velocity);
% ISAR uses stopAndGo (fast-time frozen) for fast smoke verification.

clear;
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(thisFile));
run(fullfile(repoRoot,'startup_project.m'));

% 1. Verify no duplicate local functions shadow matlab_research_lib
localFuncs = dir(fullfile(repoRoot,'functions','**','*.m'));
localStems = string(arrayfun(@(f) erase(f.name,'.m'), localFuncs, 'UniformOutput', false));
libFuncs = dir(fullfile(matlab_research_lib_dir,'**','*.m'));
libStems = string(arrayfun(@(f) erase(f.name,'.m'), libFuncs, 'UniformOutput', false));
shadowed = intersect(localStems, libStems);
assert(isempty(shadowed), ...
    'Local functions must not duplicate matlab_research_lib: %s', strjoin(shadowed, ', '));

% 2. 2-target chained RV (fastTimeExact) -> ISAR (stopAndGo) verification
FULLPATH_OVERRIDES = struct('targetIndices', [1 2], 'isarNa', 16, 'isarEchoMode', "stopAndGo");
run_RVpath = true;
run_ISARpath = true;
targetModelTypes = repmat("legacyDisk", 2, 1);
rvShowEchoTimeFrequency = false;
rvCfarPlotEnabled = false;
ISAR_PLOT_ENABLED = false;
ISAR_QUALITY_PRINT_ENABLED = false;
UAV_PLOT_GEOMETRY_ENABLED = false;

run(fullfile(repoRoot,'code','multiuav_iasr_fullpath.m'));

% Verify RV outputs (both raw and calibrated errors must be small under fastTimeExact)
assert(exist('rvMeasurement','var') == 1);
assert(all(rvMeasurement.matchedCount == 2), ...
    'Expected 2 targets matched in RV path, got %d.', rvMeasurement.matchedCount(1));
assert(rvMeasurement.maxRawRangeRelativeErrorPercent < 1.0, ...
    'RV raw range error too large: %.4f%%', ...
    rvMeasurement.maxRawRangeRelativeErrorPercent);
assert(rvMeasurement.maxRawVelocityRelativeErrorPercent < 2.0, ...
    'RV raw velocity error too large: %.4f%%', ...
    rvMeasurement.maxRawVelocityRelativeErrorPercent);
assert(rvMeasurement.maxCalibratedRangeRelativeErrorPercent < 1.0, ...
    'RV calibrated range error too large: %.4f%%', ...
    rvMeasurement.maxCalibratedRangeRelativeErrorPercent);
assert(rvMeasurement.maxCalibratedVelocityRelativeErrorPercent < 1.0, ...
    'RV calibrated velocity error too large: %.4f%%', ...
    rvMeasurement.maxCalibratedVelocityRelativeErrorPercent);

% Verify ISAR outputs and quality metrics
assert(exist('isarPost','var') == 1 && exist('isarQuality','var') == 1);
assert(numel(isarPost.valid) == 2 && all(isarPost.valid), ...
    'Both targets should produce valid ISAR images.');
assert(height(isarQuality.summary) == 2);
assert(all(isfinite(isarQuality.summary.ImageEntropy)));
assert(all(isfinite(isarQuality.summary.Hrrf3dBAfterM)));

fprintf('RV (fastTimeExact) + ISAR (stopAndGo) smoke test passed (2/2 matched, 2/2 valid).\n');
