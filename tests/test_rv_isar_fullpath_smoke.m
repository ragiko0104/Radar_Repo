%TEST_RV_ISAR_FULLPATH_SMOKE Unattended end-to-end verification of RV + ISAR chain.
% Uses 4-scatterer legacyDisk targets with stopAndGo (fast-time frozen) mode.

clear;
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(thisFile));
run(fullfile(repoRoot,'startup_project.m'));

% 1. Verify no duplicate local functions shadow matlab_research_lib
localFuncs = dir(fullfile(repoRoot,'functions','**','*.m'));
localNames = string({localFuncs.name});
assert(isequal(localNames, "run_uav_existing_pipeline.m"), ...
    'Unexpected local functions found in functions/: %s', strjoin(localNames, ', '));

% 2. 10-target RV measurement verification (4-scatterer legacyDisk, stopAndGo)
FULLPATH_OVERRIDES = struct('rvEchoMode', "stopAndGo");
run_RVpath = true;
run_ISARpath = false;
targetModelTypes = repmat("legacyDisk", 10, 1);
rvShowEchoTimeFrequency = false;
rvCfarPlotEnabled = false;
UAV_PLOT_GEOMETRY_ENABLED = false;

run(fullfile(repoRoot,'code','multiuav_iasr_fullpath.m'));

assert(exist('rvMeasurement','var') == 1);
assert(all(rvMeasurement.matchedCount == 10), ...
    'Expected all 10 targets matched in RV path, got %d.', rvMeasurement.matchedCount(1));
assert(rvMeasurement.maxCalibratedRangeRelativeErrorPercent < 1.0, ...
    'RV calibrated range error too large: %.4f%%', ...
    rvMeasurement.maxCalibratedRangeRelativeErrorPercent);
assert(rvMeasurement.maxCalibratedVelocityRelativeErrorPercent < 1.0, ...
    'RV calibrated velocity error too large: %.4f%%', ...
    rvMeasurement.maxCalibratedVelocityRelativeErrorPercent);

% 3. 2-target RV + ISAR chained imaging & quality verification (4-scatterer, stopAndGo)
clear targetModelTypes FULLPATH_OVERRIDES rvMeasurement isarPost isarQuality rvEchoMode
FULLPATH_OVERRIDES = struct('targetIndices', [1 2], 'isarNa', 16, ...
    'rvEchoMode', "stopAndGo", 'isarEchoMode', "stopAndGo");
run_RVpath = true;
run_ISARpath = true;
targetModelTypes = repmat("legacyDisk", 2, 1);
rvShowEchoTimeFrequency = false;
rvCfarPlotEnabled = false;
ISAR_PLOT_ENABLED = false;
ISAR_QUALITY_PRINT_ENABLED = false;
UAV_PLOT_GEOMETRY_ENABLED = false;

run(fullfile(repoRoot,'code','multiuav_iasr_fullpath.m'));

assert(exist('isarPost','var') == 1 && exist('isarQuality','var') == 1);
assert(numel(isarPost.valid) == 2 && all(isarPost.valid), ...
    'Both targets should produce valid ISAR images.');
assert(height(isarQuality.summary) == 2);
assert(all(isfinite(isarQuality.summary.ImageEntropy)));
assert(all(isfinite(isarQuality.summary.Hrrf3dBAfterM)));

fprintf('RV + ISAR fullpath smoke test passed (10/10 RV matched, 2/2 stopAndGo ISAR valid).\n');
