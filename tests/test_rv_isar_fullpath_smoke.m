%TEST_RV_ISAR_FULLPATH_SMOKE Verify external simCfg across runMode = "rv" / "isar" / "both".
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

% Helper base external config (2 targets, 4-scatterer legacyDisk, plots off)
baseCfg = rv_isar_default_config();
baseCfg.targets.targetIndices            = [1 2];
baseCfg.scatterer.modelTypes             = "legacyDisk";
baseCfg.isar.Na                          = 16;
baseCfg.isar.echoMode                    = "stopAndGo";
baseCfg.switches.uavPlotGeometryEnabled  = false;
baseCfg.switches.rvShowEchoTimeFrequency = false;
baseCfg.switches.rvCfarPlotEnabled       = false;
baseCfg.switches.isarPlotEnabled         = false;
baseCfg.switches.isarQualityPrintEnabled = false;

% 2. External runMode = "rv" (only RV_path runs)
clearvars -except repoRoot baseCfg
simCfg = baseCfg;
simCfg.switches.runMode = "rv";
run(fullfile(repoRoot,'code','multiuav_iasr_fullpath.m'));
assert(exist('rvMeasurement','var') == 1, 'runMode="rv" must produce rvMeasurement.');
assert(exist('isarPost','var') == 0, 'runMode="rv" must NOT run ISAR_path.');
assert(all(rvMeasurement.matchedCount == 2));
assert(rvMeasurement.maxRawRangeRelativeErrorPercent < 1.0);
assert(rvMeasurement.maxRawVelocityRelativeErrorPercent < 2.0);
fprintf('  [PASS] runMode="rv" (2/2 matched, raw R err %.4f%%, raw V err %.4f%%)\n', ...
    rvMeasurement.maxRawRangeRelativeErrorPercent, ...
    rvMeasurement.maxRawVelocityRelativeErrorPercent);

% 3. External runMode = "isar" (only ISAR_path runs with truth priors)
clearvars -except repoRoot baseCfg
simCfg = baseCfg;
simCfg.switches.runMode = "isar";
run(fullfile(repoRoot,'code','multiuav_iasr_fullpath.m'));
assert(exist('rvMeasurement','var') == 0, 'runMode="isar" must NOT run RV_path.');
assert(exist('isarPost','var') == 1 && exist('isarQuality','var') == 1);
assert(numel(isarPost.valid) == 2 && all(isarPost.valid));
assert(height(isarQuality.summary) == 2 && all(isfinite(isarQuality.summary.ImageEntropy)));
fprintf('  [PASS] runMode="isar" (2/2 valid ISAR images)\n');

% 4. External runMode = "both" (chained RV -> ISAR)
clearvars -except repoRoot baseCfg
simCfg = baseCfg;
simCfg.switches.runMode = "both";
run(fullfile(repoRoot,'code','multiuav_iasr_fullpath.m'));
assert(exist('rvMeasurement','var') == 1 && exist('isarPost','var') == 1);
assert(all(rvMeasurement.matchedCount == 2));
assert(rvMeasurement.maxCalibratedRangeRelativeErrorPercent < 1.0);
assert(rvMeasurement.maxCalibratedVelocityRelativeErrorPercent < 1.0);
assert(numel(isarPost.valid) == 2 && all(isarPost.valid));
assert(all(isfinite(isarQuality.summary.Hrrf3dBAfterM)));
fprintf('  [PASS] runMode="both" (2/2 RV matched -> 2/2 ISAR valid)\n');

fprintf('RV + ISAR external-config smoke test passed across all 3 run modes.\n');
