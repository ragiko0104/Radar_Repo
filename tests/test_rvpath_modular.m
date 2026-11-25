%TEST_RVPATH_MODULAR 功能与一致性测试：RV_path 纯计算五步函数、rvpath_run 与独立绘图。

clear;
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(thisFile));
run(fullfile(repoRoot, 'startup_project.m'));

%% 构造 2 目标 legacyDisk 测试场景（fastTimeExact，纯计算无绘图开关）
simCfg = rv_isar_default_config();
simCfg.switches.runMode      = "rv";
simCfg.targets.targetIndices = [1 2];
simCfg.scatterer.modelTypes  = "legacyDisk";
scene = rv_isar_build_scene(simCfg);

fprintf('=== 1. 逐项单元测试 RV_path 五个纯计算步骤函数 ===\n');

% Step 1: rvpath_derive_waveform
wf = rvpath_derive_waveform(scene.cfgRV);
assert(numel(wf.f0_v) == 6 && numel(wf.f1_v) == 6 && numel(wf.fc_v) == 6);
assert(numel(wf.k_v) == 6 && numel(wf.delta_tau) == 6);
assert(isequal(size(wf.fastTimeSec), [1, wf.N]));
assert(isequal(size(wf.slowTimeSec), [wf.Na, 1]));
assert(wf.IF_FITACC == false);
fprintf('  [PASS] Step 1: rvpath_derive_waveform (N=%d, Na=%d)\n', wf.N, wf.Na);

% Step 2: rvpath_generate_echo（默认关闭分组与分块 + 显式开启分组/分块一致性对比）
[rvEcho, rvEchoBaseband] = rvpath_generate_echo( ...
    scene.targets, scene.radarState, wf, scene.constants.c);
assert(isequal(size(rvEcho.total), [wf.Na, wf.N]));
assert(isequal(size(rvEchoBaseband), [wf.Na, wf.N]));
assert(isempty(rvEcho.groups), '默认 returnGroups=false 时 groups 应为空。');
assert(all(isfinite(rvEchoBaseband(:))) && any(abs(rvEchoBaseband(:)) > 0));
[rvEchoGrouped, ~] = rvpath_generate_echo( ...
    scene.targets, scene.radarState, wf, scene.constants.c, ...
    struct('returnGroups', true, 'fastTimeChunkSize', 4096, 'scattererChunkSize', 2));
assert(numel(rvEchoGrouped.groups) == 4);
assert(max(abs(rvEcho.total(:) - rvEchoGrouped.total(:))) < 1e-12, ...
    '关闭分组/分块与开启分组/分块的总回波必须数值一致。');
fprintf('  [PASS] Step 2: rvpath_generate_echo (size=%dx%d, default groups=0, grouped diff=%.2e)\n', ...
    size(rvEcho.total, 1), size(rvEcho.total, 2), max(abs(rvEcho.total(:) - rvEchoGrouped.total(:))));

% Step 3: rvpath_pulse_compress (floating + fixedPoint + idealTime via signal.LFM_wave)
pcFloat = rvpath_pulse_compress(rvEchoBaseband, wf, 'floating');
pcFixed = rvpath_pulse_compress(rvEchoBaseband, wf, 'fixedPoint');
pcIdeal = rvpath_pulse_compress(rvEchoBaseband, wf, 'idealTime');
assert(numel(pcFloat.compressed) == 6 && numel(pcFixed.compressed) == 6 && ...
    numel(pcIdeal.compressed) == 6);
assert(isequal(size(pcFloat.compressed{1}), [wf.Na, wf.N]));
assert(isequal(size(pcIdeal.compressed{1}), [wf.Na, wf.N]));
assert(strcmp(pcFloat.localOscillatorMode, 'floating'));
assert(strcmp(pcFixed.localOscillatorMode, 'fixedPoint'));
assert(strcmp(pcIdeal.localOscillatorMode, 'idealTime'));
assert(strcmp(pcIdeal.localOscillatorGenerator, 'signal.LFM_wave'));
assert(~isempty(pcFixed.localOscillatorConfigData{1}));
assert(pcIdeal.localOscillatorInfo(1).pulseSampleCount == round(wf.N_lfm));
fprintf('  [PASS] Step 3: rvpath_pulse_compress (floating & fixedPoint & idealTime)\n');

% Step 4: rvpath_cfar_peaks
rvCfar = rvpath_cfar_peaks(pcFloat, wf);
rvCfarIdeal = rvpath_cfar_peaks(pcIdeal, wf);
assert(numel(rvCfar.threshold) == 6 && numel(rvCfar.refinedPeaks) == 6);
assert(isequal(size(rvCfar.meanPeakPower), [wf.Na, 6]));
assert(all(rvCfar.meanPeakPower(1, :) > 0) && all(rvCfarIdeal.meanPeakPower(1, :) > 0));
fprintf('  [PASS] Step 4: rvpath_cfar_peaks (6 segments refined, floating & idealTime)\n');

% Step 5: rvpath_solve_rv (snrPrecision + legacyResolution + idealTime LO)
measSnr = rvpath_solve_rv(pcFloat, rvCfar, wf, scene.truth, scene.constants, ...
    struct('printSummary', false));
wfLegacyGate = wf;
wfLegacyGate.rvConsistency.mode = 'legacyResolution';
measLegacy = rvpath_solve_rv(pcFloat, rvCfar, wfLegacyGate, scene.truth, scene.constants, ...
    struct('printSummary', false));
measIdeal = rvpath_solve_rv(pcIdeal, rvCfarIdeal, wf, scene.truth, scene.constants, ...
    struct('printSummary', false));
assert(all(measSnr.matchedCount == 2) && all(measLegacy.matchedCount == 2) && ...
    all(measIdeal.matchedCount == 2));
assert(measSnr.maxCalibratedRangeRelativeErrorPercent < 1.0);
assert(measSnr.maxCalibratedVelocityRelativeErrorPercent < 1.0);
assert(measIdeal.maxCalibratedRangeRelativeErrorPercent < 1.0);
assert(measIdeal.maxCalibratedVelocityRelativeErrorPercent < 1.0);
fprintf(['  [PASS] Step 5: rvpath_solve_rv (snrPrecision, legacyResolution & idealTime: ' ...
    'rawR=%.4f%% calR=%.4f%% rawV=%.4f%% calV=%.4f%%)\n'], ...
    measIdeal.maxRawRangeRelativeErrorPercent, ...
    measIdeal.maxCalibratedRangeRelativeErrorPercent, ...
    measIdeal.maxRawVelocityRelativeErrorPercent, ...
    measIdeal.maxCalibratedVelocityRelativeErrorPercent);

%% 2. 端到端 rvpath_run 与兼容薄脚本 multiuav_rvpath.m 一致性校验
fprintf('\n=== 2. 端到端 rvpath_run 与 multiuav_rvpath.m 一致性校验 ===\n');
rvRes = rvpath_run(scene);

cfg_RVpath        = scene.cfgRV;
c                 = scene.constants.c;
Mach              = scene.constants.Mach;
targets           = scene.targets;
radarState        = scene.radarState;
R0                = scene.truth.R0;
v_radial0         = scene.truth.v_radial0;
a_radial          = scene.truth.a_radial;
rvLfmSpectrumMode = scene.switches.rvLfmSpectrumMode;

run(fullfile(repoRoot, 'code', 'multiuav_rvpath.m'));

diffEcho = max(abs(rvRes.echo.total(:) - rvEcho.total(:)));
diffPc   = max(cellfun(@(a, b) max(abs(a(:) - b(:))), ...
    rvRes.pulseCompression.compressed, rvPulseCompression.compressed));
diffCalR = max(abs(rvRes.measurement.calibratedRangeM(:) - rvMeasurement.calibratedRangeM(:)));
diffCalV = max(abs(rvRes.measurement.calibratedVelocityMps(:) - rvMeasurement.calibratedVelocityMps(:)));

assert(diffEcho < 1e-12 && diffPc < 1e-12 && diffCalR < 1e-12 && diffCalV < 1e-12);
fprintf('  [PASS] 端到端一致：maxDiff(echo)=%.2e, maxDiff(PC)=%.2e, maxDiff(calR)=%.2e, maxDiff(calV)=%.2e\n', ...
    diffEcho, diffPc, diffCalR, diffCalV);
fprintf('\n全部测试通过：RV_path 纯计算模块化验证完成。\n');
