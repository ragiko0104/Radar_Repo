%TEST_ISARPATH_MODULAR 功能与一致性测试：ISAR_path 纯计算八步函数与 isarpath_run。

clear;
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(thisFile));
run(fullfile(repoRoot, 'startup_project.m'));

%% 构造 2 目标 legacyDisk 测试场景（stopAndGo 快时间冻结，Na=16，纯计算）
simCfg = rv_isar_default_config();
simCfg.switches.runMode      = "isar";
simCfg.targets.targetIndices = [1 2];
simCfg.scatterer.modelTypes  = "legacyDisk";
simCfg.isar.Na               = 16;
simCfg.isar.echoMode         = "stopAndGo";
scene = rv_isar_build_scene(simCfg);

fprintf('=== 1. 逐项单元测试 ISAR_path 八个纯计算步骤函数 ===\n');

% Step 1: isarpath_build_scene
isarScene = isarpath_build_scene(scene);
p = isarScene.sharedParameters;
assert(isarScene.nTargets == 2 && isarScene.priorSource == "truth");
assert(isequal(size(isarScene.fastTimeSec), [1, isarScene.NHigh]));
assert(isequal(size(isarScene.pulseTimeSec), [16, 1]));
assert(isequal(size(isarScene.commonRangeAxisM), [1, p.hrrfBins]));
assert(numel(isarScene.crossAxisM) == 2 && numel(isarScene.scattererState) == 2);
fprintf('  [PASS] Step 1: isarpath_build_scene (N=%d, Na=%d, NHigh=%d)\n', ...
    isarScene.nTargets, isarScene.cfg.Na, isarScene.NHigh);

% Step 2: isarpath_generate_echo
[isarEcho, isarEchoInfo] = isarpath_generate_echo(isarScene);
assert(isequal(size(isarEcho.total), [16, isarScene.NHigh]));
assert(numel(isarEcho.perTarget) == 2);
assert(all(isfinite(isarEcho.total(:))) && any(abs(isarEcho.total(:)) > 0));
assert(strcmp(isarEchoInfo.mode, 'stopAndGo'));
fprintf('  [PASS] Step 2: isarpath_generate_echo (size=%dx%d, mode=%s)\n', ...
    size(isarEcho.total, 1), size(isarEcho.total, 2), isarEchoInfo.mode);

% Step 3: isarpath_dechirp
[isarBeat, isarBeatInfo] = isarpath_dechirp(isarEcho, isarScene);
assert(numel(isarBeat.perTarget) == 2);
assert(size(isarBeat.total, 1) == 16 && abs(isarBeat.sampleRateHz - isarScene.fsScopeHz) < 1e-6);
assert(isstruct(isarBeatInfo));
fprintf('  [PASS] Step 3: isarpath_dechirp (beat=%dx%d, fs=%.1f MHz)\n', ...
    size(isarBeat.total, 1), size(isarBeat.total, 2), isarBeat.sampleRateHz / 1e6);

% Step 4: isarpath_window_and_dispersion (target 1)
beat1  = isarBeat.perTarget(1).signal;
prior1 = isarScene.priors(1);
[segComp1, segRaw1, starts1] = isarpath_window_and_dispersion( ...
    beat1, isarScene.cfg, prior1, isarScene.pulseTimeSec, p);
assert(isequal(size(segComp1), [16, p.timeWindowSamples]));
assert(isequal(size(segRaw1),  [16, p.timeWindowSamples]));
assert(isequal(size(starts1),  [16, 1]));
fprintf('  [PASS] Step 4: isarpath_window_and_dispersion (segment=%dx%d)\n', ...
    size(segComp1, 1), size(segComp1, 2));

% Step 5: isarpath_form_hrrp
[hrrp1, cropCentres1, hrrpBefore1] = isarpath_form_hrrp( ...
    segComp1, isarScene.cfg, prior1, isarScene.pulseTimeSec, p, segRaw1);
assert(isequal(size(hrrp1), [16, p.hrrfBins]));
assert(isequal(size(hrrpBefore1), [16, p.hrrfBins]));
assert(isequal(size(cropCentres1), [16, 1]));
fprintf('  [PASS] Step 5: isarpath_form_hrrp (HRRP=%dx%d)\n', size(hrrp1, 1), size(hrrp1, 2));

% Step 6: isarpath_align_profiles (RCMC)
[aligned1, shifts1, score1] = isarpath_align_profiles(hrrp1, p);
assert(isequal(size(aligned1), [16, p.hrrfBins]));
assert(isequal(size(shifts1), [16, 1]));
assert(isequal(size(score1), [16, 2 * p.rcmcSearchRadiusBins + 1]));
fprintf('  [PASS] Step 6: isarpath_align_profiles (RCMC, max|shift|=%.2f bins)\n', max(abs(shifts1)));

% Step 7: isarpath_phase_compensate (VARIANCE_UNIT + MODEL_CV)
[phaseComp1, phaseCorr1, phaseInfo1] = isarpath_phase_compensate( ...
    aligned1, p, isarScene.cfg, prior1, isarScene.pulseTimeSec, shifts1);
assert(isequal(size(phaseComp1), [16, p.hrrfBins]));
assert(isequal(size(phaseCorr1), [16, 1]));
assert(phaseInfo1.mode == "VARIANCE_UNIT");

pModelCv = p;
pModelCv.phaseCompensationMode = "MODEL_CV";
[phaseCompCv, phaseCorrCv, phaseInfoCv] = isarpath_phase_compensate( ...
    aligned1, pModelCv, isarScene.cfg, prior1, isarScene.pulseTimeSec, shifts1);
assert(isequal(size(phaseCompCv), [16, p.hrrfBins]));
assert(max(abs(abs(phaseCorrCv) - 1)) < 1e-12 && phaseInfoCv.mode == "MODEL_CV");
fprintf('  [PASS] Step 7: isarpath_phase_compensate (VARIANCE_UNIT & MODEL_CV)\n');

% Step 8: isarpath_azimuth_fft
[img1, valid1] = isarpath_azimuth_fft(phaseComp1, p);
assert(isequal(size(img1), [16 * p.slowTimeZeroPaddingFactor, p.hrrfBins]));
assert(valid1 && all(isfinite(img1(:))));
fprintf('  [PASS] Step 8: isarpath_azimuth_fft (image=%dx%d, valid=%d)\n', ...
    size(img1, 1), size(img1, 2), valid1);

% 转发兼容聚合接口一致性验证：isarpath_range_profiles 与 isarpath_phase_image
pQual = p;
pQual.qualityEnabled = true;
[hrrpWrap, hrrpBeforeWrap, centresWrap] = isarpath_range_profiles( ...
    segRaw1, isarScene.cfg, prior1, isarScene.pulseTimeSec, pQual);
assert(max(abs(hrrpWrap(:) - hrrp1(:))) < 1e-12);
assert(max(abs(hrrpBeforeWrap(:) - hrrpBefore1(:))) < 1e-12);
assert(isequal(centresWrap, cropCentres1));

[phaseWrap, imgWrap, corrWrap, infoWrap] = isarpath_phase_image(aligned1, p);
assert(max(abs(phaseWrap(:) - phaseComp1(:))) < 1e-12);
assert(max(abs(imgWrap(:) - img1(:))) < 1e-12);
assert(max(abs(corrWrap(:) - phaseCorr1(:))) < 1e-12);
assert(infoWrap.mode == phaseInfo1.mode);

priorRefCv = isarpath_prior_phase( ...
    isarScene.cfg, prior1, isarScene.pulseTimeSec, pModelCv, shifts1, "MODEL_CV");
[phaseWrapCv, imgWrapCv, corrWrapCv] = isarpath_phase_image(aligned1, pModelCv, priorRefCv);
[imgCvDirect, ~] = isarpath_azimuth_fft(phaseCompCv, pModelCv);
assert(max(abs(phaseWrapCv(:) - phaseCompCv(:))) < 1e-12);
assert(max(abs(corrWrapCv(:) - phaseCorrCv(:))) < 1e-12);
assert(max(abs(imgWrapCv(:) - imgCvDirect(:))) < 1e-12);
fprintf('  [PASS] 兼容转发验证: isarpath_range_profiles & isarpath_phase_image (diff=0)\n');

%% 2. 端到端 isarpath_run 与 multiuav_isarpath.m 一致性验证
fprintf('\n=== 2. 端到端 isarpath_run 与 multiuav_isarpath.m 一致性校验 ===\n');
isarRes = isarpath_run(scene);
assert(numel(isarRes.post.valid) == 2 && all(isarRes.post.valid));
assert(numel(isarRes.post.hrrfBeforeDispersion) == 2);

cfg_ISARpath                = scene.cfgISAR;
targets                     = scene.targets;
radarState                  = scene.radarState;
R_isar                      = scene.truth.R0;
v_isar                      = scene.truth.v_radial0;
initialLosMotionDiagnostics = scene.truth.motionDiagnostics;
run(fullfile(repoRoot, 'code', 'multiuav_isarpath.m'));

diffEcho  = max(abs(isarRes.echo.total(:) - isarEcho.total(:)));
diffBeat  = max(abs(isarRes.beat.total(:) - isarBeat.total(:)));
diffHrrp  = max(cellfun(@(a, b) max(abs(a(:) - b(:))), isarRes.post.hrrfModelCrop, isarPost.hrrfModelCrop));
diffImg   = max(cellfun(@(a, b) max(abs(a(:) - b(:))), isarRes.post.image, isarPost.image));
assert(diffEcho < 1e-12 && diffBeat < 1e-12 && diffHrrp < 1e-12 && diffImg < 1e-12);
fprintf('  [PASS] 端到端一致：maxDiff(echo)=%.2e, maxDiff(beat)=%.2e, maxDiff(HRRP)=%.2e, maxDiff(image)=%.2e\n', ...
    diffEcho, diffBeat, diffHrrp, diffImg);

%% 3. 负拍频偏置（右侧本振）、可配置 RCMC 半径与 MODEL_CA 回归
fprintf('\n=== 3. 验证本振左右偏置、可配置 RCMC 半径与 MODEL_CA ===\n');
cfgRightLo = simCfg;
cfgRightLo.isar.frequencyBiasHz = [-20e6; -45e6];
cfgRightLo.isar.postprocess.rcmcSearchRadiusBins = 32;
sceneRightLo = rv_isar_build_scene(cfgRightLo);
assert(isequal(sceneRightLo.cfgISAR.frequencyBiasHz, [-20e6; -45e6]));
assert(sceneRightLo.isarSettings.rcmcSearchRadiusBins == 32);

isarResRight = isarpath_run(sceneRightLo);
tauEcho = 2 .* sceneRightLo.truth.R0 ./ sceneRightLo.constants.c;
assert(all(isarResRight.scene.referenceStartTimeSec > tauEcho));
assert(isarResRight.post.parameters.rcmcSearchRadiusBins == 32);
assert(size(isarResRight.post.rcmcCorrelationScore{1}, 2) == 2 * 32 + 1);
assert(all(isarResRight.post.valid));
fprintf('  [PASS] 负拍频偏置（右侧本振）与可配置 rcmcSearchRadiusBins=32 验证通过\n');

pModelCa = scene.isarSettings;
pModelCa.phaseCompensationMode = "MODEL_CA";
priorCa = isarRes.scene.priors(1);
priorCa.accelerationMps2 = scene.truth.a_radial(1);
[phaseCompCa, phaseCorrCa, phaseInfoCa] = isarpath_phase_compensate( ...
    isarRes.post.hrrfAligned{1}, pModelCa, isarRes.cfg, priorCa, ...
    isarRes.scene.pulseTimeSec, isarRes.post.rangeMigrationBins(1, :).');
assert(phaseInfoCa.mode == "MODEL_CA");
assert(abs(phaseInfoCa.accelerationUsedMps2 - priorCa.accelerationMps2) < 1e-12);
assert(max(abs(abs(phaseCorrCa) - 1)) < 1e-12 && all(isfinite(phaseCompCa(:))));
fprintf('  [PASS] isarpath_prior_phase (MODEL_CA) 回归验证通过\n');

% 3.3 自定义转轴 constantBodyRate（绕体坐标系 +z 轴）与 coneSwingSpin 场景构建检查
cfgRotZ = simCfg;
cfgRotZ.isar.bodyRotationMode = "constantBodyRate";
cfgRotZ.isar.omegaBody_B      = [0, 0, 2*pi];
sceneRotZ     = rv_isar_build_scene(cfgRotZ);
isarSceneRotZ = isarpath_build_scene(sceneRotZ);
assert(isequal(sceneRotZ.cfgISAR.omegaBody_B, [0, 0, 2*pi]));
assert(isequal(isarSceneRotZ.targetsISAR(1).payload.motion.omegaBody_B, [0, 0, 2*pi]));

cfgMicro = simCfg;
cfgMicro.isar.bodyRotationMode       = "coneSwingSpin";
cfgMicro.isar.spinFrequencyHz        = 2;
cfgMicro.isar.precessionFrequencyHz  = 1.5;
cfgMicro.isar.precessionHalfAngleDeg = 8;
cfgMicro.isar.swingFrequencyHz       = 1.2;
cfgMicro.isar.swingAngleMinDeg       = -6;
cfgMicro.isar.swingAngleMaxDeg       = 6;
sceneMicro     = rv_isar_build_scene(cfgMicro);
isarSceneMicro = isarpath_build_scene(sceneMicro);
assert(isarSceneMicro.microMotion.mode == "coneSwingSpin");
assert(abs(isarSceneMicro.microMotion.precessionHalfAngleDeg - 8) < 1e-12);
assert(abs(isarSceneMicro.microMotion.swingFrequencyHz - 1.2) < 1e-12);
assert(abs(isarSceneMicro.microMotion.swingAngleMaxDeg - 6) < 1e-12);
fprintf('  [PASS] 自定义转轴 constantBodyRate([0,0,2*pi]) 与 coneSwingSpin 参数传递验证通过\n');

fprintf('\n全部测试通过：ISAR_path 纯计算模块化验证完成。\n');
