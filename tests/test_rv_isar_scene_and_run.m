%TEST_RV_ISAR_SCENE_AND_RUN 功能与数值一致性测试：rv_isar_build_scene 与 rv_isar_run。
%
% 验证项：
%   1. rv_isar_build_scene：
%      - 默认/局部覆盖配置合并、派生波形参数与目标真值维度检查；
%      - 按需散射体实例化（legacyDisk 不实例化 uav、uav 正常实例化、异构组合与 char 输入）；
%      - 非法 runMode 与目标向量长度不匹配的错误处理。
%   2. rv_isar_run：
%      - runMode = "rv"   （传入 simCfg，仅跑 RV_path）；
%      - runMode = "isar" （传入预构建 scene，仅跑 ISAR_path 且使用真值先验）；
%      - runMode = "both" （串联 RV -> ISAR）；
%      - 工作区隔离性检查（不向调用方泄漏脚本临时变量）。
%   3. 与薄脚本 code/multiuav_iasr_fullpath.m 数值逐项对齐。

clear;
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(thisFile));
run(fullfile(repoRoot, 'startup_project.m'));

fprintf('=== 1. 测试 rv_isar_build_scene ===\n');

cfgLegacy = rv_isar_default_config();
cfgLegacy.switches.runMode      = 'rv';
cfgLegacy.targets.targetIndices = [1 2];
cfgLegacy.scatterer.modelTypes  = 'legacyDisk';
cfgLegacy.isar.Na               = 16;
cfgLegacy.isar.echoMode         = "stopAndGo";

sceneLegacy = rv_isar_build_scene(cfgLegacy);
assert(sceneLegacy.nTargets == 2);
assert(sceneLegacy.switches.runRV && ~sceneLegacy.switches.runISAR);
assert(isequal(size(sceneLegacy.truth.R0), [2, 1]));
assert(isequal(size(sceneLegacy.truth.v_radial0), [2, 1]));
assert(isequal(size(sceneLegacy.truth.a_radial), [2, 1]));
assert(isequal(size(sceneLegacy.targets), [2, 1]));
assert(all(sceneLegacy.scatterer.activeModelTypes == "legacyDisk"));
assert(~isfield(sceneLegacy.scatterer, 'uav'), ...
    'legacyDisk 模式下不应实例化 scatterer.uav。');
assert(abs(sceneLegacy.cfgISAR.B - 6e9) < 1e-6);
assert(abs(sceneLegacy.cfgISAR.fc - 15e9) < 1e-6);
assert(isequal(size(sceneLegacy.radarState.position_I([0; 1e-3])), [2, 3]));
fprintf('  [PASS] legacyDisk 场景构建与按需实例化检查\n');

% 1.1b lib 预设函数句柄 + 局部 overrides + 异构 cell 句柄检查
cfgPreset = cfgLegacy;
cfgPreset.scatterer.modelTypes = {@scatterermodel.configs.cfg_legacy_ZXT4; ...
    @scatterermodel.configs.cfg_rigid_cross4};
cfgPreset.scatterer.overrides  = {struct('radiusM', 0.45); struct()};
scenePreset = rv_isar_build_scene(cfgPreset);
assert(scenePreset.targets(1).kind == "legacyDisk" && ...
    abs(scenePreset.targets(1).metadata.radiusM - 0.45) < 1e-12);
assert(scenePreset.targets(2).kind == "rigidPointCloud" && ...
    scenePreset.targets(2).scattererCount == 4);
fprintf('  [PASS] lib 预设句柄（cfg_legacy_ZXT4 + cfg_rigid_cross4）及 overrides 检查\n');

cfgMixed = cfgLegacy;
cfgMixed.scatterer.modelTypes = ["legacyDisk"; "uav"];
sceneMixed = rv_isar_build_scene(cfgMixed);
assert(isfield(sceneMixed.scatterer.uav, 'model') && ...
    isfield(sceneMixed.scatterer.uav, 'config'));
assert(sceneMixed.targets(1).kind == "legacyDisk");
assert(sceneMixed.targets(2).kind == "uav");
fprintf('  [PASS] 异构模型 ["legacyDisk"; "uav"] 构建检查\n');

didThrow = false;
try
    badCfg = cfgLegacy;
    badCfg.switches.runMode = "invalid_mode";
    rv_isar_build_scene(badCfg);
catch ME
    didThrow = strcmp(ME.identifier, 'RV_ISAR_Path:InvalidRunMode');
end
assert(didThrow, '非法 runMode 应抛出 RV_ISAR_Path:InvalidRunMode。');

didThrow = false;
try
    badCfg = cfgLegacy;
    badCfg.targets.R0 = [10050; 20050];
    rv_isar_build_scene(badCfg);
catch ME
    didThrow = strcmp(ME.identifier, 'RV_ISAR_Path:TargetVectorLengthMismatch');
end
assert(didThrow, '目标向量长度不一致应抛出 RV_ISAR_Path:TargetVectorLengthMismatch。');
fprintf('  [PASS] 非法配置异常拦截检查\n');

fprintf('\n=== 2. 测试 rv_isar_run（rv / isar / both 三模式）===\n');

[resRV, sceneOutRV] = rv_isar_run(cfgLegacy);
assert(resRV.ranRV && ~resRV.ranISAR);
assert(resRV.priorSource == "rvMeasurement");
assert(isempty(fieldnames(resRV.isar)));
assert(all(resRV.rv.measurement.matchedCount == 2));
assert(resRV.rv.measurement.maxRawRangeRelativeErrorPercent < 1.0);
assert(resRV.rv.measurement.maxCalibratedRangeRelativeErrorPercent < 1.0);
assert(resRV.rv.measurement.maxCalibratedVelocityRelativeErrorPercent < 1.0);
assert(isfield(sceneOutRV.cfgRV, 'f0_v') && numel(sceneOutRV.cfgRV.f0_v) == 6);
assert(exist('rvEcho', 'var') == 0 && exist('rvMeasurement', 'var') == 0, ...
    'rv_isar_run 不得将脚本临时变量泄漏到调用方工作区。');
fprintf('  [PASS] runMode="rv"（2/2 匹配，校准后 R 误差 %.4f%%，V 误差 %.4f%%）\n', ...
    resRV.rv.measurement.maxCalibratedRangeRelativeErrorPercent, ...
    resRV.rv.measurement.maxCalibratedVelocityRelativeErrorPercent);

cfgIsarOnly = cfgLegacy;
cfgIsarOnly.switches.runMode = "isar";
sceneIsarOnly = rv_isar_build_scene(cfgIsarOnly);
[resISAR, ~] = rv_isar_run(sceneIsarOnly);
assert(~resISAR.ranRV && resISAR.ranISAR);
assert(resISAR.priorSource == "truth");
assert(isempty(fieldnames(resISAR.rv)));
assert(isequal(resISAR.priorMotion.rangeM, sceneIsarOnly.truth.R0));
assert(isequal(resISAR.priorMotion.velocityMps, sceneIsarOnly.truth.v_radial0));
assert(numel(resISAR.isar.post.valid) == 2 && all(resISAR.isar.post.valid));
[~, isarImgTbl] = assess_isar_image(resISAR.isar);
assert(height(isarImgTbl) == 2 && all(isfinite(isarImgTbl.ImageEntropy)));
assert(exist('isarPost', 'var') == 0 && exist('isarBeat', 'var') == 0, ...
    'rv_isar_run 不得将 ISAR 脚本临时变量泄漏到调用方工作区。');
fprintf('  [PASS] runMode="isar"（传入预构建 scene，2/2 有效成像）\n');

cfgBoth = cfgLegacy;
cfgBoth.switches.runMode = "both";
[resBoth, ~] = rv_isar_run(cfgBoth);
assert(resBoth.ranRV && resBoth.ranISAR);
assert(resBoth.priorSource == "rvMeasurement");
assert(all(resBoth.rv.measurement.matchedCount == 2));
assert(numel(resBoth.isar.post.valid) == 2 && all(resBoth.isar.post.valid));
fprintf('  [PASS] runMode="both"（RV -> ISAR 串联成功）\n');

fprintf('\n=== 3. 与薄脚本 multiuav_iasr_fullpath.m 数值逐项对齐 ===\n');
simCfg = cfgBoth;
run(fullfile(repoRoot, 'code', 'multiuav_iasr_fullpath.m'));

maxDiffCalR = max(abs(resBoth.rv.measurement.calibratedRangeM(:) - ...
    rvMeasurement.calibratedRangeM(:)));
maxDiffCalV = max(abs(resBoth.rv.measurement.calibratedVelocityMps(:) - ...
    rvMeasurement.calibratedVelocityMps(:)));
maxDiffImg1 = max(abs(resBoth.isar.post.image{1}(:) - isarPost.image{1}(:)));
maxDiffImg2 = max(abs(resBoth.isar.post.image{2}(:) - isarPost.image{2}(:)));

assert(maxDiffCalR < 1e-12, '校准距离不一致: %.3e', maxDiffCalR);
assert(maxDiffCalV < 1e-12, '校准速度不一致: %.3e', maxDiffCalV);
assert(maxDiffImg1 < 1e-12, '目标1 ISAR 图像不一致: %.3e', maxDiffImg1);
assert(maxDiffImg2 < 1e-12, '目标2 ISAR 图像不一致: %.3e', maxDiffImg2);

fprintf('  [PASS] 数值一致性对齐：maxDiff(R)=%.2e, maxDiff(V)=%.2e, maxDiff(Img)=%.2e\n', ...
    maxDiffCalR, maxDiffCalV, max(maxDiffImg1, maxDiffImg2));
fprintf('\n全部测试通过：rv_isar_build_scene 与 rv_isar_run 功能及数值一致性验证完成。\n');
