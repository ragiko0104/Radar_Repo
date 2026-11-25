%TEST_PLOT_AND_ASSESS_MODULES 验证独立评估函数（functions/assess）与绘图函数（functions/plot）。
%
% 验证项：
%   1. 纯计算链路 rv_isar_run 不自动创建任何图窗；
%   2. 全部 5 个独立评估函数（assess_rv_measurement / assess_isar_dispersion /
%      assess_isar_rcmc / assess_isar_image / assess_isar_scatterers）计算与 CSV 保存；
%   3. 全部 6 个独立绘图函数（plot_scatterer_geometry 含 legacyDisk & uav /
%      plot_target_motion_animation 含 GIF 导出 / plot_rv_time_frequency /
%      plot_rv_pulse_compression / plot_isar_hrrp_matrix / plot_isar_images）
%      无头（Visible='off'）渲染、PNG/GIF 导出及图窗释放；
%   4. 测试完成后自动清理临时导出目录。

clear;
close all;
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(thisFile));
run(fullfile(repoRoot, 'startup_project.m'));

tmpOutDir = fullfile(tempdir, ['rv_isar_plot_assess_test_' char(java.util.UUID.randomUUID)]);

try
%% 1. 跑一次 2 目标纯计算全链路，确认零弹窗
simCfg = rv_isar_default_config();
simCfg.switches.runMode      = "both";
simCfg.targets.targetIndices = [1 2];
simCfg.scatterer.modelTypes  = "legacyDisk";
simCfg.isar.Na               = 16;
simCfg.isar.echoMode         = "stopAndGo";

figCountBefore = numel(findall(groot, 'Type', 'figure'));
[result, scene] = rv_isar_run(simCfg);
figCountAfter = numel(findall(groot, 'Type', 'figure'));
assert(figCountAfter == figCountBefore, '纯计算 rv_isar_run 不得自动创建图窗。');
assert(isfield(result.isar.post, 'hrrfBeforeDispersion') && ...
    numel(result.isar.post.hrrfBeforeDispersion) == 2);
fprintf('  [PASS] 纯计算链路零弹窗，且保留 hrrfBeforeDispersion\n');

%% 2. 测试全部 5 个独立评估函数（不保存 & 保存 CSV）
saveTblOff = struct('enable', false, 'dir', tmpOutDir, 'format', 'csv');
saveTblOn  = struct('enable', true,  'dir', fullfile(tmpOutDir, 'assess'), 'format', 'csv');
assessCfg  = struct('printSummary', false, 'expectedPeakCount', 4);

[rvDet, rvTbl, p0] = assess_rv_measurement(result.rv, assessCfg, saveTblOff);
assert(p0 == "" && rvDet.matchedCount == 2 && height(rvTbl) == 2);
[~, ~, pRv] = assess_rv_measurement(result.rv, assessCfg, saveTblOn);
assert(isfile(pRv) && height(readtable(pRv)) == 2);

[dispDet, dispTbl, pDisp] = assess_isar_dispersion(result.isar, assessCfg, saveTblOn);
assert(numel(dispDet) == 2 && height(dispTbl) == 2 && isfile(pDisp));

[rcmcDet, rcmcTbl, pRcmc] = assess_isar_rcmc(result.isar, assessCfg, saveTblOn);
assert(numel(rcmcDet) == 2 && height(rcmcTbl) == 2 && isfile(pRcmc));

[imgDet, imgTbl, pImg] = assess_isar_image(result.isar, assessCfg, saveTblOn);
assert(numel(imgDet) == 2 && height(imgTbl) == 2 && all(isfinite(imgTbl.ImageEntropy)) && isfile(pImg));

[scatDet, scatTbl, pScat] = assess_isar_scatterers(result.isar, assessCfg, saveTblOn);
assert(numel(scatDet) == 2 && height(scatTbl) == 2 && isfile(pScat));
fprintf('  [PASS] 全部 5 个 assess_* 评估函数与 CSV 导出验证通过\n');

%% 3. 测试全部 6 个独立绘图函数（无头 Visible='off' + PNG/GIF 导出）
saveFigOn = struct('enable', true, 'dir', fullfile(tmpOutDir, 'figures'), ...
    'format', 'png', 'dpi', 150);

% 3.1 散射体几何（legacyDisk 惯性系 + UAV 体坐标系）
[f1, pf1] = plot_scatterer_geometry(result.isar.scene, ...
    struct('targetIndex', 1, 'visible', 'off'), saveFigOn);
assert(isgraphics(f1, 'figure') && isfile(pf1));
close(f1);

cfgUav = simCfg;
cfgUav.targets.targetIndices = 1;
cfgUav.scatterer.modelTypes  = "uav";
sceneUav = rv_isar_build_scene(cfgUav);
[f1u, pf1u] = plot_scatterer_geometry(sceneUav, ...
    struct('targetIndex', 1, 'visible', 'off'), ...
    struct('enable', true, 'dir', fullfile(tmpOutDir, 'figures'), ...
    'fileName', 'uav_geometry', 'format', 'png', 'dpi', 150));
assert(isgraphics(f1u, 'figure') && isfile(pf1u));
close(f1u);

% 3.2 慢时间 3D 转动动图（导出 GIF，步长 4 加快测试）
[f2, pf2] = plot_target_motion_animation(result.isar.scene, ...
    struct('targetIndex', 1, 'frameStep', 4, 'frameDelaySec', 0.05, 'visible', 'off'), ...
    struct('enable', true, 'dir', fullfile(tmpOutDir, 'figures')));
assert(isgraphics(f2, 'figure') && isfile(pf2) && endsWith(pf2, '.gif'));
close(f2);

% 3.3 RV 时频图
[f3, pf3] = plot_rv_time_frequency(result.rv, ...
    struct('pulseIndex', 1, 'visible', 'off'), saveFigOn);
assert(isgraphics(f3, 'figure') && isfile(pf3));
close(f3);

% 3.4 RV 3x2 脉压与 CFAR 图
[f4, pf4] = plot_rv_pulse_compression(result.rv, ...
    struct('pulseIndex', 1, 'c', scene.constants.c, 'visible', 'off'), saveFigOn);
assert(isgraphics(f4, 'figure') && isfile(pf4));
close(f4);

% 3.5 ISAR HRRP-慢时间矩阵
[f5, pf5] = plot_isar_hrrp_matrix(result.isar, ...
    struct('targetIndex', 1, 'visible', 'off'), saveFigOn);
assert(isgraphics(f5, 'figure') && isfile(pf5));
close(f5);

% 3.6 ISAR 二维成像（复用 plotting.draw_isar_img）
[f6, pf6] = plot_isar_images(result.isar, ...
    struct('targetIndices', [1 2], 'dynamicRangeDb', 30, 'visible', 'off'), saveFigOn);
assert(isgraphics(f6, 'figure') && isfile(pf6));
close(f6);

fprintf('  [PASS] 全部 6 个 plot_* 绘图函数（含 UAV 几何、GIF 动图、PNG 导出）验证通过\n');
fprintf('\n全部测试通过：独立评估模块与独立绘图模块验证完成。\n');
catch ME
    if exist(tmpOutDir, 'dir') == 7
        rmdir(tmpOutDir, 's');
    end
    rethrow(ME);
end
if exist(tmpOutDir, 'dir') == 7
    rmdir(tmpOutDir, 's');
end
