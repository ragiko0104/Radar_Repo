%% plot_rv_isar_results.m
% RV_ISAR_Path 独立绘图总入口脚本（计算与可视化彻底解耦）。
%
% 使用前提：
%   - 只画散射体几何/转动动图：工作区只需 scene（必要时本脚本自动调 isarpath_build_scene）；
%   - 画 RV 脉压/时频图或 ISAR HRRP/成像图：工作区需有 main.m / rv_isar_run 产出的 result。
% 使用方式：
%   1. 在顶部统一设置图片保存目录 saveDir、是否保存 enableSave、目标编号及各图开关；
%   2. 可整脚本运行，也可使用 MATLAB 分节运行（Ctrl+Enter）只绘制单张图。

if isempty(which('plot_isar_images'))
    run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'startup_project.m'));
end

%% ==================== 顶部统一保存与绘图配置入口 ====================
projectRoot = fileparts(fileparts(mfilename('fullpath')));
saveDir     = fullfile(projectRoot, 'results', 'figures'); % 图片/动图保存目录（可自行修改）
enableSave  = false;                                       % true: 自动导出图片/GIF；false: 仅屏幕显示
saveFormat  = 'png';                                       % 'png' | 'pdf' | 'fig' | 'jpg'
saveDpi     = 300;

targetIndex   = 1;     % 单目标观察图（几何、转动动图、HRRP 矩阵）的目标索引
targetIndices = [];    % 多目标 ISAR 成像图的目标索引（空表示全部目标）
rvPulseIndex  = 1;     % RV 时频图/脉压图的脉冲索引
figVisible    = 'on';  % 交互看图用 'on'，批处理仅导出可用 'off'

saveOptsBase = struct( ...
    'enable', enableSave, ...
    'dir',    saveDir, ...
    'format', saveFormat, ...
    'dpi',    saveDpi);

% 各图独立开关（按需设为 true / false）
doPlotScattererGeometry = true;
doPlotMotionAnimation   = true;
doPlotRvTimeFrequency   = true;
doPlotRvPulseCompress   = true;
doPlotIsarHrrpMatrix    = true;
doPlotIsarImages        = true;

%% ==================== 1. 散射体 3D 模型、LOS 与转轴示意图 ====================
if doPlotScattererGeometry
    if exist('result', 'var') && isfield(result, 'ranISAR') && result.ranISAR
        geomInput = result.isar.scene;
    elseif exist('isarScene', 'var') && isstruct(isarScene)
        geomInput = isarScene;
    elseif exist('scene', 'var') && isstruct(scene)
        geomInput = isarpath_build_scene(scene);
    else
        error('RV_ISAR_Path:MissingPlotInput', '请先运行 main 或 rv_isar_build_scene 生成 scene / result。');
    end
    [figGeom, geomSavedPath] = plot_scatterer_geometry(geomInput, struct( ...
        'targetIndex', targetIndex, 'showLabels', true, 'visible', figVisible), saveOptsBase);
end

%% ==================== 2. 目标慢时间 3D 转动动图（可选导出 GIF）====================
if doPlotMotionAnimation
    if exist('result', 'var') && isfield(result, 'ranISAR') && result.ranISAR
        animInput = result.isar.scene;
    elseif exist('isarScene', 'var') && isstruct(isarScene)
        animInput = isarScene;
    elseif exist('scene', 'var') && isstruct(scene)
        animInput = isarpath_build_scene(scene);
    else
        error('RV_ISAR_Path:MissingPlotInput', '请先生成 scene / result。');
    end
    [figAnim, animSavedPath] = plot_target_motion_animation(animInput, struct( ...
        'targetIndex', targetIndex, 'frameStep', 1, 'frameDelaySec', 0.08, ...
        'showTrails', true, 'visible', figVisible), saveOptsBase);
end

%% ==================== 3. RV 去载频回波时频图 ====================
if doPlotRvTimeFrequency && exist('result', 'var') && isfield(result, 'ranRV') && result.ranRV
    [figRvTf, rvTfSavedPath] = plot_rv_time_frequency(result.rv, struct( ...
        'pulseIndex', rvPulseIndex, 'visible', figVisible), saveOptsBase);
end

%% ==================== 4. RV 六段脉冲压缩与 CFAR 阈值图 ====================
if doPlotRvPulseCompress && exist('result', 'var') && isfield(result, 'ranRV') && result.ranRV
    cLight = 3e8;
    if exist('scene', 'var') && isfield(scene, 'constants')
        cLight = scene.constants.c;
    end
    [figRvPc, rvPcSavedPath] = plot_rv_pulse_compression(result.rv, struct( ...
        'pulseIndex', rvPulseIndex, 'c', cLight, 'visible', figVisible), saveOptsBase);
end

%% ==================== 5. ISAR 单目标 HRRP-慢时间矩阵与 RCMC 位移 ====================
if doPlotIsarHrrpMatrix && exist('result', 'var') && isfield(result, 'ranISAR') && result.ranISAR
    [figHrrp, hrrpSavedPath] = plot_isar_hrrp_matrix(result.isar, struct( ...
        'targetIndex',    targetIndex, ...
        'stages',         ["beforeDispersion", "afterDispersion", "afterRcmc", "afterPhaseComp"], ...
        'showRcmcCurve',  true, ...
        'dynamicRangeDb', 40, ...
        'visible',        figVisible), saveOptsBase);
end

%% ==================== 6. ISAR 逐目标二维成像结果 ====================
if doPlotIsarImages && exist('result', 'var') && isfield(result, 'ranISAR') && result.ranISAR
    [figIsar, isarSavedPath] = plot_isar_images(result.isar, struct( ...
        'targetIndices',  targetIndices, ...
        'dynamicRangeDb', 30, ...
        'colormap',       'turbo', ...
        'visible',        figVisible), saveOptsBase);
end
