function [fig, savedPath] = plot_scatterer_geometry(sceneOrIsarScene, plotCfg, saveOpts)
%PLOT_SCATTERER_GEOMETRY 绘制目标散射体 3D 模型、雷达 LOS 方向与旋转轴示意图。
%
% 功能说明：
%   1. 若所选目标为 UAV 模型且以体坐标系视图查看，直接复用 lib 中的
%      uavmodel.plotGeometry 绘制完整无人机机身/旋翼几何与法向；
%   2. 对 legacyDisk / rigidPointCloud / point 或传入 isarScene 的惯性系视图，
%      绘制 t=0（首脉冲）相对质心的 3D 散射点分布、雷达 LOS 方向矢量及有效转轴矢量；
%   3. 支持通过 saveOpts 将图形导出至指定路径。
%
% 输入参数：
%   sceneOrIsarScene - 1×1 struct，由 rv_isar_build_scene 产出的 scene，
%                      或由 isarpath_build_scene 产出的 isarScene（或 result.isar.scene）
%   plotCfg          - (可选) 1×1 struct，绘图配置项：
%     .targetIndex      1×1 double，目标索引（默认 1）
%     .useUavLibPlot    1×1 logical，UAV 模型是否优先调用 uavmodel.plotGeometry（默认 true）
%     .timeSec          1×1 double，UAV 几何采样时刻 (s，默认 0)
%     .showLabels       1×1 logical，是否标注散射点编号与权重（默认 true）
%     .visible          char/string，图窗可见性 'on'（默认）或 'off'
%   saveOpts         - (可选) 1×1 struct，保存配置（enable/dir/fileName/format/dpi）
%
% 输出参数：
%   fig       - 1×1 Figure 句柄
%   savedPath - 1×1 string，导出文件路径（未保存时为 ""）

if nargin < 2 || isempty(plotCfg), plotCfg = struct(); end
if nargin < 3, saveOpts = struct(); end

targetIdx     = get_cfg(plotCfg, 'targetIndex', 1);
useUavLibPlot = logical(get_cfg(plotCfg, 'useUavLibPlot', true));
timeSec       = double(get_cfg(plotCfg, 'timeSec', 0));
showLabels    = logical(get_cfg(plotCfg, 'showLabels', true));
figVisible    = char(string(get_cfg(plotCfg, 'visible', 'on')));

% 解析目标、散射点相对坐标、LOS 与转轴
hasIsarScene = isfield(sceneOrIsarScene, 'targetsISAR') && isfield(sceneOrIsarScene, 'scattererState');
if hasIsarScene
    targetObj = sceneOrIsarScene.targetsISAR(targetIdx);
else
    targetObj = sceneOrIsarScene.targets(targetIdx);
end

if useUavLibPlot && strcmpi(string(targetObj.kind), "uav") && ...
        isfield(targetObj, 'payload') && isfield(targetObj.payload, 'model')
    fig = figure('Name', sprintf('UAV Scatterer Geometry (Target %d)', targetIdx), ...
        'Color', 'w', 'Visible', figVisible);
    ax = axes(fig);
    uavmodel.plotGeometry(targetObj.payload.model, timeSec, struct( ...
        'Axes', ax, 'ShowLabels', showLabels, 'ShowFrame', true, 'ShowRotorNormals', true));
    title(ax, sprintf('Target %d (%s) Geometry at t = %.3f s', ...
        targetIdx, targetObj.name, timeSec), 'Interpreter', 'none');
    savedPath = save_output_artifact(fig, "figure", saveOpts, ...
        sprintf('scatterer_geometry_target%02d', targetIdx));
    return;
end

% 通用 3D 散射点 + LOS + 转轴绘制
if hasIsarScene
    st = sceneOrIsarScene.scattererState{targetIdx};
    ptsRel = st.position_I(:, :, 1) - st.centerPosition_I(1, :);
    alphaMag = abs(st.alpha(:, 1));
    geom = sceneOrIsarScene.crossRangeGeometry(targetIdx);
    losHat = extract_vec(geom, {'losHat_I', 'losHat', 'initialLosHat'}, [1, 0, 0]);
    rotAxis = extract_vec(geom, {'effectiveRotationVector_I', 'omega_I', 'rotationAxis_I'}, [0, 0, 1]);
    frameLabel = 'Inertial / Radar Frame (m, relative to center)';
else
    if isfield(targetObj.payload, 'positionTemplate_B')
        ptsRel = targetObj.payload.positionTemplate_B;
        alphaMag = abs(targetObj.payload.alpha0(:));
    elseif isfield(targetObj.payload, 'model') && isfield(targetObj.payload.model, 'scatterers')
        ptsRel = targetObj.payload.model.scatterers.positionTemplate_B;
        alphaMag = abs(targetObj.payload.model.scatterers.baseComplexCoefficient(:));
    else
        ptsRel = [0, 0, 0];
        alphaMag = 1;
    end
    losHat = [1, 0, 0];
    if isfield(sceneOrIsarScene, 'truth') && isfield(sceneOrIsarScene.truth, 'losHat_R')
        losHat = sceneOrIsarScene.truth.losHat_R;
    end
    rotAxis = [0, 1, 0];
    if isfield(sceneOrIsarScene, 'cfgISAR') && isfield(sceneOrIsarScene.cfgISAR, 'omegaBody_B')
        rotAxis = sceneOrIsarScene.cfgISAR.omegaBody_B;
    end
    frameLabel = 'Body Frame (m)';
end

scaleLen = max([0.5, max(vecnorm(ptsRel, 2, 2)) * 1.25]);
losHat = losHat(:).' ./ max(norm(losHat), eps);
if norm(rotAxis) > 0
    rotAxisHat = rotAxis(:).' ./ norm(rotAxis);
else
    rotAxisHat = [0, 0, 0];
end

fig = figure('Name', sprintf('Scatterer Geometry (Target %d)', targetIdx), ...
    'Color', 'w', 'Visible', figVisible);
ax = axes(fig);
hold(ax, 'on');
grid(ax, 'on');

markerSizes = 45 + 110 * (alphaMag ./ max([alphaMag(:); eps]));
scatter3(ax, ptsRel(:, 1), ptsRel(:, 2), ptsRel(:, 3), markerSizes, alphaMag, ...
    'filled', 'MarkerEdgeColor', 'k');
plot3(ax, 0, 0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 1.5);

quiver3(ax, 0, 0, 0, scaleLen * losHat(1), scaleLen * losHat(2), scaleLen * losHat(3), ...
    0, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.4);
legendEntries = {'Scatterers', 'Center', 'Radar +LOS'};
if any(rotAxisHat ~= 0)
    quiver3(ax, 0, 0, 0, scaleLen * rotAxisHat(1), scaleLen * rotAxisHat(2), ...
        scaleLen * rotAxisHat(3), 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.4);
    legendEntries{end+1} = 'Rotation Axis';
end

if showLabels && size(ptsRel, 1) <= 16
    for k = 1:size(ptsRel, 1)
        text(ax, ptsRel(k, 1), ptsRel(k, 2), ptsRel(k, 3), ...
            sprintf('  P%d (|a|=%.2f)', k, alphaMag(k)), 'FontSize', 9);
    end
end

axis(ax, 'equal');
view(ax, 38, 24);
xlabel(ax, 'X (m)');
ylabel(ax, 'Y (m)');
zlabel(ax, 'Z (m)');
legend(ax, legendEntries, 'Location', 'best');
title(ax, sprintf('Target %d (%s) — %s', targetIdx, targetObj.name, frameLabel), ...
    'Interpreter', 'none');

savedPath = save_output_artifact(fig, "figure", saveOpts, ...
    sprintf('scatterer_geometry_target%02d', targetIdx));
end

function val = get_cfg(s, name, defaultVal)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)), val = s.(name); else, val = defaultVal; end
end

function vec = extract_vec(s, candidateNames, defaultVec)
vec = defaultVec;
if ~isstruct(s), return; end
for k = 1:numel(candidateNames)
    nm = candidateNames{k};
    if isfield(s, nm) && ~isempty(s.(nm)) && numel(s.(nm)) == 3
        v = double(s.(nm)(:).');
        if all(isfinite(v)) && norm(v) > 0
            vec = v;
            return;
        end
    end
end
end
