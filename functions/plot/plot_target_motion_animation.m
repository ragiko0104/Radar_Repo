function [fig, savedPath] = plot_target_motion_animation(isarSceneOrRes, plotCfg, saveOpts)
%PLOT_TARGET_MOTION_ANIMATION 绘制目标散射点在慢时间内的 3D 转动动图（可选导出 GIF）。
%
% 功能说明：
%   基于 isarScene.scattererState{k}.position_I（Ns×3×Na）与质心位置 centerPosition_I
%   （Na×3），逐脉冲展示目标散射点绕质心的三维微动轨迹、雷达 LOS 方向与转轴方向；
%   当 saveOpts.enable == true 时自动将动画保存为 GIF 文件。
%
% 输入参数：
%   isarSceneOrRes - 1×1 struct，isarpath_build_scene 输出 isarScene 或 isarpath_run 输出 isarRes
%   plotCfg        - (可选) 1×1 struct，动图配置项：
%     .targetIndex    1×1 double，目标索引（默认 1）
%     .frameStep      1×1 double，帧步长（默认 1，逐脉冲播放）
%     .frameDelaySec  1×1 double，帧间隔时间 (s，默认 0.08；visible='off' 时不暂停)
%     .showTrails     1×1 logical，是否绘制散射点历史轨迹线（默认 true）
%     .viewAngles     1×2 double，三维视角 [az, el]（默认 [38, 24]）
%     .visible        char/string，图窗可见性 'on'（默认）或 'off'
%   saveOpts       - (可选) 1×1 struct，保存配置（enable/dir/fileName，导出为 .gif）
%
% 输出参数：
%   fig       - 1×1 Figure 句柄（停留在最后一帧）
%   savedPath - 1×1 string，导出的 GIF 完整路径（未保存时为 ""）

if nargin < 2 || isempty(plotCfg), plotCfg = struct(); end
if nargin < 3, saveOpts = struct(); end

if isfield(isarSceneOrRes, 'scene') && isfield(isarSceneOrRes.scene, 'scattererState')
    isarScene = isarSceneOrRes.scene;
else
    isarScene = isarSceneOrRes;
end

targetIdx  = get_cfg(plotCfg, 'targetIndex', 1);
frameStep  = max(1, round(double(get_cfg(plotCfg, 'frameStep', 1))));
frameDelay = double(get_cfg(plotCfg, 'frameDelaySec', 0.08));
showTrails = logical(get_cfg(plotCfg, 'showTrails', true));
viewAngles = double(get_cfg(plotCfg, 'viewAngles', [38, 24]));
figVisible = char(string(get_cfg(plotCfg, 'visible', 'on')));

st         = isarScene.scattererState{targetIdx};
pulseTimes = isarScene.pulseTimeSec(:);
nPulses    = size(st.position_I, 3);
nScat      = size(st.position_I, 1);

% 相对质心的 3D 散射点轨迹 (Ns × 3 × Na)
relPos = zeros(size(st.position_I));
for pIdx = 1:nPulses
    relPos(:, :, pIdx) = st.position_I(:, :, pIdx) - st.centerPosition_I(pIdx, :);
end
maxRadius = max([0.2, max(vecnorm(reshape(permute(relPos, [1 3 2]), [], 3), 2, 2))]);
axisLim   = maxRadius * 1.35;

geom   = isarScene.crossRangeGeometry(targetIdx);
losHat = [1, 0, 0];
if isfield(geom, 'losHat_I') && numel(geom.losHat_I) == 3
    losHat = double(geom.losHat_I(:).') ./ max(norm(geom.losHat_I), eps);
end

fig = figure('Name', sprintf('Target %d Slow-Time Rotation', targetIdx), ...
    'Color', 'w', 'Visible', figVisible);
ax = axes(fig);
hold(ax, 'on');
grid(ax, 'on');
axis(ax, 'equal');
xlim(ax, [-axisLim, axisLim]);
ylim(ax, [-axisLim, axisLim]);
zlim(ax, [-axisLim, axisLim]);
view(ax, viewAngles(1), viewAngles(2));
xlabel(ax, 'X_I - X_c (m)');
ylabel(ax, 'Y_I - Y_c (m)');
zlabel(ax, 'Z_I - Z_c (m)');

% 固定参考：质心与雷达 LOS 方向
plot3(ax, 0, 0, 0, 'k+', 'MarkerSize', 9, 'LineWidth', 1.4);
quiver3(ax, 0, 0, 0, axisLim * 0.85 * losHat(1), axisLim * 0.85 * losHat(2), ...
    axisLim * 0.85 * losHat(3), 0, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.35);

hPts = scatter3(ax, relPos(:, 1, 1), relPos(:, 2, 1), relPos(:, 3, 1), ...
    70, 'b', 'filled', 'MarkerEdgeColor', 'k');

hTrails = gobjects(0, 1);
if showTrails && nScat <= 32
    hTrails = gobjects(nScat, 1);
    for sIdx = 1:nScat
        hTrails(sIdx) = plot3(ax, relPos(sIdx, 1, 1), relPos(sIdx, 2, 1), ...
            relPos(sIdx, 3, 1), '-', 'LineWidth', 0.9, 'Color', [0.35, 0.55, 0.85]);
    end
end
hTitle = title(ax, '');

% 准备 GIF 导出路径
enableSave = isstruct(saveOpts) && isfield(saveOpts, 'enable') && logical(saveOpts.enable);
savedPath  = "";
if enableSave
    if ~isfield(saveOpts, 'dir') || isempty(saveOpts.dir)
        error('RV_ISAR_Path:MissingSaveDir', 'saveOpts.enable=true 时必须提供 saveOpts.dir。');
    end
    outDir = char(string(saveOpts.dir));
    if exist(outDir, 'dir') ~= 7
        mkdir(outDir);
    end
    if isfield(saveOpts, 'fileName') && ~isempty(saveOpts.fileName)
        [~, baseName, ~] = fileparts(char(string(saveOpts.fileName)));
    else
        baseName = sprintf('target_motion_target%02d', targetIdx);
    end
    gifPath = fullfile(outDir, [baseName, '.gif']);
    savedPath = string(gifPath);
end

frameIndices = 1:frameStep:nPulses;
if frameIndices(end) ~= nPulses
    frameIndices(end+1) = nPulses;
end

for k = 1:numel(frameIndices)
    pIdx = frameIndices(k);
    set(hPts, 'XData', relPos(:, 1, pIdx), 'YData', relPos(:, 2, pIdx), 'ZData', relPos(:, 3, pIdx));
    if ~isempty(hTrails)
        for sIdx = 1:nScat
            set(hTrails(sIdx), ...
                'XData', squeeze(relPos(sIdx, 1, 1:pIdx)), ...
                'YData', squeeze(relPos(sIdx, 2, 1:pIdx)), ...
                'ZData', squeeze(relPos(sIdx, 3, 1:pIdx)));
        end
    end
    set(hTitle, 'String', sprintf('Target %d (%s) | Pulse %d/%d | t = %.3f ms', ...
        targetIdx, isarScene.targetsISAR(targetIdx).name, pIdx, nPulses, ...
        pulseTimes(pIdx) * 1e3), 'Interpreter', 'none');
    drawnow;

    if enableSave
        frame = getframe(fig);
        [imInd, cm] = rgb2ind(frame2im(frame), 256);
        if k == 1
            imwrite(imInd, cm, gifPath, 'gif', 'LoopCount', inf, 'DelayTime', max(frameDelay, 0.02));
        else
            imwrite(imInd, cm, gifPath, 'gif', 'WriteMode', 'append', 'DelayTime', max(frameDelay, 0.02));
        end
    elseif strcmpi(figVisible, 'on') && frameDelay > 0
        pause(frameDelay);
    end
end
end

function val = get_cfg(s, name, defaultVal)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)), val = s.(name); else, val = defaultVal; end
end
