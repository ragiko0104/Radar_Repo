function [result, scene] = rv_isar_run(simCfgOrScene, options)
%RV_ISAR_RUN RV 测距测速 + ISAR 成像全链路函数化总入口。
%
% 功能说明：
%   1. 接收配置结构体 simCfg（或已由 rv_isar_build_scene 构建好的 scene 结构体），
%      必要时调用 rv_isar_build_scene 完成场景构建；
%   2. 根据 scene.switches.runMode（"rv" | "isar" | "both"）顺序调用模块化链路函数：
%        - "rv"   : 仅调用 rvpath_run(scene, rvOptions)；
%        - "isar" : 仅调用 isarpath_run(scene, priorMotion, isarOptions)（使用真值先验）；
%        - "both" : 先调用 rvpath_run，再将校准后的距速估计作为先验送入 isarpath_run；
%   3. 返回结构化结果 result 与更新后的 scene，不污染调用方工作区。
%
% 输入参数：
%   simCfgOrScene - (可选) 1×1 struct，支持两种形式：
%                     (a) 仿真配置结构体 simCfg（完整或局部覆盖字段，格式同
%                         rv_isar_default_config）；
%                     (b) 已由 rv_isar_build_scene 生成的 scene 结构体。
%                   省略或为空时使用默认配置自动构建 scene。
%   options       - (可选) 1×1 struct，透传给子链路的可选控制项：
%     .rvOptions    1×1 struct，透传给 rvpath_run 的 options
%     .isarOptions  1×1 struct，透传给 isarpath_run 的 options
%
% 输出参数：
%   result - 1×1 struct，全链路运行结果容器，字段如下：
%     .runMode       1×1 string      实际执行的模式（"rv" | "isar" | "both"）
%     .ranRV         1×1 logical     本次是否执行了 RV 链路
%     .ranISAR       1×1 logical     本次是否执行了 ISAR 链路
%     .priorSource   1×1 string      ISAR 先验来源（"rvMeasurement" | "truth" | "none"）
%     .priorMotion   1×1 struct      传给 ISAR 的先验运动参数（未跑 ISAR 时仍记录可用先验）：
%       .rangeM        N×1 double    各目标先验初始斜距 (m)
%       .velocityMps   N×1 double    各目标先验径向速度 (m/s)
%     .rv            1×1 struct      rvpath_run 返回的 RV 链路结果（未运行 RV 时为空 struct()）
%     .isar          1×1 struct      isarpath_run 返回的 ISAR 链路结果（未运行 ISAR 时为空 struct()）
%   scene  - 1×1 struct，本次运行使用的完整场景结构体（详见 rv_isar_build_scene）。

if nargin < 2 || isempty(options)
    options = struct();
end

if nargin < 1 || isempty(simCfgOrScene)
    scene = rv_isar_build_scene();
elseif is_scene_struct(simCfgOrScene)
    scene = simCfgOrScene;
else
    scene = rv_isar_build_scene(simCfgOrScene);
end

rvOpts   = get_opt(options, 'rvOptions', struct());
isarOpts = get_opt(options, 'isarOptions', struct());

runRV   = logical(scene.switches.runRV);
runISAR = logical(scene.switches.runISAR);

rvRes   = struct();
isarRes = struct();

if runRV
    rvRes = rvpath_run(scene, rvOpts);
    scene.cfgRV = rvRes.cfg;
end

if runRV
    priorRangeM      = rvRes.measurement.calibratedRangeM(:);
    priorVelocityMps = rvRes.measurement.calibratedVelocityMps(:);
    priorSource      = "rvMeasurement";
else
    priorRangeM      = scene.truth.R0(:);
    priorVelocityMps = scene.truth.v_radial0(:);
    if runISAR
        priorSource = "truth";
    else
        priorSource = "none";
    end
end

priorMotion = struct( ...
    'rangeM',      priorRangeM, ...
    'velocityMps', priorVelocityMps);

if runISAR
    isarRes = isarpath_run(scene, priorMotion, isarOpts);
end

result = struct( ...
    'runMode',     string(scene.switches.runMode), ...
    'ranRV',       runRV, ...
    'ranISAR',     runISAR, ...
    'priorSource', priorSource, ...
    'priorMotion', priorMotion, ...
    'rv',          rvRes, ...
    'isar',        isarRes);
end

function tf = is_scene_struct(candidate)
tf = isstruct(candidate) && isscalar(candidate) && ...
    isfield(candidate, 'targets') && ...
    isfield(candidate, 'cfgRV') && ...
    isfield(candidate, 'cfgISAR') && ...
    isfield(candidate, 'radarState') && ...
    isfield(candidate, 'truth') && ...
    isfield(candidate, 'switches');
end

function val = get_opt(s, name, defaultVal)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    val = s.(name);
else
    val = defaultVal;
end
end
