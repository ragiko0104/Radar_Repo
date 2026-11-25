function rvRes = rvpath_run(scene, options)
%RVPATH_RUN 测距测速链路（RV_path）纯计算五步调度总入口。
%
% 功能说明：
%   接收由 rv_isar_build_scene 构建的统一场景结构体 scene，按顺序调用 5 个纯计算函数：
%     1. rvpath_derive_waveform  : 六段 LFM 派生波形参数与快/慢时间轴；
%     2. rvpath_generate_echo    : fastTimeExact 通用散射体回波生成与下变频；
%     3. rvpath_pulse_compress   : 在线/理想时域 LFM 本振与六段频域脉冲压缩；
%     4. rvpath_cfar_peaks       : 六段 MEAN-CFAR 检测与精确寻峰；
%     5. rvpath_solve_rv         : 三段一致性匹配、距速解算与偏置校准。
%   注：本函数不执行任何绘图；可视化请调用 functions/plot/ 或 code/plot_rv_isar_results.m，
%       指标表格评估请调用 functions/assess/ 或 code/assess_rv_isar_results.m。
%
% 输入参数：
%   scene   - 1×1 struct，由 rv_isar_build_scene 产出的统一场景结构体，至少包含：
%     .cfgRV       1×1 struct，RV 波形与一致性门限配置
%     .constants   1×1 struct，物理常量（c、Mach）
%     .switches    1×1 struct，计算开关（rvLfmSpectrumMode）
%     .radarState  1×1 struct，雷达状态
%     .truth       1×1 struct，目标初始径向真值（R0/v_radial0/a_radial）
%     .targets     Ntarget×1 struct，统一散射体目标数组
%   options - (可选) 1×1 struct，可选控制项：
%     .echoOptions   1×1 struct，透传给 rvpath_generate_echo 的 options
%     .cfarOverrides 1×1 struct，透传给 rvpath_cfar_peaks 的 cfarOverrides
%     .solveOptions  1×1 struct，透传给 rvpath_solve_rv 的 options（如 printSummary）
%
% 输出参数：
%   rvRes   - 1×1 struct，RV 链路纯计算输出容器：
%     .cfg              1×1 struct，含全部派生波形字段的 RV 配置（wf）
%     .echo             1×1 struct，射频回波容器（total Na×N、groups 等）
%     .echoBaseband     Na×N complex double，去载频后的基带回波
%     .pulseCompression 1×1 struct，六段频域脉压结果（compressed 为 6×1 cell）
%     .cfar             1×1 struct，六段 CFAR 阈值与精确寻峰结果
%     .measurement      1×1 struct，距速解算、匹配、偏置校准与误差统计结果

if nargin < 2 || isempty(options)
    options = struct();
end

echoOpts      = get_opt(options, 'echoOptions', struct());
cfarOverrides = get_opt(options, 'cfarOverrides', struct());
solveOpts     = get_opt(options, 'solveOptions', struct('printSummary', false));

%% 步骤 1：六段 LFM 派生波形参数与时间轴
wf = rvpath_derive_waveform(scene.cfgRV);

%% 步骤 2：六段 LFM 回波生成（固定 fastTimeExact）与下变频
[rvEcho, rvEchoBaseband] = rvpath_generate_echo( ...
    scene.targets, scene.radarState, wf, scene.constants.c, echoOpts);

%% 步骤 3：在线/理想时域 LFM 本振与六段频域脉冲压缩
spectrumMode = 'floating';
if isfield(scene, 'switches') && isfield(scene.switches, 'rvLfmSpectrumMode') && ...
        ~isempty(scene.switches.rvLfmSpectrumMode)
    spectrumMode = scene.switches.rvLfmSpectrumMode;
end
rvPulseCompression = rvpath_pulse_compress(rvEchoBaseband, wf, spectrumMode);

%% 步骤 4：六段 CFAR 检测与精确寻峰
rvCfar = rvpath_cfar_peaks(rvPulseCompression, wf, cfarOverrides);

%% 步骤 5：三段联立一致性匹配、距速解算与偏置校准
rvMeasurement = rvpath_solve_rv( ...
    rvPulseCompression, rvCfar, wf, scene.truth, scene.constants, solveOpts);

rvRes = struct( ...
    'cfg',              wf, ...
    'echo',             rvEcho, ...
    'echoBaseband',     rvEchoBaseband, ...
    'pulseCompression', rvPulseCompression, ...
    'cfar',             rvCfar, ...
    'measurement',      rvMeasurement);
end

function val = get_opt(s, name, defaultVal)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    val = s.(name);
else
    val = defaultVal;
end
end
