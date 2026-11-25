function [isarBeat, isarBeatInfo] = isarpath_dechirp(isarEcho, isarScene, options)
%ISARPATH_DECHIRP 生成逐目标局部 LFM 本振并完成去斜、滤波与降采样，获得窄带拍频。
%
% 功能说明：
%   1. 使用各目标的先验本振起始时刻 referenceStartTimeSec，调用 lib 中的
%      radar.lfm_local_reference_beat 对逐目标高采样回波执行本振去斜；
%   2. 默认按 ANALYTIC_IF 模式完成镜像频率剔除、抗混叠低通滤波、交叠串扰隔离、
%      去直流及降采样至示波器/FPGA 侧采样率 fs_scope；
%   3. 校验输出的逐目标拍频维度与采样率一致性。
%
% 输入参数：
%   isarEcho  - 1×1 struct，由 isarpath_generate_echo 产出的高采样率回波
%               （要求包含 .perTarget(1..N).total，每项为 Na×NHigh complex）
%   isarScene - 1×1 struct，由 isarpath_build_scene 产出的 ISAR 场景结构体
%               （使用 cfg/fastTimeSec/referenceStartTimeSec/tauRefSec/fsScopeHz/dechirpOptions）
%   options   - (可选) 1×1 struct，用于覆盖 isarScene.dechirpOptions 中的字段
%               （如 Mode、RemoveMirrorFrequency、RemoveSamplingAlias、
%               RemoveTimeOverlapCrosstalk、RemoveDc 等）
%
% 输出参数：
%   isarBeat     - 1×1 struct，低采样率窄带拍频回波容器：
%     .total         Na×NBeat complex double，场景叠加拍频信号
%     .perTarget     N×1 struct，逐目标独立拍频（每项 .signal 为 Na×NBeat complex）
%     .sampleRateHz  1×1 double，输出采样率，应等于 fs_scope
%   isarBeatInfo - 1×1 struct，radar.lfm_local_reference_beat 返回的去斜与降采样元信息

cfgISAR = isarScene.cfg;
nTarget = isarScene.nTargets;
Na_ISAR = cfgISAR.Na;

isarTargetEcho = {isarEcho.perTarget.total};

beatRequest = struct( ...
    'fastTimeSec',           isarScene.fastTimeSec, ...
    'referenceStartTimeSec', isarScene.referenceStartTimeSec, ...
    'f0Hz',                  cfgISAR.fL, ...
    'f1Hz',                  cfgISAR.fH, ...
    'pulseWidthSec',         cfgISAR.Tr, ...
    'referenceDelaySec',     isarScene.tauRefSec, ...
    'outputSampleRateHz',    isarScene.fsScopeHz);

beatOptions = isarScene.dechirpOptions;
if nargin >= 3 && ~isempty(options)
    if ~(isstruct(options) && isscalar(options))
        error('RV_ISAR_Path:InvalidDechirpOptions', 'options 必须是 1×1 结构体。');
    end
    overrideNames = fieldnames(options);
    for k = 1:numel(overrideNames)
        beatOptions.(overrideNames{k}) = options.(overrideNames{k});
    end
end

[isarBeat, isarBeatInfo] = radar.lfm_local_reference_beat( ...
    isarTargetEcho, beatRequest, beatOptions);

if ~isfield(isarBeat, 'perTarget') || numel(isarBeat.perTarget) ~= nTarget
    error('RadarProject:ISARPerTargetBeatMissing', ...
        '逐目标 ISAR 后处理要求 isarBeat.perTarget 与 targets 数量一致。');
end
for targetIndex = 1:nTarget
    if ~isfield(isarBeat.perTarget(targetIndex), 'signal') || ...
            size(isarBeat.perTarget(targetIndex).signal, 1) ~= Na_ISAR
        error('RadarProject:ISARBeatPulseCountMismatch', ...
            '目标 %d 的逐目标拍频脉冲数与 cfg.Na=%d 不一致。', targetIndex, Na_ISAR);
    end
end
if abs(isarBeat.sampleRateHz - isarScene.fsScopeHz) > ...
        max(1, eps(isarScene.fsScopeHz) * max(1, isarScene.fsScopeHz))
    error('RadarProject:ISARBeatSampleRateMismatch', ...
        '拍频采样率 %.12g Hz 与配置的 fs_scope %.12g Hz 不一致。', ...
        isarBeat.sampleRateHz, isarScene.fsScopeHz);
end
end
