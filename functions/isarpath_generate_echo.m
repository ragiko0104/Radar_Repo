function [isarEcho, isarEchoInfo] = isarpath_generate_echo(isarScene, options)
%ISARPATH_GENERATE_ECHO 生成 ISAR 高采样率宽带复 LFM 模拟回波。
%
% 功能说明：
%   基于 isarpath_build_scene 构建的微动目标数组 targetsISAR 与 ISAR 载频雷达状态，
%   调用 lib 中的 radar.echo_generate_target_lfm（按单目标 AutoFastTimeGate）生成宽带复 LFM 回波。
%   支持 fastTimeExact（脉内精确演化）和 stopAndGo（快时间冻结）两种回波模式；
%   默认关闭 ReturnGroups 与 chunk 分块以大幅降低峰值内存与回调开销（仿 UAV 等数百散射点时
%   可通过 cfgISAR 或 options 开启分块）。
%
% 输入参数：
%   isarScene - 1×1 struct，由 isarpath_build_scene 产出的 ISAR 场景结构体
%               （使用 cfg/targetsISAR/radarState/pulseTimeSec/fastTimeSec/tauTxSec）
%   options   - (可选) 1×1 struct，可选覆盖项：
%     .targetAmplitude     N×1 或 N×Na double，逐目标/逐脉冲外部幅度（默认全 1）
%     .echoMode            char/string，覆盖 cfg.echoMode（"fastTimeExact"|"stopAndGo"）
%     .returnGroups        1×1 logical，是否按散射体分组返回（默认读取 cfg.returnGroups 或 false）
%     .fastTimeChunkSize   1×1 double/[]/Inf，快时间分块点数（默认 Inf 即不分块）
%     .scattererChunkSize  1×1 double/[]/Inf，散射点分块点数（默认 Inf 即不分块）
%
% 输出参数：
%   isarEcho     - 1×1 struct，宽带高采样率复回波容器：
%     .total       Na×NHigh complex double，所有目标相干叠加后的场景总回波
%     .groups      Ng×1 struct，returnGroups=true 时保存分组回波（关闭时为空 struct([])）
%     .perTarget   N×1 struct，逐目标独立回波（每项含 .total Na×NHigh）
%   isarEchoInfo - 1×1 struct，radar.echo_generate_target_lfm 返回的回波元信息

if nargin < 2 || isempty(options)
    options = struct();
end

cfgISAR = isarScene.cfg;
nTarget = isarScene.nTargets;
Na_ISAR = cfgISAR.Na;

if isfield(options, 'targetAmplitude') && ~isempty(options.targetAmplitude)
    targetAmp = double(options.targetAmplitude);
    if isvector(targetAmp)
        targetAmp = targetAmp(:) .* ones(1, Na_ISAR);
    end
else
    targetAmp = ones(nTarget, Na_ISAR);
end

echoMode = cfgISAR.echoMode;
if isfield(options, 'echoMode') && ~isempty(options.echoMode)
    echoMode = options.echoMode;
end

returnGroups = false;
if isfield(cfgISAR, 'returnGroups') && ~isempty(cfgISAR.returnGroups)
    returnGroups = logical(cfgISAR.returnGroups);
end
if isfield(options, 'returnGroups') && ~isempty(options.returnGroups)
    returnGroups = logical(options.returnGroups);
end

fastTimeChunkSize  = resolve_chunk_opt(cfgISAR, options, 'fastTimeChunkSize');
scattererChunkSize = resolve_chunk_opt(cfgISAR, options, 'scattererChunkSize');

isarEchoOptions = struct( ...
    'Mode',               echoMode, ...
    'ReturnGroups',       returnGroups, ...
    'ReturnPerTarget',    true, ...
    'AutoFastTimeGate',   true, ...
    'FastTimeChunkSize',  fastTimeChunkSize, ...
    'ScattererChunkSize', scattererChunkSize);

isarEchoRequest = struct( ...
    'pulseTimeSec',     isarScene.pulseTimeSec, ...
    'fastTimeSec',      isarScene.fastTimeSec, ...
    'f0Hz',             cfgISAR.fL, ...
    'f1Hz',             cfgISAR.fH, ...
    'fLoHz',            0, ...
    'pulseWidthSec',    cfgISAR.Tr, ...
    'deltaTauSec',      isarScene.tauTxSec, ...
    'propagationSpeed', cfgISAR.propagationSpeed, ...
    'waveformType',     'LFM', ...
    'pathLossMode',     'none', ...
    'targetAmplitude',  targetAmp);

[isarEcho, isarEchoInfo] = radar.echo_generate_target_lfm( ...
    isarScene.targetsISAR, isarScene.radarState, isarEchoRequest, isarEchoOptions);
end

function chunkVal = resolve_chunk_opt(cfgStruct, optStruct, fieldName)
chunkVal = Inf; % 默认关闭分块
if isfield(cfgStruct, fieldName) && ~isempty(cfgStruct.(fieldName))
    chunkVal = double(cfgStruct.(fieldName));
end
if isfield(optStruct, fieldName) && ~isempty(optStruct.(fieldName))
    chunkVal = double(optStruct.(fieldName));
end
if isscalar(chunkVal) && chunkVal == 0
    chunkVal = Inf;
end
end
