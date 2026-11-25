function [rvEcho, rvEchoBaseband] = rvpath_generate_echo(targets, radarState, wf, c, options)
%RVPATH_GENERATE_ECHO 生成通用散射体目标的六段 LFM 测距测速回波并下变频至基带。
%
% 功能说明：
%   1. 按六段 LFM 波形依次调用 radar.echo_generate_target_lfm（按单目标 AutoFastTimeGate）
%      生成多目标回波；
%   2. 测速链路刚需快时间多普勒效应，回波模式固定为 'fastTimeExact'（不提供快时间冻结）；
%   3. 默认关闭 ReturnGroups 与 chunk 分块以获得最高速度；如需查看 body/rotor 分组或
%      仿真 UAV 等数百散射点目标，可通过 wf 或 options 开启对应选项；
%   4. 对总射频回波乘以复数本振 exp(-1j*2*pi*fc*fastTimeSec) 得到基带回波 rvEchoBaseband。
%
% 输入参数：
%   targets    - Ntarget×1 struct，统一散射体目标数组（支持 UAV/legacyDisk/点云/单点混排）
%   radarState - 1×1 struct，雷达位置/速度/极化/多普勒符号状态结构体
%   wf         - 1×1 struct，由 rvpath_derive_waveform 产出的派生 RV 波形结构体
%                （使用 Na/N/fc/Tr/f0_v/f1_v/delta_tau/fastTimeSec/slowTimeSec，
%                可选 returnGroups/fastTimeChunkSize/scattererChunkSize）
%   c          - 1×1 double，光速 (m/s)
%   options    - (可选) 1×1 struct，可选覆盖项：
%     .targetAmplitude     Ntarget×1 或 Ntarget×Na double，目标外部幅度（默认全 1）
%     .returnGroups        1×1 logical，是否按散射体分组返回（默认读取 wf.returnGroups 或 false）
%     .fastTimeChunkSize   1×1 double/[]/Inf，快时间分块大小（默认 Inf 即不分块）
%     .scattererChunkSize  1×1 double/[]/Inf，散射点分块大小（默认 Inf 即不分块）
%
% 输出参数：
%   rvEcho         - 1×1 struct，射频域六段叠加回波容器：
%     .total       Na×N complex double，所有目标相干求和后的总回波
%     .groups      Ng×1 struct，returnGroups=true 时按分组保存的回波（关闭时为空 struct([])）
%     .groupKeys   Ng×1 string，各分组键名
%     .groupNames  Ng×1 string，各分组显示名称
%     .body/.rotor (可选) Na×N complex double，启用 returnGroups 且存在对应分组时提供
%   rvEchoBaseband - Na×N complex double，去除射频载频 fc 后的六段基带回波

if nargin < 5 || isempty(options)
    options = struct();
end

Na = wf.Na;
N  = wf.N;
nTargets = numel(targets);

if isfield(options, 'targetAmplitude') && ~isempty(options.targetAmplitude)
    targetAmp = double(options.targetAmplitude);
    if iscolumn(targetAmp) || isrow(targetAmp)
        targetAmp = targetAmp(:) .* ones(1, Na);
    end
else
    targetAmp = ones(nTargets, Na);
end

returnGroups = false;
if isfield(wf, 'returnGroups') && ~isempty(wf.returnGroups)
    returnGroups = logical(wf.returnGroups);
end
if isfield(options, 'returnGroups') && ~isempty(options.returnGroups)
    returnGroups = logical(options.returnGroups);
end

fastTimeChunkSize  = resolve_chunk_opt(wf, options, 'fastTimeChunkSize');
scattererChunkSize = resolve_chunk_opt(wf, options, 'scattererChunkSize');

rvEcho = struct( ...
    'total',      complex(zeros(Na, N)), ...
    'groups',     struct([]), ...
    'groupKeys',  strings(0, 1), ...
    'groupNames', strings(0, 1));

% 测速链路刚需脉内快时间相位演化，固定使用 fastTimeExact。
rvEchoOptions = struct( ...
    'Mode',               'fastTimeExact', ...
    'ReturnGroups',       returnGroups, ...
    'ReturnPerTarget',    false, ...
    'AutoFastTimeGate',   true, ...
    'FastTimeChunkSize',  fastTimeChunkSize, ...
    'ScattererChunkSize', scattererChunkSize);

for segmentIndex = 1:numel(wf.f0_v)
    radarStateSegment = radarState;
    segmentCarrierHz = wf.fc + 0.5 .* (wf.f0_v(segmentIndex) + wf.f1_v(segmentIndex));
    radarStateSegment.carrierFrequencyHz = segmentCarrierHz;
    radarStateSegment.wavelengthM        = c ./ segmentCarrierHz;

    echoRequestRV = struct( ...
        'pulseTimeSec',     wf.slowTimeSec, ...
        'fastTimeSec',      wf.fastTimeSec, ...
        'waveformType',     'LFM', ...
        'f0Hz',             wf.f0_v(segmentIndex), ...
        'f1Hz',             wf.f1_v(segmentIndex), ...
        'fLoHz',            wf.fc, ...
        'pulseWidthSec',    wf.Tr, ...
        'deltaTauSec',      wf.delta_tau(segmentIndex), ...
        'propagationSpeed', c, ...
        'pathLossMode',     'none', ...
        'targetAmplitude',  targetAmp);

    segmentEcho = radar.echo_generate_target_lfm( ...
        targets, radarStateSegment, echoRequestRV, rvEchoOptions);

    if returnGroups
        if isempty(rvEcho.groupKeys)
            rvEcho.groupKeys  = segmentEcho.groupKeys;
            rvEcho.groupNames = segmentEcho.groupNames;
            rvEcho.groups = repmat(struct( ...
                'key', "", 'name', "", 'signal', complex(zeros(Na, N))), ...
                numel(segmentEcho.groups), 1);
            for groupIndex = 1:numel(rvEcho.groups)
                rvEcho.groups(groupIndex).key  = segmentEcho.groups(groupIndex).key;
                rvEcho.groups(groupIndex).name = segmentEcho.groups(groupIndex).name;
            end
        elseif ~isequal(lower(rvEcho.groupKeys), lower(segmentEcho.groupKeys))
            error('multiuav:rvpath:groupMismatch', ...
                'Scatterer echo groups changed between LFM segments.');
        end
        for groupIndex = 1:numel(rvEcho.groups)
            rvEcho.groups(groupIndex).signal = rvEcho.groups(groupIndex).signal + ...
                segmentEcho.groups(groupIndex).signal;
        end
    else
        rvEcho.total = rvEcho.total + segmentEcho.total;
    end
end

if returnGroups
    for groupIndex = 1:numel(rvEcho.groups)
        rvEcho.total = rvEcho.total + rvEcho.groups(groupIndex).signal;
    end
    bodyIndex  = find(strcmpi(rvEcho.groupKeys, 'body'), 1);
    rotorIndex = find(strcmpi(rvEcho.groupKeys, 'rotor'), 1);
    if ~isempty(bodyIndex),  rvEcho.body  = rvEcho.groups(bodyIndex).signal;  end
    if ~isempty(rotorIndex), rvEcho.rotor = rvEcho.groups(rotorIndex).signal; end
end

rvCarrierRemoval = exp(-1j * 2 * pi * wf.fc .* wf.fastTimeSec);
rvEchoBaseband   = rvEcho.total .* rvCarrierRemoval;
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
