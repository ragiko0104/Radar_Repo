function [segmentComp, segmentRaw, segmentStarts] = isarpath_window_and_dispersion( ...
    beatSignal, cfgISAR, prior, pulseTimeSec, p)
%ISARPATH_WINDOW_AND_DISPERSION 对单目标窄带拍频执行时域截窗与多普勒色散补偿。
%
% 功能说明：
%   1. 调用已有共享算子 isarpath_extract_segment，按固定首脉冲窗（FIXED_FIRST_PULSE）
%      或随动窗（MOVING_PULSE）从窄带拍频中截取长度为 timeWindowSamples 的快时间片段；
%   2. 当 p.dispersionCompensationEnabled == true 时，调用 lib 中的
%      radar.doppler_dispersion_cmp（'freq' 模式）按先验速度补偿脉内多普勒展宽色散；
%   3. 同时保留补偿前截窗矩阵 segmentRaw，便于后续对比色散补偿前后的 HRRP 展宽指标。
%
% 输入参数：
%   beatSignal   - Na×NBeat complex double，单个目标的窄带去斜拍频信号
%   cfgISAR      - 1×1 struct，ISAR 波形配置（使用 fs_scope、fL、kr、propagationSpeed）
%   prior        - 1×1 struct，单目标先验结构体（epochSec、velocityMps、loStartSec）
%   pulseTimeSec - Na×1 double，ISAR 绝对慢时间轴 (s)
%   p            - 1×1 struct，由 isarpath_settings 生成的后处理参数结构体
%                  （使用 timeWindowSamples、timeMarginSamples、timeWindowMode、
%                  dispersionCompensationEnabled）
%
% 输出参数：
%   segmentComp   - Na×Nwin complex double，完成多普勒色散补偿后的时域截窗矩阵
%                   （Nwin = p.timeWindowSamples；未开启色散补偿时等于 segmentRaw）
%   segmentRaw    - Na×Nwin complex double，色散补偿前的原始时域截窗矩阵
%   segmentStarts - Na×1 double，每脉冲时域截窗在 beatSignal 中的 1-based 起始样本索引

[segmentRaw, segmentStarts] = isarpath_extract_segment( ...
    beatSignal, cfgISAR, prior, pulseTimeSec, p);

if p.dispersionCompensationEnabled
    segmentComp = radar.doppler_dispersion_cmp( ...
        segmentRaw, cfgISAR.fs_scope, cfgISAR.fL, cfgISAR.kr, ...
        prior.velocityMps * ones(1, size(segmentRaw, 1)), 'freq');
else
    segmentComp = segmentRaw;
end
end
