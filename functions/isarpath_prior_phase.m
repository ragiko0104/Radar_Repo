function reference = isarpath_prior_phase(cfg, prior, t, p, shifts, mode)
%ISARPATH_PRIOR_PHASE 在 ISAR 原始去斜/色散/截取处理规范下计算解析运动学先验相位。
%
% 功能说明：
%   用于 MODEL_CV / MODEL_CA 先验相位聚焦模式：不依赖真实散射点几何或微动真值，
%   仅依据外部给定的先验位置/速度/加速度和已知处理坐标，在每脉冲有效接收孔径中点
%   解析求取去斜正中频相位，并扣除色散补偿相位与随动 FFT 截取/RCMC 对齐引入的
%   坐标相位项，从而避免从合成 HRRP 零点取相时出现 π 跳变。
%
%   总先验相位（单位：周 cycles）由三部分叠加而成：
%     cycles(j) = positiveIfCycles(j) - dispersionCycles(j)
%                 - (alignedFftFreqHz(j) + fixedGateOffsetHz) * localWindowTimeSec(j)
%   最终返回相对第 1 脉冲归一化的单位模值复相位补偿因子：
%     phaseCorrection(j) = exp(-1j * 2*pi * (cycles(j) - cycles(1)))
%
% 输入参数：
%   cfg    - 1×1 struct，ISAR 波形配置，使用：
%              .fs, .fs_scope, .Tr, .fL, .kr, .propagationSpeed
%   prior  - 1×1 struct，单目标运动先验，使用：
%              .epochSec          1×1 double，先验参考时刻 (s)
%              .rangeM            1×1 double，参考时刻斜距 (m)
%              .velocityMps       1×1 double，参考时刻径向速度 (m/s)
%              .accelerationMps2  1×1 double，径向加速度 (m/s^2，仅 MODEL_CA 使用)
%              .loStartSec        1×1 double，局部本振起始时刻 (s)
%              .beatOffsetHz      1×1 double，拍频偏置 (Hz)
%   t      - Na×1 double，各脉冲绝对慢时间轴 (s)
%   p      - 1×1 struct，后处理参数结构体（来自 isarpath_settings / scene.isarSettings）
%   shifts - Na×1 double，各脉冲 RCMC 距离对齐移位量 (bin)
%   mode   - 1×1 string/char，先验运动学模型："MODEL_CV"（常速）或 "MODEL_CA"（匀加速）
%
% 输出参数：
%   reference - 1×1 struct，解析先验相位参考结果：
%     .phaseCorrection      Na×1 complex double，逐脉冲单位模值复相位补偿因子
%     .phaseCycles          Na×1 double，未相对首脉冲归一化的总相位（单位：周 cycles）
%     .fixedGateOffsetHz    1×1 double，首脉冲虚拟固定子 bin 距离门相对对齐频率的频偏 (Hz)
%     .apertureMidpointSec  Na×1 double，每脉冲有效重叠孔径中点的快时间坐标 (s)
%     .method               1×1 string，相位求值方法标识字符串
%     .accelerationUsedMps2 1×1 double，本次计算实际采用的径向加速度 (m/s^2)

%% 1. 解析先验加速度（MODEL_CV 强制为 0，MODEL_CA 使用 prior.accelerationMps2）
accelerationMps2 = 0;
if string(mode) == "MODEL_CA"
    accelerationMps2 = prior.accelerationMps2;
end

pulseOffsetSec = t(:) - prior.epochSec;
loStartQuantizedSec = round(prior.loStartSec * cfg.fs) / cfg.fs;

%% 2. 计算每脉冲时域截窗起始时刻（与 isarpath_extract_segment 完全一致）
segmentStartIndex = round( ...
    (prior.loStartSec - p.timeMarginSamples / cfg.fs_scope) * cfg.fs_scope) + ...
    ones(size(pulseOffsetSec));
if p.timeWindowMode == "MOVING_PULSE"
    segmentStartIndex = round( ...
        (prior.loStartSec + 2 * prior.velocityMps * pulseOffsetSec / cfg.propagationSpeed - ...
        p.timeMarginSamples / cfg.fs_scope) * cfg.fs_scope) + 1;
end
segmentStartSec = (segmentStartIndex - 1) / cfg.fs_scope;

%% 3. 不动点迭代求取回波门、本振门与截窗三者的有效重叠孔径中点
echoStartSec = loStartQuantizedSec + zeros(size(pulseOffsetSec));
echoStopSec  = echoStartSec + cfg.Tr;
for iteration = 1:4
    echoStartSec = 2 * range_at(pulseOffsetSec + echoStartSec) / cfg.propagationSpeed;
    echoStopSec  = cfg.Tr + 2 * range_at(pulseOffsetSec + echoStopSec) / cfg.propagationSpeed;
end

overlapStartSec = max(max(echoStartSec, loStartQuantizedSec), segmentStartSec);
overlapStopSec  = min(min(echoStopSec - 1 / cfg.fs, loStartQuantizedSec + cfg.Tr - 1 / cfg.fs), ...
    segmentStartSec + (p.timeWindowSamples - 1) / cfg.fs_scope);
assert(all(overlapStopSec > overlapStartSec), 'Prior-predicted echo has no LO overlap.');

apertureMidpointFastSec = (overlapStartSec + overlapStopSec) / 2;
localWindowTimeSec      = apertureMidpointFastSec - segmentStartSec;

%% 4. 在孔径中点计算去斜正中频相位与瞬时拍频
rangeAtMidpointM     = range_at(pulseOffsetSec + apertureMidpointFastSec);
velocityAtMidpointMps = prior.velocityMps + ...
    accelerationMps2 * (pulseOffsetSec + apertureMidpointFastSec);

roundTripOffsetSec = 2 * rangeAtMidpointM / cfg.propagationSpeed - loStartQuantizedSec;
loElapsedSec       = apertureMidpointFastSec - loStartQuantizedSec;
positiveIfCycles   = roundTripOffsetSec .* (cfg.fL + cfg.kr * loElapsedSec) - ...
    0.5 * cfg.kr * roundTripOffsetSec.^2;

%% 5. 扣除脉内多普勒色散补偿引入的相位与频率偏移
dispersionCycles = zeros(size(pulseOffsetSec));
dispersionHz     = zeros(size(pulseOffsetSec));
if p.dispersionCompensationEnabled
    dispersionScale  = 2 * prior.velocityMps / cfg.propagationSpeed;
    dispersionCycles = dispersionScale * (cfg.kr * localWindowTimeSec.^2 + cfg.fL * localWindowTimeSec);
    dispersionHz     = dispersionScale * (2 * cfg.kr * localWindowTimeSec + cfg.fL);
end

%% 6. 扣除随动截取中心 + RCMC 对齐对应的 FFT 坐标相位，生成单位模值复补偿因子
modelBeatFreqHz = prior.beatOffsetHz + ...
    2 * cfg.kr * prior.velocityMps * pulseOffsetSec / cfg.propagationSpeed;
cropCenterBins  = round(modelBeatFreqHz / p.frequencyBinSpacingHz) + 1 + ...
    p.cropCenterOffsetBinsEffective;
alignedFftFreqHz = (cropCenterBins - 1 + shifts(:)) * p.frequencyBinSpacingHz;

instantaneousBeatHz = cfg.kr * roundTripOffsetSec + ...
    2 * velocityAtMidpointMps / cfg.propagationSpeed .* ...
    (cfg.fL + cfg.kr * (apertureMidpointFastSec - 2 * rangeAtMidpointM / cfg.propagationSpeed)) - ...
    dispersionHz;

% 在首脉冲处固定一个虚拟子 bin 距离门偏置，并在整个 CPI 内保持恒定（对应真实 HRRP 同一列）
fixedGateOffsetHz = instantaneousBeatHz(1) - alignedFftFreqHz(1);
totalPhaseCycles  = positiveIfCycles - dispersionCycles - ...
    (alignedFftFreqHz + fixedGateOffsetHz) .* localWindowTimeSec;
phaseCorrection   = exp(-1j * 2 * pi * (totalPhaseCycles - totalPhaseCycles(1)));

reference = struct( ...
    'phaseCorrection',      phaseCorrection, ...
    'phaseCycles',          totalPhaseCycles, ...
    'fixedGateOffsetHz',    fixedGateOffsetHz, ...
    'apertureMidpointSec',  apertureMidpointFastSec, ...
    'method',               "prior polynomial + positive IF + dispersion + crop/alignment coordinates", ...
    'accelerationUsedMps2', accelerationMps2);

    function rangeM = range_at(deltaTimeSec)
        rangeM = prior.rangeM + prior.velocityMps * deltaTimeSec + ...
            0.5 * accelerationMps2 * deltaTimeSec.^2;
    end
end
