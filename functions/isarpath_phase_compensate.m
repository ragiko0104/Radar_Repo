function [phaseCompensatedHrrf, phaseCorrection, phaseInfo] = isarpath_phase_compensate( ...
    alignedHrrp, p, cfgOrReference, prior, pulseTimeSec, rangeMigrationBins)
%ISARPATH_PHASE_COMPENSATE 对包络对齐后的 HRRP 执行脉冲间相位聚焦（特显点/先验模型补偿）。
%
% 功能说明：
%   按 p.phaseCompensationMode 对已对齐的 HRRP-慢时间矩阵施加逐脉冲复相位补偿：
%     - "VARIANCE_UNIT" : 调用 lib 中的 radar.phase_cmp（归一化幅度方差最小特显点法）
%                         估计逐脉冲复相位因子；若特显点无效则安全回退为无补偿并标记
%                         phaseInfo.mode = "NONE_INVALID_REFERENCE"；
%     - "MODEL_CV" / "MODEL_CA" : 使用预计算的先验相位参考结构体（3 参数调用）或调用
%                         isarpath_prior_phase（6 参数调用），基于先验运动学轨迹施加
%                         精确处理规范下的脉冲间复相位校正因子；
%     - "NONE"          : 不施加相位校正（校正因子全 1）。
%
% 输入参数：
%   alignedHrrp        - Na×Nr（或 Na×Nr×Nc）complex double，RCMC 对齐后的 HRRP 矩阵
%   p                  - 1×1 struct，后处理参数结构体（phaseCompensationMode、phaseThreshold）
%   cfgOrReference     - (可选) 在 MODEL_CV / MODEL_CA 模式下支持两种传参形式：
%                          (a) 3 参数形式：直接传入由 isarpath_prior_phase 产出的 reference 结构体；
%                          (b) 6 参数形式：传入 ISAR 波形配置 cfgISAR（同时提供后 3 个参数）。
%   prior              - (可选，6 参数形式必填) 1×1 struct，单目标先验结构体
%   pulseTimeSec       - (可选，6 参数形式必填) Na×1 double，绝对慢时间轴 (s)
%   rangeMigrationBins - (可选，6 参数形式必填) Na×1 double，RCMC 距离移位量 (bin)
%
% 输出参数：
%   phaseCompensatedHrrf - Na×Nr（与 alignedHrrp 同形）complex double，相位补偿后 HRRP
%   phaseCorrection      - Na×1 complex double，各脉冲施加的单位模值复相位校正因子
%   phaseInfo            - 1×1 struct，相位补偿诊断信息（含 mode、referenceBin 等）

sceneHrrp = alignedHrrp(:, :, 1);
nPulses   = size(sceneHrrp, 1);
phaseMode = upper(string(p.phaseCompensationMode));

phaseInfo = struct( ...
    'mode',                           phaseMode, ...
    'referenceBin',                   NaN, ...
    'minimumReferenceAmplitudeRatio', NaN);
appliedCorrection = ones(nPulses, 1);

if phaseMode == "VARIANCE_UNIT"
    changed = radar.phase_cmp(sceneHrrp, p.phaseThreshold);
    if all(isfinite(changed(:)))
        rowPower = sum(abs(sceneHrrp).^2, 2);
        goodRows = rowPower > 0;
        appliedCorrection(goodRows) = sum( ...
            changed(goodRows, :) .* conj(sceneHrrp(goodRows, :)), 2) ./ rowPower(goodRows);
    else
        phaseInfo.mode = "NONE_INVALID_REFERENCE";
    end
elseif ismember(phaseMode, ["MODEL_CV", "MODEL_CA"])
    if nargin >= 3 && isstruct(cfgOrReference) && isfield(cfgOrReference, 'phaseCorrection')
        priorRef = cfgOrReference;
    elseif nargin >= 6 && ~isempty(cfgOrReference) && ~isempty(prior) && ~isempty(pulseTimeSec)
        priorRef = isarpath_prior_phase( ...
            cfgOrReference, prior, pulseTimeSec, p, rangeMigrationBins, phaseMode);
    else
        error('RV_ISAR_Path:MissingPriorPhaseInputs', ...
            'MODEL_CV / MODEL_CA 模式必须提供 priorReference 结构体，或完整提供 cfg/prior/pulseTime/shifts。');
    end
    appliedCorrection = priorRef.phaseCorrection;
    assert(isequal(size(appliedCorrection), [nPulses, 1]));
    assert(max(abs(abs(appliedCorrection) - 1)) < 1e-12);
    phaseInfo.referenceMethod      = priorRef.method;
    phaseInfo.accelerationUsedMps2 = priorRef.accelerationUsedMps2;
elseif phaseMode ~= "NONE"
    error('RV_ISAR_Path:UnsupportedPhaseMode', ...
        '不支持的 phaseCompensationMode="%s"。', phaseMode);
end

phaseCompensatedHrrf = alignedHrrp .* appliedCorrection;
phaseCorrection      = appliedCorrection;
end
