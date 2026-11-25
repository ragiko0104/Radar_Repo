function [aligned, s, scores] = isarpath_align_profiles(H, p, forcedShift)
%ISARPATH_ALIGN_PROFILES 对单目标 HRRP-慢时间矩阵执行包络距离徙动校正（RCMC）。
%
% 功能说明：
%   按 p.rangeMigrationMode 与 p.rangeAlignmentMode 完成脉冲间距离像包络对齐：
%     - "FFT2048_SINGLE"      : 以第 1 脉冲为参考做 2*Nr 点互相关搜索（±rcmcSearchRadiusBins），
%                               FRACTIONAL 模式下叠加 3 点抛物线子像素峰拟合；
%     - "XCORR"/"MODEL_XCORR" : 调用 lib 中的 radar.range_migration_alignment；
%     - "NONE"                : 不施加平移（移位量全 0）。
%
% 输入参数：
%   H           - Na×Nr（或 Na×Nr×Nc）complex double，输入 HRRP-慢时间矩阵
%   p           - 1×1 struct，由 isarpath_settings / scene.isarSettings 生成的后处理参数结构体
%                 （使用 rangeMigrationMode、rangeAlignmentMode、
%                 rcmcCorrelationPoints、rcmcSearchRadiusBins）
%   forcedShift - (可选) Na×1 double，外部指定的各脉冲强制距离移位量 (bin)；
%                 省略或为空时由互相关自动估计
%
% 输出参数：
%   aligned - Na×Nr（与 H 同形）complex double，距离对齐后的 HRRP 矩阵
%   s       - Na×1 double，各脉冲相对参考脉冲的距离徙动补偿量 (bin)
%   scores  - Na×(2*R+1) double 或 []，FFT2048_SINGLE 模式下的互相关得分矩阵
%             （R = p.rcmcSearchRadiusBins；其余模式返回 []）

np = size(H, 1);
nr = size(H, 2);
scores = [];

if nargin >= 3 && ~isempty(forcedShift)
    s = forcedShift(:);
elseif p.rangeMigrationMode == "FFT2048_SINGLE"
    assert(p.rcmcCorrelationPoints == 2 * nr);
    candidates = -p.rcmcSearchRadiusBins:p.rcmcSearchRadiusBins;
    indices    = mod(candidates, p.rcmcCorrelationPoints) + 1;
    spectrum   = fft([abs(double(H(:, :, 1))), zeros(np, nr)], p.rcmcCorrelationPoints, 2);
    corr       = real(ifft(spectrum .* conj(spectrum(1, :)), p.rcmcCorrelationPoints, 2));
    scores     = corr(:, indices);
    s          = zeros(np, 1);
    for j = 2:np
        [~, pk] = max(scores(j, :));
        s(j) = candidates(pk);
        if p.rangeAlignmentMode == "FRACTIONAL" && pk > 1 && pk < numel(candidates)
            yy  = scores(j, pk-1:pk+1);
            den = yy(1) - 2 * yy(2) + yy(3);
            if den < 0
                s(j) = s(j) + max(-0.5, min(0.5, 0.5 * (yy(1) - yy(3)) / den));
            end
        end
    end
elseif ismember(p.rangeMigrationMode, ["XCORR", "MODEL_XCORR"])
    assert(ismatrix(H));
    wave = nr * ifft(ifftshift(H, 2), nr, 2);
    [~, aligned, s] = radar.range_migration_alignment( ...
        wave, H, 'freq', 1:nr, p.rangeAlignmentMode == "FRACTIONAL");
    return;
else
    s = zeros(np, 1);
end

if p.rangeAlignmentMode == "FRACTIONAL"
    aligned = isar_physical_fractional_shift(H, s);
else
    aligned = H;
    for channel = 1:size(H, 3)
        for j = 1:np
            aligned(j, :, channel) = H(j, mod((0:nr-1) + round(s(j)), nr) + 1, channel);
        end
    end
end
end
