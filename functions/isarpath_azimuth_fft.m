function [isarImage, isValid] = isarpath_azimuth_fft(phaseCompensatedHrrf, p)
%ISARPATH_AZIMUTH_FFT 对相位补偿后的 HRRP-慢时间矩阵执行方位维加窗零填充 FFT 成像。
%
% 功能说明：
%   1. 当输入相位补偿 HRRP 存在非零能量时，调用 lib 中的 radar.ISAR_imaging
%      执行慢时间加窗（p.windowType）、慢时间零填充（倍率 p.slowTimeZeroPaddingFactor）
%      与方位维 FFT 并输出归一化幅度图像；
%   2. 若输入全零或结果含非有限值（NaN/Inf），安全返回零矩阵并置 isValid = false。
%
% 输入参数：
%   phaseCompensatedHrrf - Na×Nr（或 Na×Nr×Nc）complex double，相位补偿后的 HRRP 矩阵
%   p                    - 1×1 struct，由 isarpath_settings 生成的后处理参数结构体
%                          （使用 slowTimeZeroPaddingFactor、windowType）
%
% 输出参数：
%   isarImage - (Na*Kzp)×Nr double，归一化方位-距离二维 ISAR 幅度像
%               （Kzp = p.slowTimeZeroPaddingFactor，第 1 维为方位/多普勒，第 2 维为距离）
%   isValid   - 1×1 logical，图像包含有限正幅值时为 true，否则为 false

sceneHrrp  = phaseCompensatedHrrf(:, :, 1);
nCrossBins = size(sceneHrrp, 1) * p.slowTimeZeroPaddingFactor;
nRangeBins = size(sceneHrrp, 2);

if any(abs(sceneHrrp) > 0, 'all')
    isarImage = radar.ISAR_imaging( ...
        sceneHrrp, p.slowTimeZeroPaddingFactor, char(p.windowType));
else
    isarImage = zeros(nCrossBins, nRangeBins);
end

if any(~isfinite(isarImage(:)))
    warning('RadarProject:ISARImageNotFinite', ...
        'ISAR 图像包含 NaN/Inf，已将其替换为零。');
    isarImage = zeros(nCrossBins, nRangeBins);
end

isValid = any(isarImage(:) > 0);
end
