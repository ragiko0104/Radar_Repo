function [corrected, image, correction, info] = isarpath_phase_image(H, p, reference)
%ISARPATH_PHASE_IMAGE 兼容聚合接口：顺序执行脉冲间相位补偿与方位维 FFT 成像。
%
% 功能说明：
%   本函数作为共享/兼容聚合入口（供 MD_ZT 等外部链路直接调用），内部直接转发至
%   isarpath_phase_compensate 与 isarpath_azimuth_fft，消除重复实现。
%
% 输入参数：
%   H         - Na×Nr（或 Na×Nr×Nc）complex double，RCMC 对齐后的 HRRP 矩阵
%   p         - 1×1 struct，由 isarpath_settings 生成的后处理参数结构体
%   reference - (可选) 1×1 struct，由 isarpath_prior_phase 产出的先验相位参考结构体
%               （MODEL_CV / MODEL_CA 模式必填）
%
% 输出参数：
%   corrected  - Na×Nr（与 H 同形）complex double，相位补偿后 HRRP
%   image      - (Na*Kzp)×Nr double，归一化方位-距离二维 ISAR 幅度像
%   correction - Na×1 complex double，逐脉冲复相位校正因子
%   info       - 1×1 struct，相位补偿诊断信息

if nargin < 3
    reference = [];
end

[corrected, correction, info] = isarpath_phase_compensate(H, p, reference);
image = isarpath_azimuth_fft(corrected, p);
end
