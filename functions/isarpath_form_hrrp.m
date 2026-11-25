function [hrrp, cropCentres, hrrpBeforeDispersion] = isarpath_form_hrrp( ...
    segmentComp, cfgISAR, prior, pulseTimeSec, p, segmentRaw)
%ISARPATH_FORM_HRRP 将时域截窗拍频经距离维 FFT 与先验随动截取形成 HRRP-慢时间矩阵。
%
% 功能说明：
%   1. 对色散补偿后的时域截窗矩阵 segmentComp 沿快时间维（第 2 维）执行
%      p.rangeFftPoints 点 FFT，并除以 p.timeWindowSamples 保持幅值尺度一致；
%   2. 按先验拍频偏置 beatOffsetHz 与径向速度 velocityMps 计算各脉冲的随动截取中心
%      bin（cropCentres），调用 lib 中的 signal.crop_fft_bins 截取 p.hrrfBins 点，
%      形成 Na×Nr 复一维距离像-慢时间矩阵 hrrp；
%   3. 若同时传入色散补偿前的 segmentRaw（且需要输出第 3 参数），以相同截取窗口同步生成
%      hrrpBeforeDispersion，用于量化色散补偿前后的 HRRP 3dB 展宽与 PSLR 改善。
%
% 输入参数：
%   segmentComp  - Na×Nwin complex double，色散补偿后的时域截窗矩阵
%   cfgISAR      - 1×1 struct，ISAR 波形配置（使用 kr、propagationSpeed）
%   prior        - 1×1 struct，单目标先验结构体（epochSec、velocityMps、beatOffsetHz）
%   pulseTimeSec - Na×1 double，ISAR 绝对慢时间轴 (s)
%   p            - 1×1 struct，由 isarpath_settings 生成的后处理参数结构体
%                  （使用 timeWindowSamples、rangeFftPoints、hrrfBins、
%                  frequencyBinSpacingHz、cropCenterOffsetBinsEffective）
%   segmentRaw   - (可选) Na×Nwin complex double，色散补偿前的截窗矩阵（不传时
%                  hrrpBeforeDispersion 返回空矩阵 []）
%
% 输出参数：
%   hrrp                 - Na×Nr complex double，色散补偿并随动截取后的 HRRP-慢时间矩阵
%                          （Nr = p.hrrfBins）
%   cropCentres          - Na×1 double，各脉冲在距离 FFT 谱中的 1-based 截取中心 bin 索引
%   hrrpBeforeDispersion - Na×Nr complex double 或 []，色散补偿前的同窗口 HRRP 矩阵

nPulses = size(segmentComp, 1);
needBefore = (nargout >= 3) && (nargin >= 6) && ~isempty(segmentRaw);

if needBefore
    beforeFft = fft(segmentRaw, p.rangeFftPoints, 2) / p.timeWindowSamples;
end

F = fft(segmentComp, p.rangeFftPoints, 2) / p.timeWindowSamples;
% beatOffsetHz > 0 表示本振在目标左侧（提前）；< 0 表示本振在目标右侧（滞后）。
% 经过 ANALYTIC_IF 实混频与正频解析支路保留后，拍频均折叠至正频轴 |fSigned|。
fSigned = prior.beatOffsetHz + ...
    2 * cfgISAR.kr * prior.velocityMps * (pulseTimeSec(:) - prior.epochSec) / ...
    cfgISAR.propagationSpeed;
beatSign = sign(prior.beatOffsetHz);
if beatSign == 0
    beatSign = 1;
end
cropCentres = round(abs(fSigned) / p.frequencyBinSpacingHz) + 1 + ...
    beatSign * p.cropCenterOffsetBinsEffective;

hrrp = complex(zeros(nPulses, p.hrrfBins, 'like', segmentComp));
if needBefore
    hrrpBeforeDispersion = complex(zeros(nPulses, p.hrrfBins, 'like', segmentRaw));
else
    hrrpBeforeDispersion = [];
end

for pulseIndex = 1:nPulses
    hrrp(pulseIndex, :) = signal.crop_fft_bins( ...
        F(pulseIndex, :), cropCentres(pulseIndex), p.hrrfBins);
    if needBefore
        hrrpBeforeDispersion(pulseIndex, :) = signal.crop_fft_bins( ...
            beforeFft(pulseIndex, :), cropCentres(pulseIndex), p.hrrfBins);
    end
end
end
