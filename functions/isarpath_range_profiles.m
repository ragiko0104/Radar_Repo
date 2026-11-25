function [H, before, centres] = isarpath_range_profiles(segment, cfg, prior, pulseTime, p)
%ISARPATH_RANGE_PROFILES 兼容聚合接口：对时域截窗拍频执行色散补偿、距离维 FFT 与随动截取。
%
% 功能说明：
%   本函数作为共享/兼容聚合入口（供 MD_ZT 等外部链路直接调用），内部统一转发至
%   radar.doppler_dispersion_cmp 与 isarpath_form_hrrp，避免维护两套重复实现。
%
% 输入参数：
%   segment   - Na×Nwin complex double，由 isarpath_extract_segment 截取的时域拍频矩阵
%   cfg       - 1×1 struct，ISAR 波形配置（使用 fs_scope、fL、kr、propagationSpeed）
%   prior     - 1×1 struct，单目标先验结构体（epochSec、velocityMps、beatOffsetHz）
%   pulseTime - Na×1 double，ISAR 绝对慢时间轴 (s)
%   p         - 1×1 struct，由 isarpath_settings 生成的后处理参数结构体
%
% 输出参数：
%   H       - Na×Nr complex double，色散补偿并随动截取后的 HRRP-慢时间矩阵
%   before  - Na×Nr complex double 或 []，p.qualityEnabled==true 时为色散补偿前 HRRP，否则为 []
%   centres - Na×1 double，各脉冲在距离 FFT 谱中的 1-based 截取中心 bin 索引

if p.dispersionCompensationEnabled
    segmentComp = radar.doppler_dispersion_cmp( ...
        segment, cfg.fs_scope, cfg.fL, cfg.kr, ...
        prior.velocityMps * ones(1, size(segment, 1)), 'freq');
else
    segmentComp = segment;
end

if p.qualityEnabled
    [H, centres, before] = isarpath_form_hrrp( ...
        segmentComp, cfg, prior, pulseTime, p, segment);
else
    [H, centres] = isarpath_form_hrrp( ...
        segmentComp, cfg, prior, pulseTime, p);
    before = [];
end
end
