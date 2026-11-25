function wf = rvpath_derive_waveform(cfgRV)
%RVPATH_DERIVE_WAVEFORM 计算测距测速链路（RV_path）六段 LFM 派生波形参数与时间轴。
%
% 功能说明：
%   1. 校验 RV 链路基础配置（要求 Ta >= T）；
%   2. 根据三组带宽 B(1:3)、起始频偏 f1_up、频带间隔 f_space 及时间间隔 T_space，
%      派生六段上下扫频 LFM 的起止基带频率、等效中心射频频率、调频斜率与段起始时延；
%   3. 构造 RV 链路专用快时间轴 fastTimeSec（1×N）和慢时间轴 slowTimeSec（Na×1），
%      并将所有派生字段附加到原 cfgRV 结构体后输出。
%
% 输入参数：
%   cfgRV - 1×1 struct，RV 链路基础波形配置，必须包含：
%     .fs          1×1 double，基带采样率 (Hz)
%     .Tr          1×1 double，单段 LFM 脉宽 (s)
%     .T_space     1×1 double，相邻 LFM 段间保护时间间隔 (s)
%     .T           1×1 double，单帧复合波形快时间记录时长 (s)
%     .Ta          1×1 double，脉冲重复间隔 PRI (s)，要求 Ta >= T
%     .fc          1×1 double，雷达射频载频 (Hz)
%     .f_space     1×1 double，相邻 LFM 频带间隔 (Hz)
%     .B           1×3 或 3×1 double，三组频带带宽 [B1, B2, B3] (Hz)
%     .f1_up       1×1 double，频带 1 上扫起始基带频率 (Hz)
%     .Na          1×1 double，慢时间脉冲数
%
% 输出参数：
%   wf    - 1×1 struct，保留 cfgRV 全部输入字段，并新增以下派生字段：
%     .B1, .B2, .B3         1×1 double，三组带宽标量 (Hz)
%     .f0_v                 1×6 double，六段 LFM 基带起始频率 (Hz)，顺序 [1上,1下,2上,2下,3上,3下]
%     .f1_v                 1×6 double，六段 LFM 基带终止频率 (Hz)
%     .fc_v                 1×6 double，六段 LFM 等效射频中心频率 (Hz)
%     .k_v                  1×6 double，六段 LFM 调频斜率 (Hz/s，上扫为正、下扫为负)
%     .delta_tau            1×6 double，六段 LFM 在单帧内的起始时延偏移 (s)
%     .fmax_RV, .fmin_RV    1×1 double，六段基带频率最高/最低边界 (Hz)
%     .N                    1×1 double，单脉冲快时间采样点数 round(fs*T)
%     .N_lfm                1×1 double，单段 LFM 采样点数 Tr*fs
%     .N_exp2               1×1 double，nextpow2(N)
%     .N1, .N2              1×1 double，FFT 因子分解指数
%     .IF_FITACC            1×1 logical，慢时间脉冲数 Na > 1 时为 true
%     .fastTimeSec          1×N double，快时间轴 (s)
%     .slowTimeIntervalSec  1×1 double，慢时间间隔 (s)，等于 T
%     .slowTimeSec          Na×1 double，慢时间轴 (s)

if ~(isstruct(cfgRV) && isscalar(cfgRV))
    error('RV_ISAR_Path:InvalidRVConfig', 'cfgRV 必须是 1×1 结构体。');
end

T       = double(cfgRV.T);
fs      = double(cfgRV.fs);
Tr      = double(cfgRV.Tr);
T_space = double(cfgRV.T_space);
Ta      = double(cfgRV.Ta);
fc      = double(cfgRV.fc);
f_space = double(cfgRV.f_space);
f1_up   = double(cfgRV.f1_up);
B       = double(cfgRV.B(:).');
Na      = double(cfgRV.Na);

if Ta < T
    error('RV_ISAR_Path:InvalidRVTiming', ...
        'In RVpath, Ta (%.6g) should be at least the same with T (%.6g) or larger.', Ta, T);
end
if numel(B) ~= 3
    error('RV_ISAR_Path:InvalidRVBandwidth', 'cfgRV.B 必须包含 3 个频带带宽 [B1, B2, B3]。');
end

B1 = B(1);
B2 = B(2);
B3 = B(3);

f1_down = f1_up   + 2*B1 + f_space;
f2_up   = f1_down + f_space;
f2_down = f2_up   + 2*B2 + f_space;
f3_up   = f2_down + f_space;
f3_down = f3_up   + 2*B3 + f_space;

k1 = B1 / Tr;
k2 = B2 / Tr;
k3 = B3 / Tr;

fc1 = f1_up   + B1/2 + fc;
fc2 = f2_up   + B2/2 + fc;
fc3 = f3_up   + B3/2 + fc;
fc4 = f1_down - B1/2 + fc;
fc5 = f2_down - B2/2 + fc;
fc6 = f3_down - B3/2 + fc;

fc_v      = [fc1, fc4, fc2, fc5, fc3, fc6];
k_v       = [k1, -k1, k2, -k2, k3, -k3];
delta_tau = [0, Tr + T_space, Tr*2 + 2*T_space, Tr*2 + 3*T_space + Tr, ...
    Tr*4 + T_space*4, Tr*4 + T_space*5 + Tr];
f0_v      = [f1_up, f1_down, f2_up, f2_down, f3_up, f3_down];
f1_v      = [f1_up+B1, f1_down-B1, f2_up+B2, f2_down-B2, f3_up+B3, f3_down-B3];
fmax_RV   = max([f0_v, f1_v]);
fmin_RV   = min([f0_v, f1_v]);

N      = round(fs * T);
N_lfm  = Tr * fs;
N_exp2 = nextpow2(N);
N1     = round(N_exp2 / 2);
N2     = N_exp2 - N1;
IF_FITACC = (Na > 1);

fastTimeSec         = (0:N-1) ./ fs;
slowTimeIntervalSec = T;
slowTimeSec         = (0:Na-1).' .* slowTimeIntervalSec;

wf = cfgRV;
wf.B1                  = B1;
wf.B2                  = B2;
wf.B3                  = B3;
wf.fc_v                = fc_v;
wf.k_v                 = k_v;
wf.delta_tau           = delta_tau;
wf.f0_v                = f0_v;
wf.f1_v                = f1_v;
wf.fmax_RV             = fmax_RV;
wf.fmin_RV             = fmin_RV;
wf.N                   = N;
wf.N_lfm               = N_lfm;
wf.N_exp2              = N_exp2;
wf.N1                  = N1;
wf.N2                  = N2;
wf.IF_FITACC           = IF_FITACC;
wf.fastTimeSec         = fastTimeSec;
wf.slowTimeSec         = slowTimeSec;
wf.slowTimeIntervalSec = slowTimeIntervalSec;
end
