function rvPulseCompression = rvpath_pulse_compress(rvEchoBaseband, wf, spectrumMode)
%RVPATH_PULSE_COMPRESS 生成六段 LFM 本振参考并执行频域脉冲压缩。
%
% 功能说明：
%   支持三种本振参考生成模式，均通过 lib 中的 signal.MF_Freq 在频域完成六段 LFM
%   （LFM11/12/21/22/31/32）脉冲压缩并截取前 N 点作为有效脉压复信号：
%     1. 'floating'  （默认）：调用 radar.lfm_online_spectrum 在线生成浮点频域本振，
%                              经 IFFT 转至时域后送入 signal.MF_Freq；
%     2. 'fixedPoint'        ：调用 radar.lfm_online_spectrum 在线生成 FPGA 定点量化
%                              频域本振（含 192-bit 配置字），经 IFFT 后送入 signal.MF_Freq；
%     3. 'idealTime'         ：调用 lib 中的 signal.LFM_wave 生成理想时域 LFM 本振，
%                              经 signal.keep_positive_frequency 取正频解析支路并由
%                              signal.sig_Delay_and_PhaseShift 施加段起始时延 delta_tau，
%                              直接送入 signal.MF_Freq 完成脉冲压缩（支持别名
%                              'timeDomain' / 'ideal'）。
%
% 输入参数：
%   rvEchoBaseband - Na×N complex double，去载频后的六段复合基带回波
%   wf             - 1×1 struct，由 rvpath_derive_waveform 产出的派生 RV 波形结构体
%                    （使用 fs、Tr、N、N_lfm、f0_v、f1_v、delta_tau、fastTimeSec）
%   spectrumMode   - (可选) char 或 1×1 string，本振参考生成模式：
%                    'floating'（默认） | 'fixedPoint' | 'idealTime'
%                    （'timeDomain' / 'ideal' 自动归一为 'idealTime'）
%
% 输出参数：
%   rvPulseCompression - 1×1 struct，六段频域脉压结果：
%     .segmentKeys                6×1 string，段键名 ["lfm11";...;"lfm32"]
%     .segmentNames               6×1 string，段中文名 ["频带1上扫";...;"频带3下扫"]
%     .fftLength                  1×1 double，FFT 长度 2^nextpow2(N)
%     .sampleCount                1×1 double，有效快时间点数 N
%     .fastTimeSec                1×N double，快时间轴 (s)
%     .localOscillatorGenerator   char，参考生成器标识：
%                                 'radar.lfm_online_spectrum' 或 'signal.LFM_wave'
%     .localOscillatorMode        char，生效模式 'floating' | 'fixedPoint' | 'idealTime'
%     .frequencyDomainMethod      char，匹配滤波实现名 'signal.MF_Freq'
%     .localOscillatorTime        6×1 cell，每元素为 1×fftLength complex 时域本振
%     .localOscillatorSpectrum    6×1 cell，每元素为 1×fftLength complex 频域本振
%     .localOscillatorInfo        6×1 struct，每段本振生成元信息
%     .localOscillatorConfigData  6×1 cell，每段 FPGA 192-bit 配置字（非定点模式为空字符串）
%     .compressed                 6×1 cell，每元素为 Na×N complex 脉压复信号

if nargin < 3 || isempty(spectrumMode)
    spectrumMode = 'floating';
end
modeNormalized = normalize_lo_mode(spectrumMode);

N      = wf.N;
nSeg   = numel(wf.f0_v);
fftLen = 2^nextpow2(N);

if strcmp(modeNormalized, 'idealTime')
    loGeneratorName = 'signal.LFM_wave';
else
    loGeneratorName = 'radar.lfm_online_spectrum';
end

rvPulseCompression = struct();
rvPulseCompression.segmentKeys = ...
    ["lfm11"; "lfm12"; "lfm21"; "lfm22"; "lfm31"; "lfm32"];
rvPulseCompression.segmentNames = ...
    ["频带1上扫"; "频带1下扫"; "频带2上扫"; "频带2下扫"; "频带3上扫"; "频带3下扫"];
rvPulseCompression.fftLength                 = fftLen;
rvPulseCompression.sampleCount               = N;
rvPulseCompression.fastTimeSec               = wf.fastTimeSec;
rvPulseCompression.localOscillatorGenerator  = loGeneratorName;
rvPulseCompression.localOscillatorMode       = modeNormalized;
rvPulseCompression.frequencyDomainMethod     = 'signal.MF_Freq';
rvPulseCompression.localOscillatorTime       = cell(nSeg, 1);
rvPulseCompression.localOscillatorSpectrum   = cell(nSeg, 1);
rvPulseCompression.localOscillatorInfo       = struct([]);
rvPulseCompression.localOscillatorConfigData = cell(nSeg, 1);
rvPulseCompression.compressed                = cell(nSeg, 1);

for segIdx = 1:nSeg
    if strcmp(modeNormalized, 'idealTime')
        pulseReal = signal.LFM_wave( ...
            wf.f0_v(segIdx), wf.f1_v(segIdx), wf.Tr, 0, wf.fs).';
        loPadded = [pulseReal, zeros(1, fftLen - numel(pulseReal))];
        loAnalytic = signal.keep_positive_frequency(loPadded);
        loTime = signal.sig_Delay_and_PhaseShift( ...
            loAnalytic, wf.fs, wf.delta_tau(segIdx), 0, 'BB');
        loSpectrum   = fft(loTime);
        loConfigData = '';
        loInfo = struct( ...
            'mode',                 'idealTime', ...
            'sampleRateHz',         wf.fs, ...
            'startFrequencyHz',     wf.f0_v(segIdx), ...
            'endFrequencyHz',       wf.f1_v(segIdx), ...
            'fftLength',            fftLen, ...
            'pulseSampleCount',     numel(pulseReal), ...
            'pulseDurationSec',     wf.Tr, ...
            'transformDurationSec', fftLen ./ wf.fs, ...
            'chirpSlopeHzPerSec',   (wf.f1_v(segIdx) - wf.f0_v(segIdx)) ./ wf.Tr, ...
            'timeBiasSec',          wf.delta_tau(segIdx), ...
            'configDataHex',        loConfigData);
    else
        [loSpectrum, loInfo, loConfigData] = radar.lfm_online_spectrum( ...
            wf.fs, wf.f0_v(segIdx), wf.f1_v(segIdx), fftLen, ...
            round(wf.N_lfm), wf.delta_tau(segIdx), struct('Mode', modeNormalized));
        loTime = ifft(loSpectrum);
    end

    compressedFull = signal.MF_Freq(rvEchoBaseband, loTime);

    rvPulseCompression.localOscillatorTime{segIdx}       = loTime;
    rvPulseCompression.localOscillatorSpectrum{segIdx}   = loSpectrum;
    if segIdx == 1
        rvPulseCompression.localOscillatorInfo = repmat(loInfo, nSeg, 1);
    else
        rvPulseCompression.localOscillatorInfo(segIdx) = loInfo;
    end
    rvPulseCompression.localOscillatorConfigData{segIdx} = loConfigData;
    rvPulseCompression.compressed{segIdx}                = compressedFull(:, 1:N);
end
end

function mode = normalize_lo_mode(rawMode)
if ~(ischar(rawMode) && (isrow(rawMode) || isempty(rawMode))) && ...
        ~(isstring(rawMode) && isscalar(rawMode))
    error('RV_ISAR_Path:InvalidLfmSpectrumMode', ...
        'spectrumMode 必须为 ''floating''、''fixedPoint'' 或 ''idealTime''。');
end
switch lower(strtrim(char(string(rawMode))))
    case {'floating', 'float', 'double'}
        mode = 'floating';
    case {'fixedpoint', 'fixed-point', 'fixed', 'bitaccurate', 'bit-accurate'}
        mode = 'fixedPoint';
    case {'idealtime', 'ideal-time', 'ideal_time', 'timedomain', 'time-domain', 'time_domain', 'ideal'}
        mode = 'idealTime';
    otherwise
        error('RV_ISAR_Path:InvalidLfmSpectrumMode', ...
            '不支持的本振模式 "%s"，请使用 ''floating''、''fixedPoint'' 或 ''idealTime''。', ...
            char(string(rawMode)));
end
end
