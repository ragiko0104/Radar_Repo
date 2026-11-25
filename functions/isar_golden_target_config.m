function cfg = isar_golden_target_config(configPath)
%ISAR_GOLDEN_TARGET_CONFIG Per-target parameters from the approved golden.
%   CFG = ISAR_GOLDEN_TARGET_CONFIG(CONFIGPATH) reads the approved golden
%   config.json and returns the ten-target table together with the waveform
%   constants that the processing has to share with it.
%
%   Target identifiers in this repository are ZERO-BASED (the golden writes
%   them as target00..target09); MATLAB indices below are one-based, so
%   cfg.targetId = (1:10) - 1.
%
%   Fields:
%     targetId          0..9, as recorded by the golden
%     rangeM            centre range per target
%     velocityMps       radial (line-of-sight) velocity per target
%     accelerationMps2  radial acceleration per target
%     beatOffsetHz      per-target local-oscillator offset, alternates
%                       20 MHz / 45 MHz and must be read, never assumed
%     loStartSec        2*range/c - beatOffset/chirpRate, same offset
%     carrierFrequencyHz, lfmStartHz, lfmStopHz, bandwidthHz, chirpRateHzPerS
%     adcSampleRateHz, pulsePeriodS, pulseWidthS, pulseCount
%     sourcePath
%
%   The beat offset and the LO start time are derived from the SAME recorded
%   offset; a caller must not mix one target's offset with another's.

if nargin < 1 || isempty(configPath)
    root = fileparts(fileparts(mfilename('fullpath')));
    configPath = fullfile(root, 'config', 'golden_10targets_config.json');
end
if ~isfile(configPath)
    error('ISAR:Golden:MissingConfig', 'Golden config not found: %s', configPath);
end
j = jsondecode(fileread(configPath));
required = {'range_m', 'radial_velocity_mps', 'radial_acceleration_mps2', ...
    'frequency_bias_hz', 'carrier_frequency_hz', 'lfm_start_hz', ...
    'lfm_stop_hz', 'bandwidth_hz', 'chirp_rate_hz_per_s', ...
    'adc_sample_rate_hz', 'pulse_period_s', 'pulse_width_s', 'pulse_count'};
for k = 1:numel(required)
    if ~isfield(j, required{k})
        error('ISAR:Golden:MissingField', ...
            'Golden config lacks field %s: %s', required{k}, configPath);
    end
end
n = numel(j.range_m);
if numel(j.frequency_bias_hz) ~= n || numel(j.radial_velocity_mps) ~= n ...
        || numel(j.radial_acceleration_mps2) ~= n
    error('ISAR:Golden:LengthMismatch', ...
        'Golden per-target arrays disagree in length (%d targets).', n);
end
cfg = struct();
cfg.targetId = (1:n) - 1;
cfg.rangeM = double(j.range_m(:).');
cfg.velocityMps = double(j.radial_velocity_mps(:).');
cfg.accelerationMps2 = double(j.radial_acceleration_mps2(:).');
cfg.beatOffsetHz = double(j.frequency_bias_hz(:).');
cfg.carrierFrequencyHz = double(j.carrier_frequency_hz);
cfg.lfmStartHz = double(j.lfm_start_hz);
cfg.lfmStopHz = double(j.lfm_stop_hz);
cfg.bandwidthHz = double(j.bandwidth_hz);
cfg.chirpRateHzPerS = double(j.chirp_rate_hz_per_s);
cfg.adcSampleRateHz = double(j.adc_sample_rate_hz);
cfg.pulsePeriodS = double(j.pulse_period_s);
cfg.pulseWidthS = double(j.pulse_width_s);
cfg.pulseCount = double(j.pulse_count);
cfg.propagationSpeed = 3e8;
cfg.loStartSec = 2 * cfg.rangeM / cfg.propagationSpeed ...
    - cfg.beatOffsetHz / cfg.chirpRateHzPerS;
cfg.sourcePath = configPath;
end
