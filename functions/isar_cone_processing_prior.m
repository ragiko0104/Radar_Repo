function prior = isar_cone_processing_prior(rec, k, aHat, beatOffsetHz, cfg)
%ISAR_CONE_PROCESSING_PRIOR Frame prior built only from a measurement record.
%   PRIOR = ISAR_CONE_PROCESSING_PRIOR(REC,K,AHAT,BEATOFFSET,C) builds the
%   processing prior for frame K. Its inputs are the measurement record, the
%   frame index, the fitted acceleration and the recorded local-oscillator
%   offset - there is deliberately NO scene-truth argument, so the processing
%   side cannot silently fall back to R or V from the simulation.
%
%   Window placement uses rec.rangeM(k); the MODEL_CA phase model uses
%   rec.rangeM(k), rec.velocityMps(k) and aHat. The local oscillator uses
%   beatOffsetHz for BOTH the beat offset and loStartSec, so a target's
%   offset can never be mixed with another target's.
%
%   PRIOR fields: epochSec, rangeM, velocityMps, accelerationMps2,
%   beatOffsetHz, loStartSec, sources (labels describing each source).

required = {'t', 'rangeM', 'velocityMps', 'mode'};
if ~isstruct(rec) || ~all(isfield(rec, required)) ...
        || ~ismember(string(rec.mode), ["measured", "exact"]) ...
        || numel(rec.rangeM) ~= numel(rec.t)
    error('ISAR:Chain:BadRecord', ...
        'rec must come from isar_cone_measurement_record (fields t, rangeM, velocityMps, mode).');
end
if k < 1 || k > numel(rec.t)
    error('ISAR:Chain:FrameIndex', 'Frame %d outside the record (%d frames).', ...
        k, numel(rec.t));
end
if ~isscalar(aHat) || ~isfinite(aHat)
    error('ISAR:Chain:Acceleration', 'aHat must be a finite scalar.');
end
if ~isscalar(beatOffsetHz) || ~isfinite(beatOffsetHz)
    error('ISAR:Chain:BeatOffset', 'beatOffsetHz must be a finite scalar.');
end
rangeAtFrame = rec.rangeM(k);
prior = struct();
prior.epochSec = rec.t(k);
prior.rangeM = rangeAtFrame;
prior.velocityMps = rec.velocityMps(k);
prior.accelerationMps2 = aHat;
prior.beatOffsetHz = beatOffsetHz;
prior.loStartSec = 2 * rangeAtFrame / cfg.propagationSpeed ...
    - beatOffsetHz / cfg.kr;
prior.sources = struct( ...
    'range', sprintf('record(%s)', rec.mode), ...
    'velocity', sprintf('record(%s)', rec.mode), ...
    'acceleration', 'range-fit', ...
    'beatOffset', 'golden-config');
end
