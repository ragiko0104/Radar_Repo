function response=replay_isar_point_response(beat,cfg,post,targetIndex,velocityPrior)
%REPLAY_ISAR_POINT_RESPONSE Apply the measured scene operators to a point.
% No point-specific alignment, autofocus, or motion-truth correction is fit.
% This linear replay permits an honest per-scatterer PSF/linearity audit.
p=post.parameters;
if p.rangeMigrationMode~="FFT2048_SINGLE" && p.rangeMigrationMode~="NONE" && ...
        p.rangeMigrationMode~="MODEL"
    error('RadarProject:ISARReplay:UnsupportedAlignment', ...
        'Replay supports the default integer-shift FFT2048_SINGLE chain.');
end
nPulse=size(beat,1); nRange=p.hrrfBins; segment=zeros(nPulse,p.timeWindowSamples,'like',beat);
for j=1:nPulse
    first=post.segmentStartSamples(targetIndex,j);
    source=max(1,first):min(size(beat,2),first+p.timeWindowSamples-1);
    segment(j,source-first+1)=beat(j,source);
end
if p.dispersionCompensationEnabled
    segment=radar.doppler_dispersion_cmp(segment,p.inputSampleRateHz, ...
        cfg.fL,cfg.kr,velocityPrior*ones(1,nPulse),'freq');
end
spectra=fft(segment,p.rangeFftPoints,2)/p.timeWindowSamples;
H=complex(zeros(nPulse,nRange,'like',beat));
for j=1:nPulse
    cropped=signal.crop_fft_bins(spectra(j,:),post.cropCenterBins(targetIndex,j),nRange);
    H(j,:)=cropped;
end
alignment="INTEGER";
if isfield(p,'rangeAlignmentMode'),alignment=p.rangeAlignmentMode;end
if alignment=="FRACTIONAL"
    H=isar_physical_fractional_shift(H,post.rangeMigrationBins(targetIndex,:));
else
    for j=1:nPulse
        shift=post.rangeMigrationBins(targetIndex,j);
        H(j,:)=H(j,mod((0:nRange-1)+shift,nRange)+1);
    end
end
H=H.*post.phaseCorrection{targetIndex};
response=struct('hrrf',H,'image',radar.ISAR_imaging(H,p.slowTimeZeroPaddingFactor,char(p.windowType)), ...
    'rangeAxisM',post.rangeAxisM{targetIndex},'crossAxisM',post.crossAxisM{targetIndex});
end
