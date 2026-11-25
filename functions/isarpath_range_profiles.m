function [H,before,centres]=isarpath_range_profiles(segment,cfg,prior,pulseTime,p)
%ISARPATH_RANGE_PROFILES Original dispersion compensation, FFT, moving crop.
% The CV envelope model is deliberately identical for EVERY phase branch.
% CA is introduced only at the post-alignment phase-compensation position.
before=[];
if p.qualityEnabled,beforeFft=fft(segment,p.rangeFftPoints,2)/p.timeWindowSamples;end
if p.dispersionCompensationEnabled
    segment=radar.doppler_dispersion_cmp(segment,cfg.fs_scope,cfg.fL,cfg.kr, ...
        prior.velocityMps*ones(1,size(segment,1)),'freq');
end
F=fft(segment,p.rangeFftPoints,2)/p.timeWindowSamples;
f=prior.beatOffsetHz+2*cfg.kr*prior.velocityMps*(pulseTime(:)-prior.epochSec)/cfg.propagationSpeed;
centres=round(f/p.frequencyBinSpacingHz)+1+p.cropCenterOffsetBinsEffective;
H=complex(zeros(size(segment,1),p.hrrfBins,'like',segment));
if p.qualityEnabled,before=H;end
for j=1:size(H,1)
    H(j,:)=signal.crop_fft_bins(F(j,:),centres(j),p.hrrfBins);
    if p.qualityEnabled,before(j,:)=signal.crop_fft_bins(beforeFft(j,:),centres(j),p.hrrfBins);end
end
end
