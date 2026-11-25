function [aligned,s,scores]=isarpath_align_profiles(H,p,forcedShift)
%ISARPATH_ALIGN_PROFILES Original single-pass first-pulse envelope reference.
% FRACTIONAL adds a 3-point parabolic peak fit; it uses the original positive
% fast-time origin. No scalar 'shift phase recovery' factor is introduced.
np=size(H,1);nr=size(H,2);scores=[];
if nargin>=3
    s=forcedShift(:);
elseif p.rangeMigrationMode=="FFT2048_SINGLE"
    assert(p.rcmcCorrelationPoints==2*nr);
    candidates=-p.rcmcSearchRadiusBins:p.rcmcSearchRadiusBins;
    indices=mod(candidates,p.rcmcCorrelationPoints)+1;
    spectrum=fft([abs(double(H(:,:,1))) zeros(np,nr)],p.rcmcCorrelationPoints,2);
    corr=real(ifft(spectrum.*conj(spectrum(1,:)),p.rcmcCorrelationPoints,2));
    scores=corr(:,indices);s=zeros(np,1);
    for j=2:np
        [~,pk]=max(scores(j,:));s(j)=candidates(pk);
        if p.rangeAlignmentMode=="FRACTIONAL" && pk>1 && pk<numel(candidates)
            yy=scores(j,pk-1:pk+1);den=yy(1)-2*yy(2)+yy(3);
            if den<0,s(j)=s(j)+max(-.5,min(.5,.5*(yy(1)-yy(3))/den));end
        end
    end
elseif ismember(p.rangeMigrationMode,["XCORR","MODEL_XCORR"])
    assert(ismatrix(H));
    wave=nr*ifft(ifftshift(H,2),nr,2);
    [~,aligned,s]=radar.range_migration_alignment(wave,H,'freq',1:nr,p.rangeAlignmentMode=="FRACTIONAL");
    return
else
    s=zeros(np,1);
end
if p.rangeAlignmentMode=="FRACTIONAL"
    aligned=isar_physical_fractional_shift(H,s);
else
    aligned=H;
    for channel=1:size(H,3)
        for j=1:np,aligned(j,:,channel)=H(j,mod((0:nr-1)+round(s(j)),nr)+1,channel);end
    end
end
end
