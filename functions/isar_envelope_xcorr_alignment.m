function [aligned,shiftBins,info]=isar_envelope_xcorr_alignment(profiles,searchBins,mode)
%ISAR_ENVELOPE_XCORR_ALIGNMENT Shared three-pass envelope alignment only.
% Channel 1 is the actual received scene; replay ONE estimate to every channel.
assert(ismember(string(mode),["integer","fractional","physical"]));
np=size(profiles,1);nr=size(profiles,2);work=profiles(:,:,1);
shiftBins=zeros(np,1);passShifts=zeros(np,3);boundaryCount=0;
for iteration=1:3
    magnitude=abs(work);
    if iteration==1,reference=magnitude(1,:);else,reference=mean(magnitude,1);end
    spectrum=fft([magnitude zeros(np,nr)],2*nr,2);
    refSpectrum=fft([reference zeros(1,nr)],2*nr);
    correlation=real(ifft(spectrum.*conj(refSpectrum),2*nr,2));
    candidates=-searchBins:searchBins;
    [~,best]=max(correlation(:,mod(candidates,2*nr)+1),[],2);
    step=reshape(candidates(best),[],1);
    if string(mode)~="integer"
        for pulse=1:np
            if abs(step(pulse))==searchBins,boundaryCount=boundaryCount+1;continue;end
            y=correlation(pulse,mod(step(pulse)+(-1:1),2*nr)+1);
            denominator=y(1)-2*y(2)+y(3);
            if denominator<0
                offset=.5*(y(1)-y(3))/denominator;
                assert(abs(offset)<=.5+1e-10);step(pulse)=step(pulse)+offset;
            end
        end
        work=fractional_shift(work,step,mode);
    else
        work=integer_shift(work,step);
    end
    passShifts(:,iteration)=step;shiftBins=shiftBins+step;
end
if string(mode)~="integer",aligned=fractional_shift(profiles,shiftBins,mode);
else,aligned=integer_shift(profiles,shiftBins);end
error=norm(aligned(:,:,1)-work,'fro')/max(norm(work,'fro'),eps);
assert(error<1e-12);
sceneError=norm(sum(aligned(:,:,2:end),3)-aligned(:,:,1),'fro')/max(norm(aligned(:,:,1),'fro'),eps);
assert(sceneError<1e-12);
info=struct('mode',string(mode),'passShifts',passShifts,'boundaryCount',boundaryCount, ...
    'replayRelativeError',error,'sceneSumRelativeError',sceneError,'extraPhaseFactorApplied',false);
end

function output=fractional_shift(input,s,mode)
if string(mode)=="physical",output=isar_physical_fractional_shift(input,s);
else,output=isar_symmetric_fractional_shift(input,s);end
end

function out=integer_shift(H,s)
[np,nr,nc]=size(H);indices=mod((0:nr-1)+s,nr)+1;out=H;
for channel=1:nc
    A=H(:,:,channel);out(:,:,channel)=A(sub2ind([np nr],repmat((1:np)',1,nr),indices));
end
end
