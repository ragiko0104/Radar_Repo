function output=isar_physical_fractional_shift(input,shiftBins)
%ISAR_PHYSICAL_FRACTIONAL_SHIFT H(k+s), preserving the existing u>=0 IF origin.
% Inverse range FFT recovers positive physical fast time, not signed wrap.
% No midpoint subtraction and no additional per-pulse scalar phase factor.
[np,nr,nc]=size(input);s=double(shiftBins(:));assert(numel(s)==np);
multiplier=exp(-1j*2*pi*s.*(0:nr-1)/nr);
output=complex(zeros(size(input),'like',input));
for channel=1:nc,output(:,:,channel)=fft(ifft(input(:,:,channel),[],2).*multiplier,[],2);end
end
