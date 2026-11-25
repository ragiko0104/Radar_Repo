function output=isar_symmetric_fractional_shift(input,shiftBins)
%ISAR_SYMMETRIC_FRACTIONAL_SHIFT Periodic H(k+s), using signed dual frequencies.
% No fast-time-midpoint rephasing or additional scalar phase factor.
% Integer s is precisely the same operator as circshift(H,-s,2).
[np,nr,nc]=size(input);s=double(shiftBins(:));assert(numel(s)==np&&all(isfinite(s)));
dual=(-floor(nr/2):ceil(nr/2)-1)/nr;
multiplier=exp(1j*2*pi*s.*dual);
output=complex(zeros(size(input),'like',input));
for channel=1:nc
    spectrum=fftshift(fft(input(:,:,channel),[],2),2);
    output(:,:,channel)=ifft(ifftshift(spectrum.*multiplier,2),[],2);
end
end
