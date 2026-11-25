function test_isar_fractional_alignment()
% Independent Fourier-series, integer-limit and sub-bin-estimation checks.
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'functions'));
n=256;k=0:n-1;harmonics=[-37 -9 0 23 51];a=[.3+.2j -.1+.4j 1 .2-.3j -.2j];
s=[0;1;-2;.37;-1.29;3.41];H=zeros(numel(s),n);expected=H;
for j=1:numel(a)
    H=H+a(j)*exp(1j*2*pi*harmonics(j)*k/n);
    expected=expected+a(j)*exp(1j*2*pi*harmonics(j)*(k+s)/n);
end
actual=isar_symmetric_fractional_shift(H,s);
assert(max(abs(actual-expected),[],'all')<1e-12,'Signed fractional-shift convention failed.');
back=isar_symmetric_fractional_shift(actual,-s);assert(max(abs(back-H),[],'all')<1e-12);
for j=1:3,assert(max(abs(actual(j,:)-circshift(H(j,:),-s(j),2)))<1e-12);end
template=exp(-((k-110)/4).^2/2).*exp(1j*.11*2*pi*k);
drift=[0;.27;.83;1.41;-.36;-.91];single=isar_symmetric_fractional_shift(repmat(template,6,1),-drift);
profiles=cat(3,single,single*.6,single*.4);
[~,estimate,info]=isar_envelope_xcorr_alignment(profiles,8,"fractional");
% A constant final frame offset is free; compare relative displacements.
error=max(abs((estimate-estimate(1))-(drift-drift(1))));
assert(error<.02&&~info.extraPhaseFactorApplied);
nf=137;n=256;r=31.27;shift=.43;phi=.37;m=0:nf-1;bins=0:n-1;
H=fft(exp(1j*phi).*exp(1j*2*pi*r*m/n),n)/nf;
actual=isar_physical_fractional_shift(H,shift);d=r-shift-bins;
expected=exp(1j*phi).*exp(1j*pi*(nf-1)*d/n).*sin(pi*nf*d/n)./(nf*sin(pi*d/n));
assert(max(abs(actual-expected))<1e-12,'Original-time-origin shift disagrees with exact finite-aperture DFT.');
fprintf('PASS: analytic signed Fourier shift, inverse, integer limit, xcorr relative error %.6f bins.\n',error);
fprintf('PASS: original-time-origin fractional shift agrees with analytic finite-aperture point response.\n');
end
