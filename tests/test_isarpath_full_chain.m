function test_isarpath_full_chain
% Independent frozen legacy expressions check the refactored default chain.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'functions')); ensure_md_zt_paths();
cfg=struct('fs',1e9,'fs_scope',200e6,'T',1e-3,'Tr',100e-6, ...
    'fL',11.8e9,'fH',17.8e9,'kr',6e13,'propagationSpeed',3e8);
p=isarpath_settings(cfg);p.qualityEnabled=true;
prior=struct('epochSec',.7,'rangeM',12454.5,'velocityMps',3470, ...
    'accelerationMps2',100,'loStartSec',82e-6,'beatOffsetHz',20e6);
t=.7+(0:7)'/4000;rng(41);beat=randn(8,40000)+1j*randn(8,40000);
[segment,starts]=isarpath_extract_segment(beat,cfg,prior,t,p);
expected=beat(:,round((prior.loStartSec-128/200e6)*200e6)+(1:20000));
assert(isequal(segment,expected)&&all(starts==16273));
[H,before,centres]=isarpath_range_profiles(segment,cfg,prior,t,p);
wave=radar.doppler_dispersion_cmp(expected,200e6,cfg.fL,cfg.kr,3470*ones(1,8),'freq');
F=fft(wave,32768,2)/20000;B=fft(expected,32768,2)/20000;legacy=H;
for j=1:8
    centre=round((20e6+2*cfg.kr*3470*(t(j)-t(1))/3e8)/(200e6/32768))+111;
    assert(centre==centres(j));legacy(j,:)=signal.crop_fft_bins(F(j,:),centre,1024);
    assert(isequal(before(j,:),signal.crop_fft_bins(B(j,:),centre,1024)));
end
assert(isequal(H,legacy));
[A,s]=isarpath_align_profiles(H,p);old=H;
for j=2:8
    cc=real(ifft(fft([abs(H(j,:)) zeros(1,1024)]).*conj(fft([abs(H(1,:)) zeros(1,1024)]))));
    candidates=-16:16;[~,k]=max(cc(mod(candidates,2048)+1));shift=candidates(k);
    assert(s(j)==shift);old(j,:)=H(j,mod((0:1023)+shift,1024)+1);
end
assert(isequal(A,old));
[C,I]=isarpath_phase_image(A,p);direct=radar.phase_cmp(old,.12);
assert(norm(C-direct,'fro')/norm(direct,'fro')<1e-14);
assert(norm(I-radar.ISAR_imaging(direct,8,'NO'),'fro')<1e-12);
% The fractional operator must agree with integer shifts and preserve energy.
p.rangeAlignmentMode="FRACTIONAL";
fraction=isarpath_align_profiles(H,p,s);assert(norm(fraction-A,'fro')/norm(A,'fro')<2e-14);
fraction=isarpath_align_profiles(H,p,s+.37);
assert(abs(norm(fraction,'fro')/norm(H,'fro')-1)<1e-14);
% A synthetic known reference verifies the new post-alignment phase branch.
p.phaseCompensationMode="MODEL_CA";reference=ones(8,1024).*exp(1j*(0:7)'.^2);
known=struct('phaseCorrection',exp(-1j*(0:7)'.^2),'method',"test known phase",'accelerationUsedMps2',100);
[C,~,q]=isarpath_phase_image(reference,p,known);
assert(max(abs(C(:)-1))<1e-14&&max(abs(abs(q)-1))<1e-14);
fprintf('ISARpath legacy parity, fractional integer limit/energy, prior-phase position: PASS\n');
end
