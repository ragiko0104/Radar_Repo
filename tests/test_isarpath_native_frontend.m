function test_isarpath_native_frontend
[b,sensor,base]=isar_grid_setup;
cfg=struct('fs',1e9,'fs_scope',200e6,'T',1e-3,'Tr',100e-6,'fL',11.8e9, ...
    'fH',17.8e9,'kr',6e13,'propagationSpeed',3e8);
b.bodyRotationMode="coneSpin";b.spinFrequencyHz=3;b.precessionFrequencyHz=1;b.precessionHalfAngleDeg=3;
target=isar_apply_micro_motion(base,b);p=isarpath_settings(cfg);t=[0;.03175];
prior=struct('epochSec',0,'rangeM',10050,'velocityMps',3400,'accelerationMps2',100, ...
    'beatOffsetHz',20e6,'loStartSec',2*10050/3e8-20e6/6e13);
timer=tic;[H,R,a]=isarpath_synthesize_profiles(target,sensor,cfg,prior,t,p);
fprintf('Native 2-pulse cached frontend: %.2f s, sum error %.3g\n',toc(timer),a.coherentSumError);
cfg.cacheExactGeometry=false;
[uncached,refs]=isarpath_synthesize_profiles(target,sensor,cfg,prior,t,p);
assert(norm(H(:)-uncached(:))/norm(H(:))<1e-7);
assert(isequal(R,refs));
% A body-origin target independently checks prior RF synthesis, folding sign,
% LO quantisation, filters, resampling and all subsequent moving crop phases.
centre=target;centre.payload.positionTemplate_B=zeros(4,3);centre.payload.alpha0=[1;0;0;0];
centre.boundingRadiusM=0;cfg.cacheExactGeometry=true;
[known,ref]=isarpath_synthesize_profiles(centre,sensor,cfg,prior,t,p);
actual=known(:,:,1);expected=ref(:,:,2);
error=norm(actual(:)-expected(:))/norm(actual(:));assert(error<1e-7);
p.phaseCompensationMode="MODEL_CA";
fprintf('Cached/un-cached waveform equality and independently generated CA reference: PASS (%.3g)\n',error);
% The production branch computes phase directly, without taking the phase of
% a range gate whose model envelope might cross zero under CV mismatch.
prior.accelerationMps2=100;
for alignment=["INTEGER","FRACTIONAL"]
    p.rangeAlignmentMode=alignment;[A,shift]=isarpath_align_profiles(actual,p);
    direct=isarpath_prior_phase(cfg,prior,t,p,shift,"MODEL_CA");
    [C,~]=isarpath_phase_image(A,p,direct);[~,gate]=max(abs(A(1,:)));
    residual=angle(C(:,gate).*conj(C(1,gate)));
    assert(max(abs(residual))<.02,'Direct kinematic phase does not match native centre echo.');
end
one=isarpath_prior_phase(cfg,prior,t,p,shift,"MODEL_CV");
prior.accelerationMps2=0;two=isarpath_prior_phase(cfg,prior,t,p,shift,"MODEL_CV");
assert(isequal(one.phaseCorrection,two.phaseCorrection));
fprintf('Direct prior-phase gauge, native centre phase residual <0.02 rad, CV has no acceleration leakage: PASS\n');
end
