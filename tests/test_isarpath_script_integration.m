function test_isarpath_script_integration
% Run the actual parent script, not a duplicated reconstruction of it.
[~,radarState,targets]=isar_grid_setup; %#ok<ASGLU>
cfg_ISARpath=struct('fs',1e9,'fs_scope',200e6,'T',1e-3,'Tr',100e-6, ...
    'Ta',1/4000,'Na',4,'fL',11.8e9,'fH',17.8e9,'kr',6e13,'fc',14.8e9, ...
    'lambda',3e8/14.8e9,'propagationSpeed',3e8,'startTimeSec',0,'echoMode',"fastTimeExact", ...
    'bodyRotationMode',"coneSpin",'spinFrequencyHz',3,'precessionFrequencyHz',1, ...
    'precessionHalfAngleDeg',3,'priorAccelerationMps2',100);
R_isar=10050;v_isar=3400;ISAR_PLOT_ENABLED=false;ISAR_QUALITY_ENABLED=false; %#ok<NASGU>
cfg_ISARpath.postprocess=struct('phaseCompensationMode',"MODEL_CA",'rangeAlignmentMode',"FRACTIONAL");
isarScript=which('multiuav_isarpath.m');
assert(~isempty(isarScript),'RV_ISAR_Path/code/multiuav_isarpath.m not on path.');
run(isarScript);
p=isarpath_settings(cfg_ISARpath);
prior=struct('epochSec',0,'rangeM',10050,'velocityMps',3400,'accelerationMps2',100, ...
    'beatOffsetHz',20e6,'loStartSec',2*10050/3e8-20e6/6e13);
[H,refs]=isarpath_synthesize_profiles(targetsISAR,radarStateISAR,cfg_ISARpath,prior,pulseTimeSecISAR,p);
error=norm(H(:,:,1)-isarPost.hrrfModelCrop{1},'fro')/norm(H(:,:,1),'fro');assert(error<1e-7);
[A,s]=isarpath_align_profiles(H(:,:,1),p);
ref=isarpath_prior_phase(cfg_ISARpath,prior,pulseTimeSecISAR,p,s,"MODEL_CA");
[C,I]=isarpath_phase_image(A,p,ref);
assert(norm(C-isarPost.phaseCompensatedHrrf{1},'fro')/norm(C,'fro')<1e-7);
assert(norm(I-isarPost.image{1},'fro')/norm(I,'fro')<1e-7);
fprintf('Actual multiuav_isarpath.m vs streamed full-chain experiment: PASS (%.3g)\n',error);
end
