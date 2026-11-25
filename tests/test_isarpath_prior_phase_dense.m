function checks=test_isarpath_prior_phase_dense(out)
% Independent library, single body-origin point, all 32 pulses of a CPI.
[work,sensor,base]=isar_grid_setup;work.bodyRotationMode="coneSpin";
centre=isar_apply_micro_motion(base,work);centre.scattererCount=1;centre.boundingRadiusM=0;
centre.groupId=1;centre.groupKeys="centre";centre.groupNames="centre";
centre.payload.positionTemplate_B=[0 0 0];centre.payload.alpha0=1;
centre.centerStateFcn=@(value,~,t) struct('position_I',value.payload.motion.positionFcn(t));
cfg=struct('fs',1e9,'fs_scope',200e6,'Tr',100e-6,'fL',11.8e9, ...
    'fH',17.8e9,'kr',6e13,'propagationSpeed',3e8);p=isarpath_settings(cfg);checks=table();
for epoch=[0 4.9]
    prior=struct('epochSec',epoch,'rangeM',10050+3400*epoch+50*epoch^2, ...
        'velocityMps',3400+100*epoch,'accelerationMps2',100,'beatOffsetHz',20e6);
    prior.loStartSec=2*prior.rangeM/3e8-20e6/cfg.kr;t=epoch+(0:31)'/1000;
    cfg.T=ceil((prior.loStartSec+cfg.Tr+16e-6)*cfg.fs/32000)*32000/cfg.fs;
    fast=(0:round(cfg.T*cfg.fs)-1)/cfg.fs;
    request=struct('pulseTimeSec',t,'fastTimeSec',fast,'f0Hz',cfg.fL,'f1Hz',cfg.fH, ...
        'pulseWidthSec',cfg.Tr,'propagationSpeed',cfg.propagationSpeed,'pathLossMode','none');
    opts=struct('Mode','fastTimeExact','AutoFastTimeGate',true,'ReturnGroups',false, ...
        'ReturnPerTarget',false,'FastTimeChunkSize',32768);
    echo=radar.echo_generate_target_lfm(centre,sensor,request,opts);
    br=struct('fastTimeSec',fast,'referenceStartTimeSec',prior.loStartSec,'f0Hz',cfg.fL, ...
        'f1Hz',cfg.fH,'pulseWidthSec',cfg.Tr,'outputSampleRateHz',cfg.fs_scope);
    bo=struct('Mode','ANALYTIC_IF','ReturnPerTarget',false);
    beat=radar.lfm_local_reference_beat({echo.total},br,bo);clear echo
    segment=isarpath_extract_segment(beat.total,cfg,prior,t,p);clear beat
    H=isarpath_range_profiles(segment,cfg,prior,t,p);clear segment
    for alignment=["INTEGER","FRACTIONAL"]
        p.rangeAlignmentMode=alignment;[A,s]=isarpath_align_profiles(H,p);
        [~,gate]=max(sum(abs(A).^2,1));
        for mode=["MODEL_CV","MODEL_CA"]
            p.phaseCompensationMode=mode;ref=isarpath_prior_phase(cfg,prior,t,p,s,mode);
            C=isarpath_phase_image(A,p,ref);z=C(:,gate);phase=unwrap(angle(z));x=t-mean(t);
            residual=phase-polyval(polyfit(x,phase,1),x);rms=sqrt(mean(residual.^2));
            if mode=="MODEL_CA",assert(rms<.02);else,assert(rms>2&&rms<3);end
            checks=[checks;table(epoch,alignment,mode,rms,'VariableNames', ...
                {'epochSec','alignment','phaseMode','nonlinearPhaseRmsRad'})]; %#ok<AGROW>
        end
    end
end
disp(checks);fprintf('Dense native centre control: CA <0.02 rad; CV retains unknown acceleration: PASS\n');
if nargin>0,writetable(checks,fullfile(out,'dense_centre_phase_control.csv'));end
end
