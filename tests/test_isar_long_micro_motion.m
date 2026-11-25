function test_isar_long_micro_motion()
% Independent RF-product check, intrapulse error, static tip and RCMC replay.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'functions')); ensure_md_zt_paths();
spec=struct('Ntip',1,'radiusM',.3,'phase0Rad',pi/6, ...
    'axialSeparationM',1.5,'alpha0',[1;.3;.3;.3]);
motion=struct('type','hover','p0_I',[10050 0 0],'attitude0YprRad',[-3*pi/4 0 0]);
base=scatterermodel.createLegacyDiskTarget(spec,motion);
radar=struct('position_I',[0 0 0],'velocity_I',[0 0 0], ...
    'carrierFrequencyHz',14.8e9,'wavelengthM',3e8/14.8e9,'dopplerSign',-1);
cfg=struct('sampleRateHz',200e6,'pulseWidthSec',100e-6,'propagationSpeed',3e8, ...
    'fL',11.8e9,'fH',17.8e9,'beatOffsetHz',20e6,'beatGuardHz',90e6, ...
    'rangeFftPoints',32768,'rangeBins',1024,'centerRangeM',10050, ...
    'slowPadding',8,'rcmcSearchBins',16);
motionCfg=struct('spinFrequencyHz',1,'precessionFrequencyHz',2, ...
    'precessionHalfAngleDeg',10,'swingFrequencyHz',1,'swingAngleMinDeg',-5,'swingAngleMaxDeg',5);
for mode=["spin","coneSpin","swingSpin","coneSwingSpin"]
    motionCfg.bodyRotationMode=mode;target=isar_apply_micro_motion(base,motionCfg);
    times=[0;.123;.777;4.999];
    [H,info]=isar_streamed_stretch_profiles(target,radar,times,cfg);
    assert(info.maxRangeInterpolationErrorM<1e-8&&info.maxAbsBeatHz<21e6);
    assert(max(info.loDelaySec)-min(info.loDelaySec)==0);
    % Direct full-fast-time target evaluation and separate analytic RF/LO
    % exponentials provide an independent reference for the stable IF phase.
    u=(0:round(cfg.sampleRateHz*cfg.pulseWidthSec)-1)/cfg.sampleRateHz;
    k=(cfg.fH-cfg.fL)/cfg.pulseWidthSec;
    reference=exp(1j*2*pi*(cfg.fL*u+.5*k*u.^2));
    for n=1:numel(times)
        tau=info.loDelaySec(n);
        direct=scatterermodel.evaluateTarget(target,radar,times(n)+tau+u);
        for j=1:4
            delayed=tau+u-2*direct.rangeM(j,:)/cfg.propagationSpeed;
            echo=target.payload.alpha0(j)*double(delayed>=0&delayed<=cfg.pulseWidthSec).* ...
                exp(1j*2*pi*(cfg.fL*delayed+.5*k*delayed.^2));
            beat=echo.*conj(reference).*exp(1j*2*pi*cfg.beatOffsetHz*u);
            expected=fftshift(fft(beat,cfg.rangeFftPoints))/numel(u);
            crop=cfg.rangeFftPoints/2+1+(-cfg.rangeBins/2:cfg.rangeBins/2-1);
            relative=norm(double(H(n,:,j))-expected(crop))/norm(expected(crop));
            assert(relative<2e-6,'Direct analytic RF mismatch: %s, %.3g',mode,relative);
        end
    end
    [block,~]=isar_streamed_stretch_profiles(target,radar,(0:31)'*.001,cfg);
    frame=isar_long_rd_frame(block,cfg,1000,target,radar,(0:31)'*.001);
    assert(frame.replayError<1e-12);
    if mode=="spin"
        tip=double(block(:,:,1));assert(norm(tip-tip(1,:),'fro')<1e-12);
        p=frame.metrics(frame.metrics.method=="noRCMC"&frame.metrics.pointIndex==1,:);
        assert(p.peakDopplerHz==0&&abs(p.doppler3dBWidthHz-27.6842)<.03);
    end
    % Default frequencies have an exact 1 s period; cache reuse must agree
    % with actual samples far into the requested 5 s interval.
    [repeat,~]=isar_streamed_stretch_profiles(target,radar,times+1,cfg);
    assert(norm(double(H(:))-double(repeat(:)))/norm(double(H(:)))<2e-6);
end
fprintf('Long micro-motion RF parity, intrapulse geometry and coherent RCMC tests passed.\n');
end
