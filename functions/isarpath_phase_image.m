function [corrected,image,correction,info]=isarpath_phase_image(H,p,reference)
%ISARPATH_PHASE_IMAGE Common original ISARpath phase-compensation position.
% REFERENCE is the explicit kinematic phase from isarpath_prior_phase.
% It contains only supplied priors and known processing coordinates.
scene=H(:,:,1);info=struct('mode',p.phaseCompensationMode,'referenceBin',NaN,'minimumReferenceAmplitudeRatio',NaN);
correction=ones(size(scene,1),1);
if p.phaseCompensationMode=="VARIANCE_UNIT"
    changed=radar.phase_cmp(scene,p.phaseThreshold);
    if all(isfinite(changed(:)))
        power=sum(abs(scene).^2,2);good=power>0;
        correction(good)=sum(changed(good,:).*conj(scene(good,:)),2)./power(good);
    else
        info.mode="NONE_INVALID_REFERENCE";
    end
elseif ismember(p.phaseCompensationMode,["MODEL_CV","MODEL_CA"])
    assert(isstruct(reference),'Use isarpath_prior_phase; do not take phase from a model HRRP null.');
    correction=reference.phaseCorrection;
    assert(isequal(size(correction),[size(H,1),1]));
    assert(max(abs(abs(correction)-1))<1e-12);
    info.referenceMethod=reference.method;
    info.accelerationUsedMps2=reference.accelerationUsedMps2;
end
corrected=H.*correction;
if any(abs(corrected(:,:,1))>0,'all')
    image=radar.ISAR_imaging(corrected(:,:,1),p.slowTimeZeroPaddingFactor,char(p.windowType));
else
    image=zeros(size(H,1)*p.slowTimeZeroPaddingFactor,size(H,2));
end
end
