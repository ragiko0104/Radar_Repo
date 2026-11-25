function state=isar_rigid_micro_motion_state(target,radarState,t,indices)
%ISAR_RIGID_MICRO_MOTION_STATE Vectorized analytic spin/coning/Z-swing state.
% Same state contract as the library rigid point-cloud evaluator. No
% stop-and-go approximation or far-field range approximation is used.
t=double(t(:)); nt=numel(t); m=target.payload.motion;
truth=m.conicalMicroMotion; dt=t-truth.referenceTimeSec;
ws=truth.spinRateRadps; wp=truth.coningRateRadps;
psi=ws*dt; phi=wp*dt; a=truth.coningHalfAngleRad;
sy=sin(psi); cy=cos(psi); sp=sin(phi); cp=cos(phi);
sa=sin(a); ca=cos(a); C0=truth.initialBodyToInertial;
% D = Ry(phi)*Rz(a)*Ry(psi).
D=zeros(3,3,nt);
D(1,1,:)=cp.*ca.*cy-sp.*sy; D(1,2,:)=-cp.*sa; D(1,3,:)=cp.*ca.*sy+sp.*cy;
D(2,1,:)=sa.*cy; D(2,2,:)=ca; D(2,3,:)=sa.*sy;
D(3,1,:)=-sp.*ca.*cy-cp.*sy; D(3,2,:)=sp.*sa; D(3,3,:)=-sp.*ca.*sy+cp.*cy;
C=pagemtimes(C0,D);
omegaB=[wp*sa.*cy,repmat(ws+wp*ca,nt,1),wp*sa.*sy];
omegaI=reshape(pagemtimes(C,reshape(omegaB.',3,1,nt)),3,nt).';
% Apply exactly the library's inertial-Z left rotation and add its angular
% velocity. Do not silently drop swing from fastTimeExact rigid evaluation.
if isfield(truth,'enabled')&&truth.enabled.swing
    mid=.5*(truth.swingAngleMinRad+truth.swingAngleMaxRad);
    amplitude=.5*(truth.swingAngleMaxRad-truth.swingAngleMinRad);
    beta=mid-amplitude*cos(truth.swingRateRadps*dt);
    betaDot=amplitude*truth.swingRateRadps*sin(truth.swingRateRadps*dt);
    cb=cos(beta);sb=sin(beta);
    S=zeros(3,3,nt);
    S(1,1,:)=cb;S(1,2,:)=-sb;S(2,1,:)=sb;S(2,2,:)=cb;S(3,3,:)=1;
    C=pagemtimes(S,C);
    omegaI=[cb.*omegaI(:,1)-sb.*omegaI(:,2), ...
        sb.*omegaI(:,1)+cb.*omegaI(:,2),omegaI(:,3)+betaDot];
    omegaB=reshape(pagemtimes(pagetranspose(C),reshape(omegaI.',3,1,nt)),3,nt).';
end
centre=m.positionFcn(t); centreV=m.velocityFcn(t); centreA=m.accelerationFcn(t);
kin=struct('position_I',centre,'velocity_I',centreV,'acceleration_I',centreA, ...
    'C_IB',C,'angularVelocity_I',omegaI,'angularVelocity_B',omegaB);
if isempty(indices),state=kin;return;end
if isa(radarState,'function_handle'),radarState=radarState(t);end
radarP=radar_rows(radarState,{'positionFcn','position_I','positionM'},t);
radarV=radar_rows(radarState,{'velocityFcn','velocity_I','velocityMps'},t);
np=numel(indices); points=target.payload.positionTemplate_B(indices,:);
relative=pagemtimes(points,pagetranspose(C));
position=relative+reshape(centre.',1,3,nt);
velocity=cross(repmat(reshape(omegaI.',1,3,nt),np,1,1),relative,2)+ ...
    reshape(centreV.',1,3,nt);
displacement=position-reshape(radarP.',1,3,nt);
distance=reshape(sqrt(sum(displacement.^2,2)),np,nt);
relativeV=velocity-reshape(radarV.',1,3,nt);
rate=reshape(sum(relativeV.*displacement,2),np,nt)./max(distance,eps);
if isfield(radarState,'wavelengthM'),lambda=radarState.wavelengthM;
elseif isfield(radarState,'carrierFrequencyHz'),lambda=3e8/radarState.carrierFrequencyHz;
else,lambda=Inf;end
signD=-1;if isfield(radarState,'dopplerSign'),signD=radarState.dopplerSign;end
state=struct('rangeM',distance,'alpha',repmat(target.payload.alpha0(indices),1,nt), ...
    'position_I',position,'velocity_I',velocity,'rangeRateMps',rate, ...
    'dopplerHz',signD*2*rate/lambda,'centerPosition_I',centre, ...
    'centerVelocity_I',centreV,'kinematics',kin);
end

function rows=radar_rows(radar,names,t)
raw=[0 0 0];
for k=1:numel(names)
    if isfield(radar,names{k})&&~isempty(radar.(names{k})),raw=radar.(names{k});break;end
end
if isa(raw,'function_handle'),raw=raw(t);end
if isequal(size(raw),[1 3]),rows=repmat(raw,numel(t),1);
elseif isequal(size(raw),[numel(t) 3]),rows=raw;
elseif numel(t)==1&&numel(raw)==3,rows=reshape(raw,1,3);
else,error('RadarProject:ISARMotion:RadarShape','Radar vectors must be Nt-by-3.');end
end
