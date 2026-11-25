function geometry=isar_cross_range_geometry(targets,radar,cfg,pulseTime)
%ISAR_CROSS_RANGE_GEOMETRY Target-specific local Doppler-to-metre conversion.
% Evaluate at received chirp midpoint, including apparent rotation from a
% changing LOS. Precessing 3-D motion has only a local image-plane scale;
% plane drift and rate variation are retained as explicit diagnostics.
geometry=repmat(struct(),numel(targets),1);
for k=1:numel(targets)
    target=targets(k);
    first=target.centerStateFcn(target,radar,pulseTime(1));
    radarFirst=radar_position(radar,pulseTime(1));
    echoOffset=2*norm(first.position_I(1,:)-radarFirst)/cfg.propagationSpeed+cfg.Tr/2;
    query=pulseTime(:)+echoOffset; mid=mean(query);
    [rate,crossHat,los,omega]=local_geometry(target,radar,mid);
    if ~(isfinite(rate)&&rate>1e-10)
        error('RadarProject:ISARMotion:UnobservableAspect', ...
            'Target %d has no observable rotation perpendicular to LOS.',k);
    end
    rates=zeros(numel(query),1); directions=zeros(numel(query),3);
    for j=1:numel(query)
        [rates(j),directions(j,:)]=local_geometry(target,radar,query(j));
    end
    drift=acosd(max(-1,min(1,directions*crossHat.')));
    geometry(k).referenceTimeSec=mid;
    geometry(k).effectiveRateRadps=rate;
    geometry(k).crossRangeHat_I=crossHat;
    geometry(k).losHat_I=los;
    geometry(k).effectiveAngularVelocity_I=omega;
    geometry(k).pulseEffectiveRateRadps=rates;
    geometry(k).maximumPlaneDriftDeg=max(drift);
    geometry(k).relativeRateVariation=(max(rates)-min(rates))/rate;
    geometry(k).localProjectionOnly=max(drift)>1e-6 || ...
        geometry(k).relativeRateVariation>1e-6;
    geometry(k).nominalResolutionM=cfg.lambda/(2*rate*numel(pulseTime)*cfg.Ta);
end
end

function [rate,hat,los,omega]=local_geometry(target,radar,t)
kin=target.centerStateFcn(target,radar,t);
position=radar_position(radar,t);velocity=radar_velocity(radar,t);
r=kin.position_I(1,:)-position;distance=norm(r);los=r/max(distance,eps);
v=kin.velocity_I(1,:)-velocity;
losRate=(v-dot(v,los)*los)/max(distance,eps);
omega=kin.angularVelocity_I(1,:)-cross(los,losRate);
tangent=cross(los,omega);rate=norm(tangent);hat=tangent/max(rate,eps);
end

function value=radar_position(radar,t)
value=radar_value(radar,t,{'positionFcn','position_I','positionM'});
end
function value=radar_velocity(radar,t)
value=radar_value(radar,t,{'velocityFcn','velocity_I','velocityMps'});
end
function value=radar_value(radar,t,names)
if isa(radar,'function_handle'),radar=radar(t);end
value=[0 0 0];
for k=1:numel(names)
    if isfield(radar,names{k})&&~isempty(radar.(names{k})),value=radar.(names{k});break;end
end
if isa(value,'function_handle'),value=value(t);end
value=reshape(value,1,3);
end
