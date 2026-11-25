function [segment,starts]=isarpath_extract_segment(beat,cfg,prior,pulseTime,p)
%ISARPATH_EXTRACT_SEGMENT The original fixed/moving receiver time window.
dt=pulseTime(:)-prior.epochSec;
starts=round((prior.loStartSec-p.timeMarginSamples/cfg.fs_scope)*cfg.fs_scope)+ones(size(dt));
if p.timeWindowMode=="MOVING_PULSE"
    starts=round((prior.loStartSec+2*prior.velocityMps*dt/cfg.propagationSpeed- ...
        p.timeMarginSamples/cfg.fs_scope)*cfg.fs_scope)+1;
end
segment=zeros(numel(dt),p.timeWindowSamples,'like',beat);
for j=1:numel(dt)
    source=max(1,starts(j)):min(size(beat,2),starts(j)+p.timeWindowSamples-1);
    segment(j,source-starts(j)+1)=beat(j,source);
end
end
