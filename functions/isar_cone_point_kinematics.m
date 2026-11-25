function [rangeM, rateMps] = isar_cone_point_kinematics(t, bodyPoint, p, cfg)
%ISAR_CONE_POINT_KINEMATICS Range and range-rate of an arbitrary cone point.
%   [R,DOT] = ISAR_CONE_POINT_KINEMATICS(T,BODYPOINT,P,CFG) evaluates the
%   line-of-sight projection of one scatterer under the library attitude
%
%       C(t) = C0 * Rz(beta) * Ry(thp) * Rz(alpha) * Ry(ths)
%       beta = -A*cos(2*pi*fw*t + psiW)     thp = 2*pi*fp*t + psiP
%       ths = 2*pi*fs*t + psiS              C0 = Rz(-135 deg), LOS = +X
%
%   P is [fs fp alphaDeg ampDeg fw psiP psiW psiS].
%
%   Unlike isar_tip_motion_kinematics this keeps the spin factor: the tip is
%   on the rotation axis so Ry(ths) fixes it, but a base point at radius
%   0.3 m sweeps a circle and its projection oscillates at fs (plus cross
%   terms with fp). That is the observable used to identify the spin rate.
%
%   Expanding the chain gives the compact form used here,
%
%       u = Rz(alpha) Ry(ths) b,  v = Ry(thp) u,
%       R = (v1*(sb-cb) + v2*(sb+cb))/sqrt(2),
%
%   together with its exact time derivative, so no finite differences are
%   needed and the model stays smooth in the parameters.

if nargin < 4 || isempty(cfg)
    cfg = isar_tip_motion_cfg();
end
% One scatterer only. A 4x3 body-point table silently linear-indexes to
% [P(1,1) P(2,1) P(3,1)], which is a valid-length but meaningless position,
% so reject anything that is not a single 3-vector.
if ~isnumeric(bodyPoint) || ~isreal(bodyPoint) || numel(bodyPoint) ~= 3 ...
        || ~isvector(bodyPoint) || any(~isfinite(bodyPoint(:)))
    error('ISAR:Cone:BodyPointShape', ...
        ['bodyPoint must be one 3-element position vector; received a ' ...
        '%s of size %s. Pass bodyPoints(k,:) for scatterer k.'], ...
        class(bodyPoint), mat2str(size(bodyPoint)));
end
bodyPoint = double(bodyPoint(:)).';
t = double(t(:));
px = bodyPoint(1);
py = bodyPoint(2);
pz = bodyPoint(3);
fs = p(1);
fp = p(2);
a = deg2rad(p(3));
A = deg2rad(p(4));
fw = p(5);
thp = 2 * pi * fp * t + p(6);
ths = 2 * pi * fs * t + p(8);
w = 2 * pi * fw * t + p(7);
beta = -A * cos(w);
betaDot = A * 2 * pi * fw * sin(w);
sb = sin(beta);
cb = cos(beta);
cs = cos(ths);
ss = sin(ths);
cp = cos(thp);
sp = sin(thp);
ca = cos(a);
sa = sin(a);
x1 = px * cs + pz * ss;
u1 = x1 * ca - py * sa;
u2 = x1 * sa + py * ca;
u3 = -px * ss + pz * cs;
v1 = u1 .* cp + u3 .* sp;
v2 = u2;
rangeM = (v1 .* (sb - cb) + v2 .* (sb + cb)) / sqrt(2);
if nargout > 1
    thsDot = 2 * pi * fs;
    thpDot = 2 * pi * fp;
    du1 = u3 * thsDot * ca;
    du2 = u3 * thsDot * sa;
    du3 = -x1 * thsDot;
    dv1 = du1 .* cp - u1 .* sp * thpDot + du3 .* sp + u3 .* cp * thpDot;
    dv2 = du2;
    rateMps = (dv1 .* (sb - cb) + dv2 .* (sb + cb) ...
        + betaDot .* (v1 .* (cb + sb) + v2 .* (cb - sb))) / sqrt(2);
end
end
