function [rangeM, rateMps] = isar_tip_motion_kinematics(t, fp, alphaDeg, ampDeg, fw, psiP, psiW, cfg)
%ISAR_TIP_MOTION_KINEMATICS Analytic cone-tip range and range-rate.
%   [R,DOT] = ISAR_TIP_MOTION_KINEMATICS(T,FP,ALPHA,A,FSW,PSIP,PSIW,CFG)
%   evaluates the tip offset from the cone body origin for the library
%   convention of isar_apply_micro_motion:
%
%       C(t) = C0 * Rz(beta(t)) * Ry(2*pi*fp*t) * Rz(alpha) * Ry(2*pi*fs*t)
%       beta(t) = -A*cos(2*pi*fsw*t + psiW),  C0 = Rz(-135 deg), LOS = +X
%
%   with the tip on the body -y axis at distance cfg.tipLeverM from the
%   origin. Because Rz(alpha) is applied to a body point on the y axis and
%   Ry(2*pi*fs*t) then leaves it invariant, the spin rate fs cannot appear in
%   the tip trajectory: the tip is exactly on the spin axis.
%
%   The LOS projection reduces to
%
%       R(t) = L/sqrt(2) * [ sin(a)*cos(thetaP)*(sin(beta)-cos(beta))
%                            - cos(a)*(sin(beta)+cos(beta)) ]
%
%   where a is the precession half angle, thetaP = 2*pi*fp*t + psiP.
%   R is the projection onto the line of sight relative to the cone centre,
%   so its mean is unobservable and is absorbed by the constant term of the
%   measurement model.
%
%   This closed form matches the stored library truth (tip_geometry_tracks)
%   to 3e-4 cm and 3e-3 Hz; tests/test_isar_tip_motion_mode_identification
%   guards the convention.

if nargin < 8 || isempty(cfg)
    cfg = isar_tip_motion_cfg();
end
t = double(t);
a = deg2rad(alphaDeg);
A = deg2rad(ampDeg);
thetaP = 2 * pi * fp * t + psiP;
w = 2 * pi * fw * t + psiW;
beta = -A * cos(w);
betaDot = A * 2 * pi * fw * sin(w);
sb = sin(beta);
cb = cos(beta);
st = sin(thetaP);
ct = cos(thetaP);
sa = sin(a);
ca = cos(a);
k = cfg.tipLeverM / sqrt(2);
% The tip lever is negative along +y in the body frame; the sign is folded
% into the projection below, which is the form validated against the library.
rangeM = k * (sa .* ct .* (sb - cb) - ca * (sb + cb));
if nargout > 1
    rateMps = k * (sa .* (-2 * pi * fp * st .* (sb - cb) + ct .* (cb + sb) .* betaDot) ...
        - ca * (cb - sb) .* betaDot);
end
end
