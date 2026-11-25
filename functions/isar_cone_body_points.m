function P = isar_cone_body_points()
%ISAR_CONE_BODY_POINTS Body-frame coordinates of the four cone scatterers.
%   P = ISAR_CONE_BODY_POINTS() returns the 4x3 template used by
%   scatterermodel.createLegacyDiskTarget for the unified four-scatterer cone:
%   one tip on the body -y axis and three base points on the circle of radius
%   0.3 m at y = +0.75 m, spaced 120 deg starting at 30 deg in the body x-z
%   plane.
%
%   Row 1 is the tip P1, rows 2:4 are the base points P2..P4. The tip lies on
%   the spin axis, the base points do not, which is why only the base points
%   carry the spin rate.

radiusM = 0.3;
axialM = 1.5;
phi = deg2rad(30) + (0:2)' * (2 * pi / 3);
P = [0, -axialM / 2, 0;
    radiusM * cos(phi), repmat(axialM / 2, 3, 1), radiusM * sin(phi)];
end
