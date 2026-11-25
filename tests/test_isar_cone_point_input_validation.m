function test_isar_cone_point_input_validation()
%TEST_ISAR_CONE_POINT_INPUT_VALIDATION Single-point entry must reject tables.
%   The 4x3 body-point table used to be accepted silently: MATLAB linear
%   indexing then read [P(1,1) P(2,1) P(3,1)], a valid-length but meaningless
%   position. Guard against a regression.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'functions'));
cfg = isar_tip_motion_cfg();
P = isar_cone_body_points();
p = [1, 0.7, 10, 0, 1, 0, 0, 0];
t = (0:0.01:0.2)';

% 1. a 4x3 table must be rejected, not silently linear-indexed
threw = false;
try
    isar_cone_point_kinematics(t, P, p, cfg);
catch err
    threw = strcmp(err.identifier, 'ISAR:Cone:BodyPointShape');
end
assert(threw, 'A 4x3 body-point table must raise ISAR:Cone:BodyPointShape.');

% 2. other wrong shapes are rejected too
for bad = {P(1:2, :), [1 2 3 4], 1, 'abc'}
    threw = false;
    try
        isar_cone_point_kinematics(t, bad{1}, p, cfg);
    catch err
        threw = strcmp(err.identifier, 'ISAR:Cone:BodyPointShape');
    end
    assert(threw, 'Malformed bodyPoint must be rejected.');
end

% 3. legal calls still work and agree across orientation
[rRow, vRow] = isar_cone_point_kinematics(t, P(1, :), p, cfg);
[rCol, vCol] = isar_cone_point_kinematics(t, P(1, :)', p, cfg);
assert(isequal(size(rRow), size(t)), 'Range must follow the time vector shape.');
assert(max(abs(rRow - rCol)) == 0, 'Row and column points must agree.');
assert(max(abs(vRow - vCol)) == 0, 'Row and column rates must agree.');

% 4. the rejected call would have produced a different point had it passed
%    (documents what the bug did: first column of three rows)
wrongPoint = [P(1, 1), P(2, 1), P(3, 1)];
rWrong = isar_cone_point_kinematics(t, wrongPoint, p, cfg);
assert(max(abs(rWrong - rRow)) > 1e-6, ...
    'The silently-indexed point differs from the tip, so the guard matters.');

fprintf('test_isar_cone_point_input_validation passed.\n');
end
