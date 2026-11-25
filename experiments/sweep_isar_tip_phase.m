function results = sweep_isar_tip_phase(outputDirectory)
% Blind 0.5 Hz tip fits across angle, initial phase and independent noise.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'functions'));
if nargin < 1
    outputDirectory = fullfile(root, 'results', 'isar_tip_phase');
end
if ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end
cfg = isar_tip_motion_cfg();
points = isar_cone_body_points();
t = 0.0156 + (0:49)' * cfg.sampleStepSec;
[angle, phase, trial] = ndgrid([3 5 6 8 20], [-pi+0.02 -pi/2 0 pi/2 pi-0.02], 1:2);
results = table(angle(:), phase(:), trial(:), ...
    'VariableNames', {'alphaDeg', 'phaseRad', 'trial'});
n = height(results);
results.noiseSeed = (20261122 + (0:n-1))';
results.model = strings(n, 1);
results.errFpHz = zeros(n, 1);
results.errAlphaDeg = zeros(n, 1);
previousRng = rng;
cleanup = onCleanup(@() rng(previousRng)); %#ok<NASGU>
for i = 1:n
    rng(results.noiseSeed(i), 'twister');
    p = [1 0.5 results.alphaDeg(i) 0 0 results.phaseRad(i) 0 0];
    y = isar_simulate_cone_peaks(t, points(1,:), p, cfg, 0.2, [5.15 0.14]);
    id = isar_identify_tip_mode(t, y, cfg);
    results.model(i) = id.model;
    results.errFpHz(i) = id.selected.theta(1) - 0.5;
    results.errAlphaDeg(i) = id.selected.theta(2) - results.alphaDeg(i);
end
results.ok = results.model == "H1" & abs(results.errFpHz) < 0.01 & abs(results.errAlphaDeg) < 0.1;
writetable(results, fullfile(outputDirectory, 'tip_phase.csv'));
fprintf('Tip phase sweep: %d/%d pass; max errors %.6f Hz, %.6f deg.\n', ...
    nnz(results.ok), n, max(abs(results.errFpHz)), max(abs(results.errAlphaDeg)));
disp(results(~results.ok,:));
end
