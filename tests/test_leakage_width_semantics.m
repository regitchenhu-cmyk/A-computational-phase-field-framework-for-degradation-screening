function test_leakage_width_semantics()
%TEST_LEAKAGE_WIDTH_SEMANTICS Darcy width scales Q_D, not the y-integrated Q_P.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
addpath(genpath(fullfile(root, 'Models')));
addpath(genpath(fullfile(root, 'Shapes')));
addpath(genpath(fullfile(root, 'PostProcessing')));

[qD1, qP1] = run_width_case(root, 1.0, 'Width1');
[qD2, qP2] = run_width_case(root, 2.0, 'Width2');

assert(qD1 > 0 && qP1 > 0, 'Leakage regression case returned a zero rate.');
assert(abs(qD2 / qD1 - 2.0) < 1e-10, ...
    'Darcy leakage did not scale linearly with extrusion width.');
assert(abs(qP2 / qP1 - 1.0) < 1e-10, ...
    'Aperture leakage incorrectly retained a second width factor.');

fprintf('PASS test_leakage_width_semantics: QD ratio=%.12g, QP ratio=%.12g.\n', ...
    qD2 / qD1, qP2 / qP1);
end


function [qD, qP] = run_width_case(root, width, name)
resultRoot = fullfile(root, 'Results_Regression_LeakageWidth_20260824');
opts = calibrated_degradation_options(resultRoot);
opts.sname = name;
opts.mesh = struct('Nx', 4, 'Ny', 6);
opts.nMax = 1;
opts.tmaxDays = 1;
opts.initialDt = 30;
opts.maxDt = 30;
opts.plotEvery = inf;
opts.saveEvery = inf;
opts.verboseEvery = inf;
opts.useParpool = false;
opts.maxThreads = 1;
opts.cleanResults = true;
opts.forceRestart = false;
opts.outOfPlaneWidth = width;

[~, ~, results] = main_seal(1, opts);
qD = results.leakage_Q(end);
qP = results.leakage_Q_poiseuille(end);
end
