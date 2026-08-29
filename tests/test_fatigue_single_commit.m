function test_fatigue_single_commit()
%TEST_FATIGUE_SINGLE_COMMIT Regression test for idempotent fatigue assembly.
%
% The fatigue increment must be based on the committed integration-point
% state even when Newton/staggered iterations assemble the same physical
% step repeatedly.  With the energetic term disabled and a unit pressure
% ratio, every integration point should receive exactly baseRate*dN once.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
addpath(genpath(fullfile(root, 'Models')));
addpath(genpath(fullfile(root, 'Shapes')));
addpath(genpath(fullfile(root, 'PostProcessing')));

resultRoot = fullfile(root, 'Results_Regression_FatigueCommit_20260823');
opts = calibrated_degradation_options(resultRoot);
opts.sname = 'SingleStep';
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

opts.P_max = 21e6;
opts.P_min = 0.5e6;
opts.f_cycle = 0.05;
opts.W0 = 1e30;
opts.beta = 2;
opts.fatigueP_ref = opts.P_max - opts.P_min;
opts.fatiguePressureExponent = 0;
opts.fatigueEnergyRate = 1.0;
opts.fatigueBaseRate = 1e-4;
opts.fatigueMaxRate = 1e-3;
opts.contactFatigueCoeff = 0;

[physics, ~, results] = main_seal(1, opts);
fatigue = physics.models{8};
assert(abs(fatigue.energyRateScale - 1.0) < 1e-12, ...
    'Energetic fatigue-rate scale was not propagated to the model.');
expectedCycles = opts.initialDt * opts.f_cycle;
expectedDamage = opts.fatigueBaseRate * expectedCycles;

assert(abs(fatigue.totalCycles - expectedCycles) < 1e-12, ...
    'Cycles were not committed exactly once.');
assert(max(abs(fatigue.Df_ip(:) - expectedDamage)) < 1e-10, ...
    'Fatigue history was accumulated more than once in one physical step.');
assert(max(abs(fatigue.Df_ip_trial(:) - fatigue.Df_ip(:))) < 1e-12, ...
    'Committed and trial fatigue histories differ after Commit.');
assert(abs(results.Df_max(end) - expectedDamage) < 1e-10, ...
    'Reported fatigue maximum does not match the committed history.');

committed = fatigue.Df_ip;
physics.Assemble(fatigue.Df_Step);
trial1 = fatigue.Df_ip_trial;
physics.Assemble(fatigue.Df_Step);
trial2 = fatigue.Df_ip_trial;
assert(max(abs(trial2(:) - trial1(:))) < 1e-12, ...
    'Repeated assembly changed the fatigue trial state.');
assert(max(abs(fatigue.Df_ip(:) - committed(:))) < 1e-12, ...
    'Assembly modified the committed fatigue history.');

fprintf('PASS test_fatigue_single_commit: Df=%.6g after %.6g cycles.\n', ...
    expectedDamage, expectedCycles);
end
