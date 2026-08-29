function opts = calibrated_degradation_options(resultRoot)
    %CALIBRATED_DEGRADATION_OPTIONS Prescribed paper-study parameter set.
    %
    % The historical function name is retained for interface compatibility.
    % These values are fixed/tuned numerical inputs for the comparative
    % campaign; they were not fitted to experimental seal-life or leakage
    % data. The parameterization introduces a thresholded degradation drive
    % and a pressure-assisted term to avoid immediate CL-driven saturation.

    if nargin < 1 || isempty(resultRoot)
        resultRoot = './Results_Parametric_Calibrated';
    end

    opts = struct();
    opts.saveRoot = resultRoot;
    opts.plotEvery = inf;
    opts.saveEvery = inf;
    opts.verboseEvery = inf;
    opts.useParpool = false;
    opts.cleanResults = true;
    opts.forceRestart = false;
    opts.maxThreads = 1;

    opts.CL_boundary = 0.65;
    opts.alpha_boundary = 0.05;
    opts.fluidDegrade = 0.25;
    opts.agingDegrade = 0.45;
    opts.fatigueDegrade = 0.55;

    opts.alpha_coupling = 1.2;
    opts.degradeThreshold = 0.34;
    opts.w_CL = 0.25;
    opts.w_aging = 0.35;
    opts.w_fatigue = 0.24;
    opts.w_pressure = 0.16;
    opts.P_ref = 35e6;

    % Keep aging-agent transport independent of phi in the calibrated
    % submission model.  The previous hard-coded 1000*phi multiplier is
    % available only through an explicit sensitivity override.
    opts.agingPhiDiffusionFactor = 0.0;

    opts.W0 = 6.0e7;
    opts.beta = 2.0;
    opts.fatiguePressureExponent = 1.2;
    opts.fatigueP_ref = 21e6;
    opts.fatigueEnergyRate = 1.0;
    opts.fatigueBaseRate = 1.0e-6;
    opts.fatigueMaxRate = 1.0e-5;
end
