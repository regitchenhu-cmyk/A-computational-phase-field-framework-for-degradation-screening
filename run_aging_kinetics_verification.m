function T = run_aging_kinetics_verification(outdir)
    %RUN_AGING_KINETICS_VERIFICATION Verify the n=1 reaction update.
    %
    % For a spatially uniform state without diffusion, the calibrated aging
    % law has the exact solution alpha=1-(1-alpha0)exp(-k_a*t).  This script
    % compares that solution with backward Euler at the time-step sizes used
    % for long-time integration.

    if nargin < 1 || isempty(outdir)
        outdir = './Figures_Parametric_Corrected';
    end
    if ~exist(outdir, 'dir'), mkdir(outdir); end

    temperaturesC = [80; 100; 120];
    dts = [3600; 6*3600];
    alpha0 = 0.05;
    tEnd = 20*86400;
    kAge = 2e-9;
    Ea = 50e3;
    R = 8.31446261815324;
    Tref = 293.15;

    rows = {};
    for i = 1:numel(temperaturesC)
        tempK = temperaturesC(i) + 273.15;
        kA = kAge * exp(-Ea/R * (1/tempK - 1/Tref));
        exact = 1 - (1-alpha0)*exp(-kA*tEnd);
        for j = 1:numel(dts)
            n = ceil(tEnd/dts(j));
            dt = tEnd/n;
            numeric = 1 - (1-alpha0)/(1+kA*dt)^n;
            relError = abs(numeric-exact)/max(abs(exact), eps);
            rows(end+1,:) = {temperaturesC(i), dt, exact, numeric, relError}; %#ok<AGROW>
        end
    end

    T = cell2table(rows, 'VariableNames', ...
        {'Temperature_C','TimeStep_s','AlphaExact','AlphaBackwardEuler','RelativeError'});
    writetable(T, fullfile(outdir, 'aging_kinetics_verification.csv'));
    disp(T)
end
