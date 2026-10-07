function run_one_corrected_case(caseName, resultRoot, maxThreads)
    %RUN_ONE_CORRECTED_CASE Run one submission-calibration case.
    %
    % Results are written to a new root so the archived calibration remains
    % untouched.  The corrected configuration uses constant aging-agent
    % diffusivity, normalized degradation-drive weights, a connected extent
    % threshold of phi = 0.35, and percentile-based screening outputs.

    addpath(genpath('./Models'))
    addpath(genpath('./Shapes'))
    addpath(genpath('./PostProcessing'))

    if nargin < 2 || isempty(resultRoot)
        resultRoot = './Results_Parametric_Corrected';
    end
    if nargin < 3 || isempty(maxThreads)
        maxThreads = 1;
    end
    opts = calibrated_degradation_options(resultRoot);
    opts.tmaxDays = 20;
    opts.nMax = 260;
    opts.plotEvery = inf;
    opts.saveEvery = inf;
    opts.verboseEvery = inf;
    opts.useParpool = false;
    opts.cleanResults = true;
    opts.forceRestart = false;
    opts.maxThreads = maxThreads;
    opts.degradationExtentThreshold = 0.35;
    opts.sname = char(caseName);

    [kind, values] = parse_case_name(char(caseName));
    switch kind
        case 'grid'
            opts.mesh = struct('Nx', values(1), 'Ny', values(2));
            opts.l = 0.30e-3;
        case 'length_scale'
            opts.mesh = struct('Nx', 28, 'Ny', 40);
            h = sqrt((3.5e-3/opts.mesh.Nx) * (5.0e-3/opts.mesh.Ny));
            opts.l = values(1) * h;
        case 'orthogonal'
            opts.mesh = struct('Nx', 24, 'Ny', 34);
            opts.l = 0.30e-3;
            opts.P_max = values(1) * 1e6;
            opts.T_service = values(2) + 273.15;
            opts.f_cycle = values(3) / 1000;
        otherwise
            error('Unknown corrected-study case name: %s', caseName);
    end

    main_seal(1, opts);
end

function [kind, values] = parse_case_name(name)
    tok = regexp(name, '^Grid_G(\d+)x(\d+)$', 'tokens', 'once');
    if ~isempty(tok)
        kind = 'grid';
        values = [str2double(tok{1}), str2double(tok{2})];
        return
    end

    tok = regexp(name, '^LengthScale_lh([\d.]+)$', 'tokens', 'once');
    if ~isempty(tok)
        kind = 'length_scale';
        values = str2double(tok{1});
        return
    end

    tok = regexp(name, '^Orthogonal_O\d+_P(\d+)_T(\d+)_f(\d+)$', 'tokens', 'once');
    if ~isempty(tok)
        kind = 'orthogonal';
        values = [str2double(tok{1}), str2double(tok{2}), str2double(tok{3})];
        return
    end

    kind = 'unknown';
    values = [];
end
