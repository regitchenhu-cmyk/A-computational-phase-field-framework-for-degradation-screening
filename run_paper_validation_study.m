function validation = run_paper_validation_study(overrides)
    %RUN_PAPER_VALIDATION_STUDY Grid and time-step checks for manuscript use.
    %
    % The middle-mesh baseline may be loaded from an existing optimized
    % contact-aware run; coarse/fine and time-step variants are run fresh.

    if nargin < 1 || isempty(overrides)
        overrides = struct();
    end

    sourceRoot = fileparts(mfilename('fullpath'));
    addpath(sourceRoot);
    addpath(genpath(fullfile(sourceRoot, 'Models')));
    addpath(genpath(fullfile(sourceRoot, 'Shapes')));

    root = string(get_option(overrides, 'saveRoot', './Results_Paper_OptimizedRuns/Validation'));
    if ~exist(root, 'dir'), mkdir(root); end

    baseFile = string(get_option(overrides, 'baseFile', ...
        './Results_Paper_OptimizedRuns/FixedVsContact/Standard_21MPa_80C_contactAware/end.mat'));
    scenario = get_option(overrides, 'scenario', 1);
    Q_allow = get_option(overrides, 'Q_allow', 1e-12);

    common = optimized_options(overrides, root);

    gridRows = [];
    fprintf('\n=== Validation: grid sensitivity ===\n');
    gridRows = [gridRows; run_or_load_case("G30x42", "grid", scenario, common, ...
        struct('mesh', struct('Nx', 30, 'Ny', 42), 'nameSuffix', "grid_30x42"), Q_allow)]; %#ok<AGROW>
    gridRows = [gridRows; load_case("G40x56", "grid", baseFile, Q_allow)]; %#ok<AGROW>
    gridRows = [gridRows; run_or_load_case("G50x70", "grid", scenario, common, ...
        struct('mesh', struct('Nx', 50, 'Ny', 70), 'nameSuffix', "grid_50x70"), Q_allow)]; %#ok<AGROW>
    gridSummary = struct2table(gridRows);
    gridSummary = add_relative_columns(gridSummary);
    writetable(gridSummary, fullfile(root, "grid_sensitivity_summary.csv"));

    timeRows = [];
    fprintf('\n=== Validation: time-step sensitivity ===\n');
    timeRows = [timeRows; run_or_load_case("dt20s", "time", scenario, common, ...
        struct('initialDt', 20, 'nMax', 215, 'nameSuffix', "dt_20s"), Q_allow)]; %#ok<AGROW>
    timeRows = [timeRows; load_case("dt30s_baseline", "time", baseFile, Q_allow)]; %#ok<AGROW>
    timeRows = [timeRows; run_or_load_case("dt45s", "time", scenario, common, ...
        struct('initialDt', 45, 'nMax', 170, 'nameSuffix', "dt_45s"), Q_allow)]; %#ok<AGROW>
    timeSummary = struct2table(timeRows);
    timeSummary = add_relative_columns(timeSummary);
    writetable(timeSummary, fullfile(root, "time_step_sensitivity_summary.csv"));

    fig = plot_validation(gridSummary, timeSummary);
    figFile = fullfile(root, "validation_summary.png");
    exportgraphics(fig, figFile, 'Resolution', 300);
    close(fig);

    validation.root = root;
    validation.gridSummary = gridSummary;
    validation.timeSummary = timeSummary;
    validation.figure = figFile;

    disp(gridSummary)
    disp(timeSummary)
    fprintf('Saved validation results to %s\n', root);
end

function common = optimized_options(overrides, root)
    common = overrides;
    common.saveRoot = root + "/cases";
    common.tmaxDays = get_option(common, 'tmaxDays', 3);
    common.nMax = get_option(common, 'nMax', 180);
    common.saveEvery = get_option(common, 'saveEvery', inf);
    common.plotEvery = get_option(common, 'plotEvery', inf);
    common.verboseEvery = get_option(common, 'verboseEvery', inf);
    common.useParpool = get_option(common, 'useParpool', true);
    common.cleanResults = get_option(common, 'cleanResults', true);
    common.useContactPenalty = true;
    common.alpha_coupling = get_option(common, 'alpha_coupling', 0.8);
    common.contactFatigueCoeff = get_option(common, 'contactFatigueCoeff', 0.10);
    common.leakage_k_phi = get_option(common, 'leakage_k_phi', 8e-18);
    common.leakagePhiExponent = get_option(common, 'leakagePhiExponent', 4.0);
    common.leakageAperturePhi = get_option(common, 'leakageAperturePhi', 1.2e-5);
    common.leakageApertureExponent = get_option(common, 'leakageApertureExponent', 3.0);
    common.leakageContactChi = get_option(common, 'leakageContactChi', 5.5);
end

function row = run_or_load_case(label, group, scenario, common, caseOverrides, Q_allow)
    opts = merge_struct(common, caseOverrides);
    endFile = expected_end_file(scenario, opts);
    if isfile(endFile)
        fprintf('\n--- Loading completed %s (%s) from %s ---\n', label, group, endFile);
        row = load_case(label, group, endFile, Q_allow);
        return;
    end

    lockFile = expected_lock_file(opts);
    if isfile(lockFile)
        fprintf('\n--- Waiting for running %s (%s): %s ---\n', label, group, lockFile);
        waitStart = tic;
        while isfile(lockFile) && ~isfile(endFile)
            pause(5);
            if toc(waitStart) > 3 * 3600
                error('Validation:WorkerTimeout', ...
                    'Timed out waiting for validation worker: %s', lockFile);
            end
        end
        if isfile(endFile)
            row = load_case(label, group, endFile, Q_allow);
            return;
        end
        error('Validation:WorkerFailed', ...
            'Validation worker stopped without writing %s', endFile);
    end

    fprintf('\n--- Running %s (%s) ---\n', label, group);
    [physics, tvec, results] = main_seal(scenario, opts);
    row = extract_metrics(label, group, physics, tvec, results, Q_allow);
end

function endFile = expected_end_file(scenario, opts)
if isfield(opts, 'sname') && ~isempty(opts.sname)
    baseName = string(opts.sname);
else
    switch scenario
        case 1
            baseName = "Standard_21MPa_80C";
        case 2
            baseName = "HighTemp_21MPa_120C";
        case 3
            baseName = "HighPress_35MPa_80C";
        otherwise
            error('Validation:UnknownScenario', 'Unknown scenario: %g', scenario);
    end
end
endFile = fullfile(string(opts.saveRoot), ...
    baseName + "_" + string(opts.nameSuffix), 'end.mat');
end

function lockFile = expected_lock_file(opts)
caseRoot = string(opts.saveRoot);
validationRoot = fileparts(caseRoot);
lockFile = fullfile(validationRoot, 'locks', ...
    string(opts.nameSuffix) + ".running");
end

function row = load_case(label, group, filePath, Q_allow)
    fprintf('\n--- Loading %s from %s ---\n', label, filePath);
    data = load(filePath, "physics", "tvec", "results");
    row = extract_metrics(label, group, data.physics, data.tvec, data.results, Q_allow);
end

function row = extract_metrics(label, group, physics, tvec, results, Q_allow)
    leakModel = find_model(physics, "LeakagePostProcessor");
    if ~isempty(leakModel)
        leakModel.Evaluate(physics);
    end

    row.Case = string(label);
    row.Group = string(group);
    row.FinalTime_days = physics.time / 86400;
    row.PhiMax = get_phi_max(physics);
    row.DfMax = safe_last(results, 'Df_max');
    row.LFrac_mm = safe_last(results, 'LFrac') * 1e3;
    row.DarcyQ_m3s = safe_metric(leakModel, 'Q', safe_last(results, 'leakage_Q'));
    row.PoiseuilleQ_m3s = safe_metric(leakModel, 'Q_poiseuille', safe_last(results, 'leakage_Q_poiseuille'));
    row.PathCost = safe_metric(leakModel, 'pathCost', safe_last(results, 'leakage_path_cost'));
    row.PathMinAperture_um = safe_metric(leakModel, 'pathMinAperture', safe_last(results, 'leakage_path_min_aperture')) * 1e6;
    row.TopPc_MPa = safe_last(results, 'contact_top_mean') / 1e6;
    evalTime = 3 * 86400;
    row.EvalTime_days = evalTime / 86400;
    row.DfAt3d = value_at_time(tvec, results, 'Df_max', evalTime);
    row.DarcyQAt3d_m3s = value_at_time(tvec, results, 'leakage_Q', evalTime);
    row.PoiseuilleQAt3d_m3s = value_at_time(tvec, results, ...
        'leakage_Q_poiseuille', evalTime);
    row.PhaseP99At3d = value_at_time(tvec, results, 'phi_p99', evalTime);
    row.ScreeningIndexAt3d = value_at_time(tvec, results, ...
        'screening_index', evalTime);
    row.TopPcAt3d_MPa = value_at_time(tvec, results, ...
        'contact_top_mean', evalTime) / 1e6;
    row.DfRelToBaseline_pct = NaN;
    row.QpRelToBaseline_pct = NaN;
    row.QAllow_m3s = Q_allow;

    if ~isempty(tvec)
        row.TvecEnd_days = tvec(end) / 86400;
    else
        row.TvecEnd_days = row.FinalTime_days;
    end
end

function fig = plot_validation(gridSummary, timeSummary)
    fig = figure('Color', 'w', 'Position', [80 80 1150 720]);
    tiledlayout(2,2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile
    bar(categorical(gridSummary.Case), gridSummary.DfAt3d);
    ylabel('D_f max at 3 d [-]'); title('Grid sensitivity: fatigue'); grid on

    nexttile
    semilogy(categorical(gridSummary.Case), max(gridSummary.PoiseuilleQAt3d_m3s, realmin), 'o-', 'LineWidth', 1.5);
    ylabel('Poiseuille Q at 3 d [m^3/s]'); title('Grid sensitivity: leakage'); grid on

    nexttile
    bar(categorical(timeSummary.Case), timeSummary.DfAt3d);
    ylabel('D_f max at 3 d [-]'); title('Time-step sensitivity: fatigue'); grid on

    nexttile
    semilogy(categorical(timeSummary.Case), max(timeSummary.PoiseuilleQAt3d_m3s, realmin), 'o-', 'LineWidth', 1.5);
    ylabel('Poiseuille Q at 3 d [m^3/s]'); title('Time-step sensitivity: leakage'); grid on
end

function T = add_relative_columns(T)
    idx = find(contains(T.Case, "baseline") | contains(T.Case, "G40x56"), 1, 'first');
    if isempty(idx), idx = ceil(height(T)/2); end
    baseDf = T.DfAt3d(idx);
    baseQp = T.PoiseuilleQAt3d_m3s(idx);
    T.DfRelToBaseline_pct = 100 * (T.DfAt3d - baseDf) / max(abs(baseDf), eps);
    T.QpRelToBaseline_pct = 100 * (T.PoiseuilleQAt3d_m3s - baseQp) / max(abs(baseQp), eps);
end

function value = value_at_time(tvec, results, fieldName, targetTime)
if isempty(tvec) || tvec(1) > targetTime || tvec(end) < targetTime
    error('Validation:EvaluationTime', ...
        'Saved history does not bracket %.15g s.', targetTime);
end
if ~isfield(results, fieldName) || numel(results.(fieldName)) ~= numel(tvec)
    error('Validation:EvaluationField', ...
        'History %s is absent or does not match tvec.', fieldName);
end
value = interp1(double(tvec(:)), double(results.(fieldName)(:)), ...
    targetTime, 'linear');
if ~isscalar(value) || ~isfinite(value)
    error('Validation:EvaluationValue', ...
        'Interpolated %s at %.15g s is invalid.', fieldName, targetTime);
end
end

function out = merge_struct(a, b)
    out = a;
    names = fieldnames(b);
    for k = 1:numel(names)
        out.(names{k}) = b.(names{k});
    end
end

function mdl = find_model(physics, modelName)
    mdl = [];
    for m = 1:length(physics.models)
        if isprop(physics.models{m}, 'myName') && physics.models{m}.myName == modelName
            mdl = physics.models{m};
            return;
        end
    end
end

function value = get_phi_max(physics)
    [phiType, phiStep] = physics.dofSpace.getDofType({"phi"});
    allNodes = physics.mesh.GetAllNodesForGroup(1);
    phiDofs = physics.dofSpace.getDofIndices(phiType, allNodes);
    value = max(physics.StateVec{phiStep}(phiDofs));
end

function value = safe_metric(model, propName, fallback)
    value = fallback;
    if ~isempty(model) && isprop(model, propName) && ~isempty(model.(propName))
        value = model.(propName);
    end
end

function value = safe_last(results, fieldName)
    if isfield(results, fieldName)
        s = results.(fieldName);
        if isempty(s)
            value = 0;
        elseif isscalar(s)
            value = s;
        else
            value = s(end);
        end
    else
        value = 0;
    end
end

function value = get_option(s, name, defaultValue)
    if isfield(s, name) && ~isempty(s.(name))
        value = s.(name);
    else
        value = defaultValue;
    end
end
