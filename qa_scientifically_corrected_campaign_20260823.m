function report = qa_scientifically_corrected_campaign_20260823(resultRoot, opts)
%QA_SCIENTIFICALLY_CORRECTED_CAMPAIGN_20260823 Automated 16-case QA.
%
% REPORT = QA_SCIENTIFICALLY_CORRECTED_CAMPAIGN_20260823(RESULTROOT, OPTS)
% validates the completed scientifically corrected campaign without changing
% any simulation result.  It requires the canonical 16-case manifest, loads
% end.mat only (never a numeric checkpoint), checks time/cycle/history
% integrity, and optionally verifies the canonical post-processing CSV.
%
% Optional OPTS fields:
%   OldRoot       - comparison campaign (default Results_Parametric_Corrected)
%   SummaryCsv    - canonical summary CSV to verify
%   RequireCsv    - missing SummaryCsv is a hard failure (default false)
%   FailOnError   - throw if a hard check fails (default true)
%   ReportDir     - if nonempty, write QA tables there (default empty/read-only)
%
% Example, after post-processing:
%   o = struct('RequireCsv',true,'FailOnError',true, ...
%       'SummaryCsv','./Figures_Parametric_ScientificallyCorrected_20260823/parametric_summary.csv');
%   R = qa_scientifically_corrected_campaign_20260823( ...
%       './Results_Parametric_ScientificallyCorrected_20260823', o);

root = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(resultRoot)
    resultRoot = fullfile(root, 'Results_Parametric_ScientificallyCorrected_20260823');
end
if nargin < 2 || isempty(opts), opts = struct(); end

opts.OldRoot = option(opts, 'OldRoot', fullfile(root, 'Results_Parametric_Corrected'));
opts.SummaryCsv = option(opts, 'SummaryCsv', fullfile(root, ...
    'Figures_Parametric_ScientificallyCorrected_20260823', 'parametric_summary.csv'));
opts.RequireCsv = logical(option(opts, 'RequireCsv', false));
opts.FailOnError = logical(option(opts, 'FailOnError', true));
opts.ReportDir = char(option(opts, 'ReportDir', ''));

addpath(root);
addpath(genpath(fullfile(root, 'Models')));
addpath(genpath(fullfile(root, 'Shapes')));
addpath(genpath(fullfile(root, 'PostProcessing')));

manifest = scientifically_corrected_manifest_20260823();
checks = cell(0,7);
summaries = repmat(empty_summary(), height(manifest), 1);

if ~isfolder(resultRoot)
    checks(end+1,:) = check_row('CAMPAIGN','result_root','HARD',false, ...
        resultRoot,'existing directory','Campaign result root is missing.');
else
    dirs = dir(resultRoot);
    actual = string({dirs([dirs.isdir]).name});
    actual = actual(~startsWith(actual,'.'));
    extras = setdiff(actual(:), manifest.Case);
    checks(end+1,:) = check_row('CAMPAIGN','manifest_extra_directories','HARD', ...
        isempty(extras), strjoin(extras,','),'none', ...
        'Unexpected directories can contaminate an automatically discovered summary.');
end

for i = 1:height(manifest)
    endFile = fullfile(resultRoot, manifest.Case(i), 'end.mat');
    if ~isfile(endFile)
        summaries(i) = empty_summary(manifest(i,:));
        checks(end+1,:) = check_row(manifest.Case(i),'end_mat_exists','HARD',false, ...
            endFile,'existing end.mat','A numeric checkpoint is not accepted as completion.');
        continue
    end
    [summaries(i), caseChecks] = validate_case(endFile, manifest(i,:));
    checks = [checks; caseChecks]; %#ok<AGROW>
end

caseSummary = struct2table(summaries);
checkTable = cell2table(checks, 'VariableNames', ...
    {'Case','CheckId','Severity','Status','Value','Limit','Message'});
checkTable.Case = string(checkTable.Case);
checkTable.CheckId = string(checkTable.CheckId);
checkTable.Severity = string(checkTable.Severity);
checkTable.Status = string(checkTable.Status);
checkTable.Value = string(checkTable.Value);
checkTable.Limit = string(checkTable.Limit);
checkTable.Message = string(checkTable.Message);

allEndFiles = all(isfile(fullfile(resultRoot, manifest.Case, 'end.mat')));
canonicalSummary = table();
if allEndFiles
    folders = cellstr(fullfile(resultRoot, manifest.Case));
    canonicalSummary = collect_scenario_metrics(folders, cellstr(manifest.Case));
    if height(canonicalSummary) == height(manifest)
        canonicalSummary.Study = cellstr(manifest.Study);
        canonicalSummary = movevars(canonicalSummary, 'Study', 'Before', 'Label');
        canonicalSummary = canonicalize_submission_summary(canonicalSummary);
        checks2 = check_canonical_csv(canonicalSummary, manifest, opts.SummaryCsv, opts.RequireCsv);
    else
        checks2 = check_row('CAMPAIGN','canonical_summary_height','HARD',false, ...
            height(canonicalSummary),height(manifest), ...
            'collect_scenario_metrics skipped at least one completed case.');
    end
else
    checks2 = check_row('CAMPAIGN','canonical_csv_consistency','REVIEW','NOT_ASSESSED', ...
        'campaign incomplete','16 end.mat files', ...
        'CSV consistency is evaluated only after every canonical end.mat exists.');
end
checkTable = append_checks(checkTable, checks2);

gridConvergence = build_grid_convergence(caseSummary);
l9Effects = build_l9_effects(caseSummary, manifest);
timeStepConvergence = table("NOT_ASSESSED", ...
    "The 16-case campaign contains no time-step variants; use the independent 3-day dt validation.", ...
    'VariableNames', {'Status','Reason'});
oldNew = compare_old_new(resultRoot, opts.OldRoot, manifest);

hardFail = checkTable.Severity == "HARD" & checkTable.Status == "FAIL";
reviewFlag = checkTable.Severity == "REVIEW" & checkTable.Status == "FAIL";
report = struct();
report.Manifest = manifest;
report.CaseSummary = caseSummary;
report.Checks = checkTable;
report.GridConvergence = gridConvergence;
report.TimeStepConvergence = timeStepConvergence;
report.L9Effects = l9Effects;
report.CanonicalSummary = canonicalSummary;
report.OldNew = oldNew;
report.Passed = ~any(hardFail);
report.HardFailureCount = nnz(hardFail);
report.ReviewFlagCount = nnz(reviewFlag);

fprintf('\nScientifically corrected campaign QA: %d/%d end.mat files, %d hard failure(s), %d review flag(s).\n', ...
    nnz(isfile(fullfile(resultRoot, manifest.Case, 'end.mat'))), height(manifest), ...
    report.HardFailureCount, report.ReviewFlagCount);
disp(caseSummary(:, {'Order','Case','TimeDays','Cycles','PhiP99','AgingP99', ...
    'DfP99','DfP999','DfMax','ScreeningIndex','DfIPMin','DfIPMax','HardPass'}));
if ~isempty(gridConvergence), disp(gridConvergence); end
if any(hardFail), disp(checkTable(hardFail,:)); end

if ~isempty(opts.ReportDir)
    if ~isfolder(opts.ReportDir), mkdir(opts.ReportDir); end
    writetable(caseSummary, fullfile(opts.ReportDir, 'qa_case_summary.csv'));
    writetable(checkTable, fullfile(opts.ReportDir, 'qa_checks.csv'));
    writetable(gridConvergence, fullfile(opts.ReportDir, 'qa_grid_convergence.csv'));
    writetable(timeStepConvergence, fullfile(opts.ReportDir, 'qa_time_step_status.csv'));
    writetable(l9Effects, fullfile(opts.ReportDir, 'qa_l9_effects.csv'));
    if ~isempty(oldNew)
        writetable(oldNew, fullfile(opts.ReportDir, 'qa_old_new_differences.csv'));
    end
end

if opts.FailOnError && any(hardFail)
    error('CampaignQA:HardFailure', '%d hard campaign QA check(s) failed.', nnz(hardFail));
end
end

function [summary, rows] = validate_case(endFile, expected)
summary = empty_summary(expected);
rows = cell(0,7);
caseName = expected.Case;
rows(end+1,:) = check_row(caseName,'end_mat_exists','HARD',true,endFile, ...
    'existing end.mat','Final result exists.');

try
    S = load(endFile);
catch ME
    rows(end+1,:) = check_row(caseName,'end_mat_load','HARD',false,ME.identifier, ...
        'load succeeds',ME.message);
    return
end

requiredTop = {'mesh','physics','solver','dt','tvec','results','n_max','tmax'};
missingTop = setdiff(requiredTop, fieldnames(S));
rows(end+1,:) = check_row(caseName,'top_level_variables','HARD',isempty(missingTop), ...
    strjoin(missingTop,','),strjoin(requiredTop,','),'Required saved variables.');
if ~isempty(missingTop), return; end

try
    t = double(S.tvec(:));
    r = S.results;
    p = S.physics;
    mesh = S.mesh;
    n = numel(t);

    rows(end+1,:) = check_row(caseName,'saved_controls','HARD', ...
        close_num(S.dt,30,1e-12) && S.n_max == 260 && close_num(S.tmax,20*86400,1e-12), ...
        sprintf('dt0=%.15g,nMax=%.15g,tmax=%.15g',S.dt,S.n_max,S.tmax), ...
        '30 s, 260, 20 days','Campaign controls match the declared run.');
    rows(end+1,:) = check_row(caseName,'time_finite','HARD',all(isfinite(t)), ...
        nnz(~isfinite(t)),0,'No nonfinite saved time.');
    rows(end+1,:) = check_row(caseName,'time_strictly_increasing','HARD', ...
        n >= 2 && t(1) == 0 && all(diff(t) > 0), ...
        sprintf('n=%d,minDt=%.15g',n,min(diff(t))), ...
        't(1)=0 and diff(t)>0','Time vector is complete and ordered.');
    refTime = expected_time_vector(S.dt,S.tmax,S.n_max,6*3600);
    if numel(refTime) == n
        refErr = max(abs(refTime-t));
    else
        refErr = Inf;
    end
    rows(end+1,:) = check_row(caseName,'time_schedule','HARD', ...
        numel(refTime)==n && refErr <= 1e-8, ...
        sprintf('points=%d,maxError=%.6g',n,refErr), ...
        sprintf('%d points, <=1e-8 s',numel(refTime)), ...
        'Reconstructed adaptive schedule must match the saved vector.');
    finalTimeOK = t(end) > S.tmax && t(end)-S.tmax <= 6*3600+1e-8;
    rows(end+1,:) = check_row(caseName,'time_coverage','HARD',finalTimeOK, ...
        t(end)/86400,'(20,20.25] days', ...
        'main_seal terminates after the first step beyond tmax.');
    rows(end+1,:) = check_row(caseName,'physics_time_alignment','HARD', ...
        close_num(p.time,t(end),1e-10) && close_num(p.dt,t(end)-t(end-1),1e-10), ...
        sprintf('physics.time=%.15g,physics.dt=%.15g',p.time,p.dt), ...
        sprintf('%.15g,%.15g',t(end),t(end)-t(end-1)), ...
        'Saved physics and reporting times agree.');

    g = mesh.getGroupIndex("Internal");
    nElem = size(mesh.Elementgroups{g}.Elems,1);
    nxNodes = numel(unique(mesh.Nodes(:,1)));
    nyNodes = numel(unique(mesh.Nodes(:,2)));
    meshOK = nElem == expected.Nx*expected.Ny && ...
        nxNodes == 2*expected.Nx+1 && nyNodes == 2*expected.Ny+1;
    rows(end+1,:) = check_row(caseName,'mesh_manifest','HARD',meshOK, ...
        sprintf('elements=%d,uniqueNodes=%dx%d',nElem,nxNodes,nyNodes), ...
        sprintf('%d,%dx%d',expected.Nx*expected.Ny,2*expected.Nx+1,2*expected.Ny+1), ...
        'Actual Q9 mesh matches the case name.');
    rows(end+1,:) = check_row(caseName,'mesh_finite','HARD', ...
        all(isfinite(mesh.Nodes(:))) && all(isfinite(mesh.Area(:))) && all(mesh.Area(:)>0), ...
        sprintf('nodes=%d,areaMin=%.6g',size(mesh.Nodes,1),min(mesh.Area(:))), ...
        'finite nodes and positive measures','Mesh geometry is usable.');

    phase = find_model(p,'SealPhaseFieldDamage');
    aging = find_model(p,'AgingDegradation');
    fatigue = find_model(p,'FatigueDamage');
    modelOK = ~isempty(phase) && ~isempty(aging) && ~isempty(fatigue) && ...
        numel(p.models)>=8 && isa(p.models{8},'FatigueDamage');
    rows(end+1,:) = check_row(caseName,'required_models','HARD',modelOK, ...
        strjoin(string(cellfun(@class,p.models,'UniformOutput',false)),','), ...
        'phase, aging, FatigueDamage at model 8','Post-processing model contract.');
    if ~modelOK, return; end

    configOK = close_num(phase.l,expected.LengthScaleMm*1e-3,1e-11) && ...
        close_num(phase.P_max,expected.PressureMPa*1e6,1e-11) && ...
        close_num(aging.T_service,expected.TemperatureC+273.15,1e-11) && ...
        close_num(fatigue.P_max,expected.PressureMPa*1e6,1e-11) && ...
        close_num(fatigue.f_cycle,expected.FrequencyHz,1e-11);
    rows(end+1,:) = check_row(caseName,'model_manifest','HARD',configOK, ...
        sprintf('l=%.6gmm,P=%.6gMPa,T=%.6gC,f=%.6gHz',phase.l*1e3, ...
        fatigue.P_max/1e6,aging.T_service-273.15,fatigue.f_cycle), ...
        sprintf('%.6g,%.6g,%.6g,%.6g',expected.LengthScaleMm, ...
        expected.PressureMPa,expected.TemperatureC,expected.FrequencyHz), ...
        'Saved model parameters match the canonical manifest.');
    calibrationOK = close_num(fatigue.W0,6.0e7,1e-12) && ...
        close_num(fatigue.beta,2.0,1e-12) && ...
        close_num(fatigue.energyRateScale,1.0,1e-12) && ...
        close_num(fatigue.baseRate,1.0e-6,1e-12) && ...
        close_num(fatigue.maxRate,1.0e-5,1e-12) && ...
        close_num(fatigue.pressureExponent,1.2,1e-12) && ...
        close_num(fatigue.P_ref,21e6,1e-12) && ...
        close_num(phase.kmin,5e-3,1e-12) && ...
        close_num(phase.degradeThreshold,0.34,1e-12) && ...
        close_num(phase.w_CL,0.25,1e-12) && close_num(phase.w_aging,0.35,1e-12) && ...
        close_num(phase.w_fatigue,0.24,1e-12) && close_num(phase.w_pressure,0.16,1e-12);
    rows(end+1,:) = check_row(caseName,'scientific_calibration','HARD',calibrationOK, ...
        sprintf('W0=%.6g,beta=%.6g,r0=%.6g,base=%.6g,cap=%.6g,mP=%.6g,Pref=%.6g,kmin=%.6g', ...
        fatigue.W0,fatigue.beta,fatigue.energyRateScale,fatigue.baseRate,fatigue.maxRate, ...
        fatigue.pressureExponent,fatigue.P_ref,phase.kmin), ...
        '6e7,2,1,1e-6,1e-5,1.2,21e6,0.005 and normalized phase weights', ...
        'The saved run uses the declared scientifically corrected calibration.');

    stateFinite = true; stateCommitted = true; stateErr = 0;
    for k = 1:numel(p.StateVec)
        x = p.StateVec{k}; xo = p.StateVec_Old{k};
        stateFinite = stateFinite && all(isfinite(x(:))) && all(isfinite(xo(:)));
        if ~isequal(size(x),size(xo))
            stateCommitted = false; stateErr = Inf;
        elseif ~isempty(x)
            e = max(abs(double(x(:))-double(xo(:))));
            stateErr = max(stateErr,e);
            stateCommitted = stateCommitted && e <= 1e-10*(1+max(abs(double(x(:)))));
        end
    end
    rows(end+1,:) = check_row(caseName,'state_finite','HARD',stateFinite, ...
        double(~stateFinite),0,'StateVec and StateVec_Old contain no NaN/Inf.');
    rows(end+1,:) = check_row(caseName,'state_committed','HARD',stateCommitted, ...
        stateErr,'relative 1e-10','Final StateVec equals its committed copy.');
    histOK = ~isempty(phase.Hist) && all(isfinite(phase.Hist(:))) && min(phase.Hist(:)) >= -1e-12;
    rows(end+1,:) = check_row(caseName,'history_field','HARD',histOK, ...
        sprintf('min=%.6g,max=%.6g',min(phase.Hist(:)),max(phase.Hist(:))), ...
        'finite and >=0','Phase history field is physically admissible.');
    logOK = isprop(S.solver,'convergence_log') && ~isempty(S.solver.convergence_log) && ...
        all(isfinite(S.solver.convergence_log(:)));
    rows(end+1,:) = check_row(caseName,'solver_log_finite','HARD',logOK, ...
        mat2str(size(S.solver.convergence_log)),'nonempty and finite', ...
        'The current solver throws before end.mat on nonconvergence.');

    vectorFields = ["CL_avg","CL_max","LFrac","degradation_connected", ...
        "phi_avg","phi_p99","alpha_p99","Df_p99","screening_index", ...
        "alpha_avg","alpha_max","Df_max","totalCycles","leakage_Q", ...
        "leakage_Q_poiseuille","leakage_kmax","leakage_aperture_max", ...
        "leakage_path_cost","leakage_path_mean_aperture", ...
        "leakage_path_min_aperture","leakage_connected", ...
        "leakage_failed_darcy","leakage_failed_poiseuille", ...
        "contactPressureMean","contactClosureMean","contact_top_mean", ...
        "contact_bottom_mean","contact_right_mean","contact_top_max", ...
        "contact_bottom_max","contact_right_max"];
    scalarFields = ["Q_allow","leakage_time_darcy","leakage_time_poiseuille"];
    requiredResults = [vectorFields scalarFields];
    missingResults = setdiff(requiredResults,string(fieldnames(r)));
    rows(end+1,:) = check_row(caseName,'result_fields','HARD',isempty(missingResults), ...
        strjoin(missingResults,','),strjoin(requiredResults,','),'Required time histories.');
    if ~isempty(missingResults), return; end

    lengthOK = true;
    for f = vectorFields
        lengthOK = lengthOK && numel(r.(f)) == n;
    end
    rows(end+1,:) = check_row(caseName,'result_vector_lengths','HARD',lengthOK, ...
        n,'all vector fields equal numel(tvec)','No truncated/misaligned history.');

    sentinelFields = ["leakage_path_cost","leakage_time_darcy","leakage_time_poiseuille"];
    finiteFields = setdiff(requiredResults,sentinelFields);
    nonfinite = strings(0,1);
    for f = finiteFields
        if any(~isfinite(double(r.(f)(:))))
            nonfinite(end+1,1) = f; %#ok<AGROW>
        end
    end
    rows(end+1,:) = check_row(caseName,'result_finite_whitelist','HARD',isempty(nonfinite), ...
        strjoin(nonfinite,','),'all except declared sentinels', ...
        'Only path-cost +Inf and not-crossed time NaN are allowed.');
    pcost = double(r.leakage_path_cost(:));
    pathCostOK = all(~isnan(pcost)) && all(pcost >= 0) && all(~isinf(pcost) | pcost > 0);
    rows(end+1,:) = check_row(caseName,'leakage_path_cost_sentinel','HARD',pathCostOK, ...
        sprintf('+Inf=%d,finite=%d',nnz(isinf(pcost)&pcost>0),nnz(isfinite(pcost))), ...
        'finite >=0 or +Inf','+Inf denotes no finite initial path; NaN/-Inf are invalid.');

    qAllow = double(r.Q_allow);
    [darcyOK,darcyMsg] = leakage_semantics(r.leakage_Q,r.leakage_failed_darcy, ...
        r.leakage_time_darcy,qAllow,t);
    [poisOK,poisMsg] = leakage_semantics(r.leakage_Q_poiseuille, ...
        r.leakage_failed_poiseuille,r.leakage_time_poiseuille,qAllow,t);
    rows(end+1,:) = check_row(caseName,'leakage_time_darcy_sentinel','HARD',darcyOK, ...
        darcyMsg,'NaN iff threshold is never crossed','Failure time semantics.');
    rows(end+1,:) = check_row(caseName,'leakage_time_poiseuille_sentinel','HARD',poisOK, ...
        poisMsg,'NaN iff threshold is never crossed','Failure time semantics.');

    boundedFields = ["CL_avg","CL_max","phi_avg","phi_p99","alpha_p99", ...
        "Df_p99","screening_index","alpha_avg","alpha_max","Df_max"];
    lo = Inf; hi = -Inf;
    for f = boundedFields
        lo = min(lo,min(double(r.(f)(:))));
        hi = max(hi,max(double(r.(f)(:))));
    end
    rows(end+1,:) = check_row(caseName,'reported_normalized_bounds','HARD', ...
        lo >= -1e-10 && hi <= 1+1e-10,sprintf('[%.6g,%.6g]',lo,hi),'[0,1]', ...
        'Reported physical metrics are bounded.');
    hierarchyOK = all(r.CL_avg(:) <= r.CL_max(:)+1e-10) && ...
        all(r.Df_p99(:) <= r.Df_max(:)+1e-10) && ...
        all(r.phi_avg(:) <= r.phi_p99(:)+1e-10);
    rows(end+1,:) = check_row(caseName,'reported_metric_hierarchy','HARD',hierarchyOK, ...
        'avg/P99/max','ordered','Basic summary statistics are internally consistent.');
    screenExpected = max([r.phi_p99(:),r.alpha_p99(:),r.Df_p99(:)],[],2);
    screenErr = max(abs(double(r.screening_index(:))-double(screenExpected(:))));
    rows(end+1,:) = check_row(caseName,'screening_identity','HARD',screenErr<=1e-12, ...
        screenErr,'<=1e-12','I_S=max(phi_P99,alpha_P99,Df_P99).');

    Df = double(fatigue.Df_ip(:));
    Dft = double(fatigue.Df_ip_trial(:));
    dfFinite = ~isempty(Df) && all(isfinite(Df)) && all(isfinite(Dft));
    dfBounds = dfFinite && min(Df)>=-1e-12 && max(Df)<=1+1e-12;
    rows(end+1,:) = check_row(caseName,'fatigue_ip_bounds','HARD',dfBounds, ...
        sprintf('[%.15g,%.15g]',min(Df),max(Df)),'[0,1]', ...
        'Physical fatigue bounds use integration points, not overshooting nodal projection.');
    trialErr = Inf;
    if isequal(size(fatigue.Df_ip),size(fatigue.Df_ip_trial)) && ~isempty(Df)
        trialErr = max(abs(Df-Dft));
    end
    rows(end+1,:) = check_row(caseName,'fatigue_trial_committed','HARD',trialErr<=1e-12, ...
        trialErr,'<=1e-12','Final trial history was committed exactly once.');
    shapeOK = isequal(size(fatigue.Df_ip),[nElem mesh.ipcount1D^2]);
    rows(end+1,:) = check_row(caseName,'fatigue_ip_shape','HARD',shapeOK, ...
        mat2str(size(fatigue.Df_ip)),sprintf('[%d %d]',nElem,mesh.ipcount1D^2), ...
        'Every integration point is represented.');

    cycles = double(r.totalCycles(:));
    expectedCycles = fatigue.f_cycle*t;
    initialOK = abs(cycles(1))<=1e-14 && abs(r.Df_max(1))<=1e-14 && ...
        abs(r.Df_p99(1))<=1e-14 && abs(r.screening_index(1))<=1e-14;
    rows(end+1,:) = check_row(caseName,'clean_initial_history','HARD',initialOK, ...
        sprintf('N=%.6g,DfMax=%.6g,DfP99=%.6g,I=%.6g',cycles(1), ...
        r.Df_max(1),r.Df_p99(1),r.screening_index(1)), ...
        'all zero','Rejects a contaminated restart/history at t=0.');
    cycleErr = max(abs(cycles-expectedCycles));
    cycleTol = 1e-9*max(1,expectedCycles(end));
    rows(end+1,:) = check_row(caseName,'cycles_equal_frequency_time','HARD',cycleErr<=cycleTol, ...
        cycleErr,sprintf('<=%.6g',cycleTol), ...
        'N=f*t proves one Timedep commit per saved physical step.');
    cycleModelErr = abs(fatigue.totalCycles-cycles(end));
    rows(end+1,:) = check_row(caseName,'cycles_model_result_alignment','HARD', ...
        cycleModelErr<=cycleTol,cycleModelErr,sprintf('<=%.6g',cycleTol), ...
        'Model and reporting histories agree.');
    monoTol = 1e-12;
    fatigueMonotone = all(diff(cycles)>=-monoTol) && ...
        all(diff(r.Df_max(:))>=-monoTol) && all(diff(r.Df_p99(:))>=-monoTol);
    rows(end+1,:) = check_row(caseName,'fatigue_monotonicity','HARD',fatigueMonotone, ...
        sprintf('min dN=%.6g,min dMax=%.6g,min dP99=%.6g',min(diff(cycles)), ...
        min(diff(r.Df_max(:))),min(diff(r.Df_p99(:)))), 'all >=-1e-12', ...
        'Committed fatigue cannot heal.');
    cap = min(1,double(fatigue.maxRate)*cycles);
    capExcess = max(double(r.Df_max(:))-cap);
    incExcess = max(diff(double(r.Df_max(:)))-double(fatigue.maxRate)*diff(cycles));
    capTol = 1e-10;
    rows(end+1,:) = check_row(caseName,'fatigue_per_cycle_cap','HARD', ...
        capExcess<=capTol && incExcess<=capTol, ...
        sprintf('global=%.6g,step=%.6g',capExcess,incExcess),'<=1e-10', ...
        'Detects repeated within-step accumulation despite a capped dD/dN.');
    endMaxErr = max([abs(r.Df_max(end)-max(Df)), ...
        abs(fatigue.Df_max-max(Df))]);
    endP99Err = abs(r.Df_p99(end)-prctile(Df,99));
    rows(end+1,:) = check_row(caseName,'fatigue_result_endpoint','HARD', ...
        max(endMaxErr,endP99Err)<=1e-12, ...
        sprintf('maxErr=%.6g,p99Err=%.6g',endMaxErr,endP99Err),'<=1e-12', ...
        'Reported terminal fatigue is recomputed from committed IP history.');

    nonnegativeFields = ["leakage_Q","leakage_Q_poiseuille","leakage_kmax", ...
        "leakage_aperture_max","leakage_path_mean_aperture", ...
        "leakage_path_min_aperture","contactPressureMean","contact_top_mean", ...
        "contact_bottom_mean","contact_right_mean","contact_top_max", ...
        "contact_bottom_max","contact_right_max"];
    minNonnegative = Inf;
    for f = nonnegativeFields
        minNonnegative = min(minNonnegative,min(double(r.(f)(:))));
    end
    rows(end+1,:) = check_row(caseName,'nonnegative_physical_outputs','HARD', ...
        minNonnegative>=-1e-14,minNonnegative,'>=-1e-14', ...
        'Leakage/contact measures cannot be negative.');
    sealWidth = max(mesh.Nodes(:,1))-min(mesh.Nodes(:,1));
    rows(end+1,:) = check_row(caseName,'degradation_extent_bounds','HARD', ...
        min(r.LFrac(:))>=-1e-12 && max(r.LFrac(:))<=sealWidth+1e-12, ...
        sprintf('[%.6g,%.6g]',min(r.LFrac(:)),max(r.LFrac(:))), ...
        sprintf('[0,%.6g]',sealWidth),'Connected extent stays inside the seal.');

    rows(end+1,:) = check_row(caseName,'aging_monotonicity','REVIEW', ...
        all(diff(r.alpha_p99(:))>=-1e-8),min(diff(r.alpha_p99(:))),'>=-1e-8', ...
        'Scientific trend check, not a completion invariant.');
    rows(end+1,:) = check_row(caseName,'phase_percentile_monotonicity','REVIEW', ...
        all(diff(r.phi_p99(:))>=-1e-8),min(diff(r.phi_p99(:))),'>=-1e-8', ...
        'Large decreases require an irreversibility review.');
    rows(end+1,:) = check_row(caseName,'screening_monotonicity','REVIEW', ...
        all(diff(r.screening_index(:))>=-1e-8),min(diff(r.screening_index(:))),'>=-1e-8', ...
        'The maximum of irreversible components should not materially decrease.');

    summary.TimeDays = t(end)/86400;
    summary.Steps = n-1;
    summary.Cycles = cycles(end);
    summary.CLAvg = r.CL_avg(end);
    summary.CLMax = r.CL_max(end);
    summary.PhiP99 = r.phi_p99(end);
    summary.AgingP99 = r.alpha_p99(end);
    summary.DfP99 = r.Df_p99(end);
    summary.DfP999 = prctile(Df,99.9);
    summary.DfMax = r.Df_max(end);
    summary.ScreeningIndex = r.screening_index(end);
    summary.DfIPMin = min(Df);
    summary.DfIPMax = max(Df);
    summary.DfFracGt05 = mean(Df>0.5);
    summary.EndFile = string(endFile);

catch ME
    rows(end+1,:) = check_row(caseName,'case_validation_exception','HARD',false, ...
        ME.identifier,'no exception',ME.message);
end

hard = string(rows(:,3))=="HARD" & string(rows(:,4))=="FAIL";
summary.HardPass = ~any(hard);
end

function rows = check_canonical_csv(T, manifest, summaryCsv, requireCsv)
rows = cell(0,7);
if ~isfile(summaryCsv)
    if requireCsv
        rows(end+1,:) = check_row('CAMPAIGN','canonical_csv_exists','HARD',false, ...
            summaryCsv,'existing CSV','Post-processing CSV is required for final artifact QA.');
    else
        rows(end+1,:) = check_row('CAMPAIGN','canonical_csv_consistency','REVIEW', ...
            'NOT_ASSESSED',summaryCsv,'optional until post-processing', ...
            'Re-run QA with RequireCsv=true after figures/tables are generated.');
    end
    return
end

C = readtable(summaryCsv,'TextType','string');
requiredNames = ["ScreeningIndex","PhaseP99","FatigueP99","AgingP99", ...
    "ConnectedDegradationExtentMm","TimeScreening07Days","TimeScreening1Days"];
legacyNames = ["FailureIndex","RobustFailureIndexP99","CrackNorm","CrackLengthMm", ...
    "FatigueNorm","AgingNorm","FluidNorm","TimeFI07Days","TimeFI1Days", ...
    "TimeCrackCritDays","FatigueP99IP"];
schemaOK = all(ismember(requiredNames,string(C.Properties.VariableNames))) && ...
    ~any(ismember(legacyNames,string(C.Properties.VariableNames)));
rows(end+1,:) = check_row('CAMPAIGN','canonical_csv_schema','HARD',schemaOK, ...
    strjoin(string(C.Properties.VariableNames),','), ...
    'explicit phase/screening names and no crack/failure aliases', ...
    'Submission-facing CSV uses current physical terminology.');
orderOK = height(C)==height(manifest) && ismember('Label',C.Properties.VariableNames) && ...
    isequal(string(C.Label),manifest.Case);
rows(end+1,:) = check_row('CAMPAIGN','canonical_csv_order','HARD',orderOK, ...
    strjoin(string(C.Label),','),strjoin(manifest.Case,','), ...
    'CSV rows must exactly follow the 16-case manifest.');
if ~orderOK, return; end

common = intersect(string(T.Properties.VariableNames),string(C.Properties.VariableNames),'stable');
ignore = ["Folder","File"];
common = setdiff(common,ignore,'stable');
bad = strings(0,1); maxErr = 0;
for name = common
    a = T.(name); b = C.(name);
    if isnumeric(a) || islogical(a)
        a = double(a); b = double(b);
        if ~isequal(size(a),size(b)) || ~isequal(isnan(a),isnan(b))
            bad(end+1,1)=name; %#ok<AGROW>
            continue
        end
        finite = isfinite(a) & isfinite(b);
        if any(finite(:))
            e = max(abs(a(finite)-b(finite)));
            maxErr = max(maxErr,e);
            if e > 1e-10*(1+max(abs(a(finite))))
                bad(end+1,1)=name; %#ok<AGROW>
            end
        end
        if ~isequal(isinf(a),isinf(b))
            bad(end+1,1)=name; %#ok<AGROW>
        end
    else
        if ~isequal(string(a),string(b))
            bad(end+1,1)=name; %#ok<AGROW>
        end
    end
end
rows(end+1,:) = check_row('CAMPAIGN','canonical_csv_values','HARD',isempty(bad), ...
    sprintf('bad=%s,maxError=%.6g',strjoin(unique(bad),','),maxErr), ...
    'MAT-derived canonical table within relative 1e-10', ...
    'Tables and endpoint figures must share the same canonical CSV values.');
end

function G = build_grid_convergence(S)
G = table();
if height(S)<3 || any(~S.HardPass(1:3)), return; end
metrics = ["ScreeningIndex","DfP99","DfP999","DfMax","PhiP99","AgingP99","CLAvg"];
rows = cell(numel(metrics),12);
h = [0.193373117;0.125;0.088388348];
for k = 1:numel(metrics)
    y = S.(metrics(k))(1:3);
    scale = max(abs(y(3)),0.01);
    eCM = abs(y(2)-y(1))/scale;
    eMF = abs(y(3)-y(2))/scale;
    contraction = eMF/max(eCM,eps);
    if eMF<=0.01 && contraction<=1
        status = "PASS";
    elseif eMF<=0.02
        status = "REVIEW";
    else
        status = "FAIL_CLAIM";
    end
    coarseSensitive = eCM>0.05 && eMF<=0.01;
    if coarseSensitive
        note = "medium/fine converged; coarse grid is tail-sensitive";
    elseif status=="PASS"
        note = "successive fine-grid change is below 1% and contracts";
    else
        note = "do not claim convergence without review";
    end
    rows(k,:) = {metrics(k),h(1),h(2),h(3),y(1),y(2),y(3),eCM,eMF, ...
        contraction,coarseSensitive,status+": "+note};
end
G = cell2table(rows,'VariableNames',{'Metric','hCoarseMm','hMediumMm','hFineMm', ...
    'Coarse','Medium','Fine','RelativeCoarseMedium','RelativeMediumFine', ...
    'ErrorContractionRatio','CoarseSensitive','Status'});
end

function E = build_l9_effects(S,M)
E = table();
idx = find(M.Study=="orthogonal");
if numel(idx)~=9 || any(~S.HardPass(idx)), return; end
metrics = ["ScreeningIndex","DfP99","AgingP99"];
factors = ["PressureMPa","TemperatureC","FrequencyHz"];
rows = cell(0,8);
for metric = metrics
    y = S.(metric)(idx);
    ranges = zeros(1,3);
    means = cell(1,3);
    for j = 1:3
        x = M.(factors(j))(idx);
        lev = unique(x,'sorted'); mu = zeros(1,3);
        for q = 1:3, mu(q)=mean(y(x==lev(q))); end
        means{j}=mu; ranges(j)=range(mu);
    end
    ranking = rank_effects(factors,ranges,y);
    for j = 1:3
        mu=means{j};
        rows(end+1,:)={metric,factors(j),mu(1),mu(2),mu(3),ranges(j), ...
            all(diff(mu)>=-1e-10),ranking}; %#ok<AGROW>
    end
end
E=cell2table(rows,'VariableNames',{'Response','Factor','Level1Mean','Level2Mean', ...
    'Level3Mean','Range','Nondecreasing','RangeRanking'});
end

function ranking = rank_effects(labels,ranges,response)
% Treat round-off-sized effects as ties rather than inventing a rank order.
tol = max(1e-12,1e-9*max(abs(response),[],'omitnan'));
clean = ranges;
clean(abs(clean)<=tol) = 0;
[sorted,ord] = sort(clean,'descend');
ranking = labels(ord(1));
for k = 2:numel(ord)
    if abs(sorted(k)-sorted(k-1))<=tol
        separator = " = ";
    else
        separator = " > ";
    end
    ranking = ranking + separator + labels(ord(k));
end
end

function D = compare_old_new(newRoot,oldRoot,M)
D = table();
if ~isfolder(oldRoot), return; end
metrics = ["CL_avg","CL_max","phi_p99","alpha_p99","Df_p99", ...
    "Df_max","screening_index","LFrac","leakage_Q","leakage_Q_poiseuille"];
rows = cell(0,9);
for i=1:height(M)
    nf=fullfile(newRoot,M.Case(i),'end.mat'); of=fullfile(oldRoot,M.Case(i),'end.mat');
    if ~isfile(nf) || ~isfile(of), continue; end
    a=load(of,'results','tvec'); b=load(nf,'results','tvec');
    for metric=metrics
        if ~isfield(a.results,metric) || ~isfield(b.results,metric), continue; end
        x=double(a.results.(metric)(:)); y=double(b.results.(metric)(:));
        if isequal(a.tvec,b.tvec) && numel(x)==numel(y)
            yy=y; xx=x;
        else
            tCommon=double(b.tvec(:));
            xx=interp1(double(a.tvec(:)),x,tCommon,'linear',NaN); yy=y;
        end
        valid=isfinite(xx)&isfinite(yy);
        if any(valid)
            rmse=sqrt(mean((yy(valid)-xx(valid)).^2));
            maxAbs=max(abs(yy(valid)-xx(valid)));
        else
            rmse=NaN; maxAbs=NaN;
        end
        old=x(end); new=y(end); delta=new-old;
        rel=delta/max(abs(old),metric_scale(metric));
        rows(end+1,:)={M.Case(i),metric,old,new,delta,rel,rmse,maxAbs, ...
            isequal(a.tvec,b.tvec)}; %#ok<AGROW>
    end
end
if ~isempty(rows)
    D=cell2table(rows,'VariableNames',{'Case','Metric','OldFinal','NewFinal', ...
        'Delta','RelativeDelta','CurveRMSE','CurveMaxAbs','SameTimeVector'});
end
end

function v = metric_scale(metric)
if startsWith(metric,"leakage_Q"), v=1e-30; else, v=0.01; end
end

function [ok,msg] = leakage_semantics(q,flag,stored,qAllow,t)
q=double(q(:)); flag=logical(flag(:));
expected=cummax(double(q>=qAllow))>0;
idx=find(expected,1,'first');
if isempty(idx)
    ok=~any(flag) && isscalar(stored) && isnan(stored);
    msg=sprintf('crossing=none,flag=%d,time=%g',any(flag),stored);
else
    expectedTime=t(idx)/86400;
    ok=isequal(flag,expected) && isscalar(stored) && isfinite(stored) && ...
        abs(double(stored)-expectedTime)<=1e-10*(1+expectedTime);
    msg=sprintf('crossing=%d,expected=%.15g,time=%.15g',idx,expectedTime,stored);
end
end

function t = expected_time_vector(dt0,tmax,nmax,maxDt)
t=0;
for step=1:nmax
    if step<=30
        dtk=dt0;
    else
        dtk=dt0*1.04^(min(step-30,200));
    end
    dtk=min(dtk,maxDt);
    t(end+1,1)=t(end)+dtk; %#ok<AGROW>
    if t(end)>tmax, break; end
end
end

function mdl = find_model(physics,className)
mdl=[];
for k=1:numel(physics.models)
    if isa(physics.models{k},className)
        mdl=physics.models{k}; return
    end
end
end

function tf = close_num(a,b,relTol)
tf=isscalar(a)&&isscalar(b)&&isfinite(a)&&isfinite(b)&& ...
    abs(double(a)-double(b))<=relTol*(1+abs(double(b)));
end

function s = empty_summary(row)
if nargin<1 || isempty(row)
    order=NaN; caseName=""; study="";
else
    order=row.Order; caseName=row.Case; study=row.Study;
end
s=struct('Order',order,'Case',caseName,'Study',study,'EndFile',"", ...
    'TimeDays',NaN,'Steps',NaN,'Cycles',NaN,'CLAvg',NaN,'CLMax',NaN, ...
    'PhiP99',NaN,'AgingP99',NaN,'DfP99',NaN,'DfP999',NaN,'DfMax',NaN, ...
    'ScreeningIndex',NaN,'DfIPMin',NaN,'DfIPMax',NaN,'DfFracGt05',NaN, ...
    'HardPass',false);
end

function row = check_row(caseName,id,severity,passed,value,limit,message)
if islogical(passed)
    if passed, status='PASS'; else, status='FAIL'; end
else
    status=char(string(passed));
end
row={char(string(caseName)),char(string(id)),char(string(severity)),status, ...
    value_text(value),value_text(limit),char(string(message))};
end

function out = value_text(value)
if isstring(value) || ischar(value)
    out=char(strjoin(string(value),','));
elseif isnumeric(value) || islogical(value)
    if isscalar(value), out=sprintf('%.17g',double(value)); else, out=mat2str(value); end
else
    out=char(string(value));
end
end

function T = append_checks(T,rows)
if isempty(rows), return; end
A=cell2table(rows,'VariableNames',T.Properties.VariableNames);
for name=string(T.Properties.VariableNames)
    A.(name)=string(A.(name));
end
T=[T;A];
end

function value = option(s,name,defaultValue)
if isfield(s,name) && ~isempty(s.(name)), value=s.(name); else, value=defaultValue; end
end
