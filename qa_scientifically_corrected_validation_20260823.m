function report = qa_scientifically_corrected_validation_20260823()
%QA_SCIENTIFICALLY_CORRECTED_VALIDATION_20260823 Audit six final states.
%
% The check is read-only with respect to solver states. It verifies the
% fixed/contact pair and the four independently run grid/time-step variants,
% then compares the two published validation CSV files with values re-derived
% from end.mat.

root = fileparts(mfilename('fullpath'));
paperRoot = fullfile(root, 'Results_Paper_ScientificallyCorrected_20260823');
validationRoot = fullfile(paperRoot, 'Validation');
qaRoot = fullfile(validationRoot, 'qa');
if ~isfolder(qaRoot), mkdir(qaRoot); end

addpath(root);
addpath(genpath(fullfile(root, 'Models')));
addpath(genpath(fullfile(root, 'Shapes')));

caseSpec = table( ...
    ["Fixed BC"; "Contact-aware"; "G30x42"; "G50x70"; "dt20s"; "dt45s"], ...
    ["fixed"; "contact"; "grid"; "grid"; "time"; "time"], ...
    [40;40;30;50;40;40], [56;56;42;70;56;56], ...
    [30;30;30;30;20;45], [180;180;180;180;215;170], ...
    [string(fullfile(paperRoot,'FixedVsContact','Standard_21MPa_80C_fixedBC','end.mat')); ...
     string(fullfile(paperRoot,'FixedVsContact','Standard_21MPa_80C_contactAware','end.mat')); ...
     string(fullfile(validationRoot,'cases','Standard_21MPa_80C_grid_30x42','end.mat')); ...
     string(fullfile(validationRoot,'cases','Standard_21MPa_80C_grid_50x70','end.mat')); ...
     string(fullfile(validationRoot,'cases','Standard_21MPa_80C_dt_20s','end.mat')); ...
     string(fullfile(validationRoot,'cases','Standard_21MPa_80C_dt_45s','end.mat'))], ...
    'VariableNames', {'Case','Group','Nx','Ny','InitialDt_s','NMax','EndFile'});

checks = cell(0,7);
metrics = repmat(empty_metrics(),height(caseSpec),1);
for i = 1:height(caseSpec)
    [metrics(i),caseChecks] = validate_case(caseSpec(i,:));
    checks = [checks; caseChecks]; %#ok<AGROW>
end
caseTable = struct2table(metrics);

gridCsv = fullfile(validationRoot,'grid_sensitivity_summary.csv');
timeCsv = fullfile(validationRoot,'time_step_sensitivity_summary.csv');
fixedCsv = fullfile(paperRoot,'fixed_vs_contact_summary.csv');
checks = [checks; check_summary_files(caseTable,gridCsv,timeCsv,fixedCsv)]; %#ok<AGROW>

checkTable = cell2table(checks,'VariableNames', ...
    {'Case','CheckId','Severity','Status','Value','Limit','Message'});
for name = ["Case","CheckId","Severity","Status","Value","Limit","Message"]
    checkTable.(name) = string(checkTable.(name));
end

sensitivity = build_sensitivity(caseTable);
hardFail = checkTable.Severity=="HARD" & checkTable.Status=="FAIL";
report = struct('CaseSummary',caseTable,'Checks',checkTable, ...
    'Sensitivity',sensitivity,'Passed',~any(hardFail), ...
    'HardFailureCount',nnz(hardFail));

writetable(caseTable,fullfile(qaRoot,'qa_validation_cases.csv'));
writetable(checkTable,fullfile(qaRoot,'qa_validation_checks.csv'));
writetable(sensitivity,fullfile(qaRoot,'qa_validation_sensitivity.csv'));
fprintf('Validation QA: %d/%d end.mat, %d hard failure(s), %d checks.\n', ...
    nnz(isfile(caseSpec.EndFile)),height(caseSpec),nnz(hardFail),height(checkTable));
disp(caseTable);
disp(sensitivity);
if any(hardFail)
    disp(checkTable(hardFail,:));
    error('ValidationQA:HardFailure','%d validation checks failed.',nnz(hardFail));
end
end

function [summary,rows] = validate_case(spec)
summary = empty_metrics();
summary.Case = spec.Case;
summary.Group = spec.Group;
summary.EndFile = spec.EndFile;
rows = cell(0,7);
if ~isfile(spec.EndFile)
    rows(end+1,:) = check_row(spec.Case,'end_mat_exists','HARD',false, ...
        spec.EndFile,'existing end.mat','Only final states are accepted.');
    return
end
rows(end+1,:) = check_row(spec.Case,'end_mat_exists','HARD',true, ...
    spec.EndFile,'existing end.mat','Final state exists.');
S = load(spec.EndFile);
requiredTop = {'mesh','physics','solver','dt','tvec','results','n_max','tmax'};
missing = setdiff(requiredTop,fieldnames(S));
rows(end+1,:) = check_row(spec.Case,'top_level_variables','HARD',isempty(missing), ...
    strjoin(missing,','),strjoin(requiredTop,','),'Required saved variables.');
if ~isempty(missing), return; end

t = double(S.tvec(:));
r = S.results;
summary.TimeDays = t(end)/86400;
summary.DfMax = last(r,'Df_max');
summary.DfP99 = last(r,'Df_p99');
summary.PhaseP99 = last(r,'phi_p99');
summary.ScreeningIndex = last(r,'screening_index');
summary.DarcyQ = last(r,'leakage_Q');
summary.PoiseuilleQ = last(r,'leakage_Q_poiseuille');
summary.TopPcMPa = last(r,'contact_top_mean')/1e6;
summary.SpanningPath = logical(last(r,'degradation_connected'));
summary.DfAt3d = value_at_time(t,r,'Df_max',3*86400);
summary.DarcyQAt3d = value_at_time(t,r,'leakage_Q',3*86400);
summary.PoiseuilleQAt3d = value_at_time(t,r,'leakage_Q_poiseuille',3*86400);

controlOK = close_num(S.dt,spec.InitialDt_s,1e-12) && ...
    S.n_max==spec.NMax && close_num(S.tmax,3*86400,1e-12);
rows(end+1,:) = check_row(spec.Case,'saved_controls','HARD',controlOK, ...
    sprintf('dt=%.15g,nMax=%g,tmax=%.15g',S.dt,S.n_max,S.tmax), ...
    sprintf('dt=%g,nMax=%g,tmax=259200',spec.InitialDt_s,spec.NMax), ...
    'Validation controls match the declared case.');
timeOK = numel(t)>=2 && t(1)==0 && all(isfinite(t)) && all(diff(t)>0) && ...
    t(end)>S.tmax && t(end)-S.tmax<=6*3600+1e-8;
rows(end+1,:) = check_row(spec.Case,'time_vector','HARD',timeOK, ...
    sprintf('points=%d,finalDays=%.9g',numel(t),t(end)/86400), ...
    'finite increasing, final in (3,3.25] days','Adaptive schedule is complete.');
rows(end+1,:) = check_row(spec.Case,'physics_time_alignment','HARD', ...
    close_num(S.physics.time,t(end),1e-10),S.physics.time,t(end), ...
    'Physics and reporting times agree.');

g = S.mesh.getGroupIndex("Internal");
nElem = size(S.mesh.Elementgroups{g}.Elems,1);
nxNodes = numel(unique(S.mesh.Nodes(:,1)));
nyNodes = numel(unique(S.mesh.Nodes(:,2)));
meshOK = nElem==spec.Nx*spec.Ny && nxNodes==2*spec.Nx+1 && nyNodes==2*spec.Ny+1;
rows(end+1,:) = check_row(spec.Case,'mesh_manifest','HARD',meshOK, ...
    sprintf('%d,%dx%d',nElem,nxNodes,nyNodes), ...
    sprintf('%d,%dx%d',spec.Nx*spec.Ny,2*spec.Nx+1,2*spec.Ny+1), ...
    'Actual Q9 mesh matches the case label.');

requiredResults = ["CL_avg","CL_max","phi_p99","alpha_p99","Df_p99", ...
    "screening_index","Df_max","totalCycles","leakage_Q", ...
    "leakage_Q_poiseuille","leakage_path_cost","leakage_time_darcy", ...
    "leakage_time_poiseuille","degradation_connected"];
missingResults = setdiff(requiredResults,string(fieldnames(r)));
rows(end+1,:) = check_row(spec.Case,'result_fields','HARD',isempty(missingResults), ...
    strjoin(missingResults,','),strjoin(requiredResults,','),'Required histories.');
if ~isempty(missingResults), return; end

lengthOK = true;
finiteOK = true;
sentinels = ["leakage_path_cost","leakage_time_darcy","leakage_time_poiseuille"];
for name = requiredResults
    values = double(r.(name)(:));
    if ismember(name,["leakage_time_darcy","leakage_time_poiseuille"])
        lengthOK = lengthOK && isscalar(values);
    else
        lengthOK = lengthOK && numel(values)==numel(t);
    end
    if ~ismember(name,sentinels)
        finiteOK = finiteOK && all(isfinite(values));
    end
end
rows(end+1,:) = check_row(spec.Case,'history_lengths','HARD',lengthOK, ...
    numel(t),'histories match tvec; leak times are scalar sentinels', ...
    'No truncated series.');
pathCost = double(r.leakage_path_cost(:));
leakTimes = [double(r.leakage_time_darcy(:));double(r.leakage_time_poiseuille(:))];
sentinelOK = all((isfinite(pathCost)&pathCost>=0) | (isinf(pathCost)&pathCost>0)) && ...
    all(isfinite(leakTimes) | isnan(leakTimes));
rows(end+1,:) = check_row(spec.Case,'finite_and_sentinels','HARD',finiteOK&&sentinelOK, ...
    sprintf('finite=%d,sentinel=%d',finiteOK,sentinelOK), ...
    'finite except declared +Inf/NaN sentinels','Nonfinite values are intentional.');

fatigue = find_model(S.physics,'FatigueDamage');
phase = find_model(S.physics,'SealPhaseFieldDamage');
modelOK = ~isempty(fatigue) && ~isempty(phase);
rows(end+1,:) = check_row(spec.Case,'required_models','HARD',modelOK, ...
    sprintf('fatigue=%d,phase=%d',~isempty(fatigue),~isempty(phase)), ...
    'both present','State checks require both models.');
if ~modelOK, return; end

Df = double(fatigue.Df_ip(:));
Dft = double(fatigue.Df_ip_trial(:));
trialErr = Inf;
if isequal(size(Df),size(Dft)) && ~isempty(Df)
    trialErr = max(abs(Df-Dft));
end
rows(end+1,:) = check_row(spec.Case,'fatigue_bounds','HARD', ...
    all(isfinite(Df))&&all(Df>=-1e-12)&&all(Df<=1+1e-12), ...
    sprintf('[%.6g,%.6g]',min(Df),max(Df)),'[0,1]','Bounded committed state.');
rows(end+1,:) = check_row(spec.Case,'fatigue_trial_committed','HARD',trialErr<=1e-12, ...
    trialErr,'<=1e-12','Final trial state was committed.');

cycles = double(r.totalCycles(:));
cycleErr = max(abs(cycles-fatigue.f_cycle*t));
rows(end+1,:) = check_row(spec.Case,'cycles_equal_frequency_time','HARD', ...
    cycleErr<=1e-8*(1+cycles(end)),cycleErr,'relative <=1e-8', ...
    'Cycle accumulation follows the saved time vector.');
screen = max([double(r.phi_p99(:)),double(r.alpha_p99(:)),double(r.Df_p99(:))],[],2);
screenErr = max(abs(screen-double(r.screening_index(:))));
rows(end+1,:) = check_row(spec.Case,'screening_definition','HARD',screenErr<=1e-12, ...
    screenErr,'<=1e-12','Stored screening envelope is exact.');
summary.HardPass = ~any(strcmp(string(rows(:,4)),"FAIL"));
end

function rows = check_summary_files(caseTable,gridCsv,timeCsv,fixedCsv)
rows = cell(0,7);
files = [string(gridCsv),string(timeCsv),string(fixedCsv)];
for file = files
    rows(end+1,:) = check_row('SUMMARY','summary_exists','HARD',isfile(file), ...
        file,'existing CSV','Final summaries must exist.'); %#ok<AGROW>
end
if ~all(isfile(files)), return; end
G = readtable(gridCsv,'TextType','string');
D = readtable(timeCsv,'TextType','string');
F = readtable(fixedCsv,'TextType','string');
rows(end+1,:) = compare_rows(caseTable,G,["G30x42","Contact-aware","G50x70"], ...
    ["G30x42","G40x56","G50x70"],"grid_summary_values");
rows(end+1,:) = compare_rows(caseTable,D,["dt20s","Contact-aware","dt45s"], ...
    ["dt20s","dt30s_baseline","dt45s"],"time_summary_values");
rows(end+1,:) = compare_fixed(caseTable,F);
end

function row = compare_rows(caseTable,T,sourceNames,csvNames,checkId)
ok = height(T)==3 && isequal(string(T.Case),csvNames(:));
maxErr = Inf;
if ok
    maxErr = 0;
    for i=1:3
        src = caseTable(caseTable.Case==sourceNames(i),:);
        maxErr = max(maxErr,abs(src.DfAt3d-T.DfAt3d(i)));
        maxErr = max(maxErr,abs(src.DarcyQAt3d-T.DarcyQAt3d_m3s(i)));
        maxErr = max(maxErr,abs(src.PoiseuilleQAt3d-T.PoiseuilleQAt3d_m3s(i)));
        maxErr = max(maxErr,abs(src.TimeDays-T.FinalTime_days(i)));
    end
    ok = maxErr<=1e-12;
end
row = check_row('SUMMARY',checkId,'HARD',ok,maxErr,'<=1e-12', ...
    'CSV values are re-derived from the six final states.');
end

function row = compare_fixed(caseTable,T)
ok = height(T)==2 && isequal(string(T.Case),["Fixed BC";"Contact-aware"]);
maxErr = Inf;
if ok
    maxErr=0;
    for i=1:2
        src=caseTable(caseTable.Case==string(T.Case(i)),:);
        maxErr=max(maxErr,abs(src.DfMax-T.FinalDfMax(i)));
        maxErr=max(maxErr,abs(src.PoiseuilleQ-T.FinalPoiseuille_Q_m3s(i)));
        maxErr=max(maxErr,abs(src.TimeDays-T.FinalTime_days(i)));
    end
    ok=maxErr<=1e-12;
end
row=check_row('SUMMARY','fixed_contact_summary_values','HARD',ok,maxErr, ...
    '<=1e-12','Fixed/contact CSV is re-derived from final states.');
end

function T = build_sensitivity(C)
groups = ["grid";"time"];
labels = ["G30x42 / G40x56 / G50x70";"dt20s / dt30s / dt45s"];
T = table(groups,labels,zeros(2,1),zeros(2,1),zeros(2,1),zeros(2,1), ...
    'VariableNames',{'Group','Cases','DfMaxCoarseOrSmall','DfMaxBaseline', ...
    'DfMaxFineOrLarge','MaxAbsRelativeDf_pct'});
sets = {["G30x42","Contact-aware","G50x70"], ...
        ["dt20s","Contact-aware","dt45s"]};
for i=1:2
    values=zeros(3,1);
    for j=1:3, values(j)=C.DfAt3d(C.Case==sets{i}(j)); end
    T{i,3:5}=values(:).';
    T.MaxAbsRelativeDf_pct(i)=100*max(abs(values-values(2)))/max(abs(values(2)),eps);
end
end

function mdl = find_model(physics,className)
mdl=[];
for i=1:numel(physics.models)
    if isa(physics.models{i},className), mdl=physics.models{i}; return; end
end
end

function value = last(r,name)
value=double(r.(name)(end));
end

function value = value_at_time(t,r,name,target)
if t(1)>target || t(end)<target || ~isfield(r,name) || numel(r.(name))~=numel(t)
    error('ValidationQA:EvaluationTime','Cannot evaluate %s at %.15g s.',name,target);
end
value=interp1(double(t(:)),double(r.(name)(:)),target,'linear');
end

function tf = close_num(a,b,tol)
tf=isscalar(a)&&isscalar(b)&&isfinite(a)&&isfinite(b)&&abs(a-b)<=tol*(1+abs(b));
end

function row = check_row(caseName,id,severity,ok,value,limit,message)
if islogical(ok) && isscalar(ok) && ok, status="PASS"; else, status="FAIL"; end
row={string(caseName),string(id),string(severity),status,string(value),string(limit),string(message)};
end

function s = empty_metrics()
s=struct('Case',"",'Group',"",'EndFile',"",'TimeDays',NaN,'DfMax',NaN, ...
    'DfP99',NaN,'PhaseP99',NaN,'ScreeningIndex',NaN,'DarcyQ',NaN, ...
    'PoiseuilleQ',NaN,'TopPcMPa',NaN,'SpanningPath',false, ...
    'DfAt3d',NaN,'DarcyQAt3d',NaN,'PoiseuilleQAt3d',NaN,'HardPass',false);
end
