function run_module_verification()
%RUN_MODULE_VERIFICATION Independent analytical transport/fatigue checks.
%
% This function does not call or modify the production solver. It reproduces
% the implemented numerical formulas in isolated, auditable calculations:
%   1) 1-D Fick diffusion with Dirichlet-Neumann mixed boundaries, using
%      consistent-mass linear finite elements and backward Euler.
%   2) Constant-Wplus single-integration-point fatigue updates, including
%      below-cap and rate-capped per-cycle branches.

rootDir = fileparts(mfilename('fullpath'));
outDir = fullfile(rootDir, 'outputs');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

logPath = fullfile(outDir, 'module_verification.log');
if exist(logPath, 'file')
    delete(logPath);
end
diary(logPath);
diaryCleanup = onCleanup(@() diary('off'));

fprintf('Independent analytical module verification\n');
fprintf('Started: %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
fprintf('MATLAB: %s\n', version);
fprintf('Root: module_verification (archive-relative)\n\n');

oldFontName = get(groot, 'defaultAxesFontName');
oldFontSize = get(groot, 'defaultAxesFontSize');
oldTextFontName = get(groot, 'defaultTextFontName');
fontCleanup = onCleanup(@() restore_graphics_defaults( ...
    oldFontName, oldFontSize, oldTextFontName));
set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 10);
set(groot, 'defaultTextFontName', 'Times New Roman');

write_parameter_provenance(outDir);
print_acceptance_criteria();

fprintf('=== 1-D Fick diffusion check ===\n');
fick = verify_fick_mixed_boundary(outDir);

fprintf('\n=== Constant-Wplus fatigue point check ===\n');
fatigue = verify_fatigue_point(outDir);

status = [fick.Status; fatigue.Status];
writetable(status, fullfile(outDir, 'verification_status.csv'));

allPassed = all(status.Passed);
fprintf('\n=== Overall status ===\n');
disp(status);
fprintf('OVERALL_PASS=%d\n', allPassed);
fprintf('Completed: %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));

save(fullfile(outDir, 'module_verification_all_results.mat'), ...
    'fick', 'fatigue', 'status', 'allPassed');

diary off;
clear diaryCleanup fontCleanup;

if ~allPassed
    error('run_module_verification:VerificationFailed', ...
        'One or more independent verification checks failed.');
end
end


function result = verify_fick_mixed_boundary(outDir)
% Reproduce the production diffusion discretization in one spatial dimension.
% PDE: dC/dt - D*d2C/dx2 = 0 on 0 < x < L.
% BCs: C(0,t)=Cs and dC/dx(L,t)=0; IC: C(x,0)=0.

p = struct();
p.L = 3.5e-3;          % m, seal width
p.D = 6.0e-13;         % m^2/s, matrix diffusivity D_L
p.Cs = 0.65;           % prescribed dimensionless uptake at exposed boundary
p.tEnd = 20*86400;     % s
p.nSeries = 3000;

nxValues = [14, 70];
dtHoursValues = [12, 3, 0.25];

summaryRows = struct([]);
profileTables = cell(numel(nxValues)*numel(dtHoursValues), 1);
historyTables = cell(numel(nxValues)*numel(dtHoursValues), 1);
caseResults = struct([]);
row = 0;

for iNx = 1:numel(nxValues)
    for iDt = 1:numel(dtHoursValues)
        row = row + 1;
        nx = nxValues(iNx);
        dtHours = dtHoursValues(iDt);
        dt = dtHours*3600;
        dtLabel = strrep(sprintf('%g', dtHours), '.', 'p');
        caseId = sprintf('Nx%03d_dt%sh', nx, dtLabel);

        [x, t, Cnum, CavgNum, Cfinal, M] = solve_fick_fe_be( ...
            p.L, p.D, p.Cs, p.tEnd, nx, dt);
        Cexact = fick_profile_exact(x, p.tEnd, p.L, p.D, p.Cs, p.nSeries);
        CavgExact = fick_average_exact(t, p.L, p.D, p.Cs, p.nSeries);

        err = Cfinal - Cexact;
        l2Abs = sqrt(trapz(x, err.^2)/p.L);
        l2Den = sqrt(trapz(x, Cexact.^2)/p.L);
        l2Rel = l2Abs/max(l2Den, eps);
        linfAbs = max(abs(err));
        linfRel = linfAbs/max(max(abs(Cexact)), eps);
        avgFinalNum = (ones(size(x))'*(M*Cfinal))/p.L;
        avgFinalExact = CavgExact(end);
        avgAbs = abs(avgFinalNum-avgFinalExact);
        avgRel = avgAbs/max(abs(avgFinalExact), eps);

        summaryRows(row).Case = string(caseId);
        summaryRows(row).Nx = nx;
        summaryRows(row).NodeCount = nx+1;
        summaryRows(row).h_m = p.L/nx;
        summaryRows(row).dt_s = dt;
        summaryRows(row).dt_h = dtHours;
        summaryRows(row).StepCount = numel(t)-1;
        summaryRows(row).FinalTime_s = t(end);
        summaryRows(row).L2Abs = l2Abs;
        summaryRows(row).L2Rel = l2Rel;
        summaryRows(row).LinfAbs = linfAbs;
        summaryRows(row).LinfRel = linfRel;
        summaryRows(row).AverageNumeric = avgFinalNum;
        summaryRows(row).AverageExact = avgFinalExact;
        summaryRows(row).AverageAbsError = avgAbs;
        summaryRows(row).AverageRelError = avgRel;

        nProfile = numel(x);
        profileTables{row} = table( ...
            repmat(string(caseId), nProfile, 1), ...
            repmat(nx, nProfile, 1), repmat(dtHours, nProfile, 1), ...
            x, Cfinal, Cexact, err, abs(err), ...
            'VariableNames', {'Case','Nx','dt_h','x_m','C_numeric', ...
            'C_exact','Error','AbsError'});

        nHist = numel(t);
        historyTables{row} = table( ...
            repmat(string(caseId), nHist, 1), ...
            repmat(nx, nHist, 1), repmat(dtHours, nHist, 1), ...
            t, t/86400, CavgNum, CavgExact, CavgNum-CavgExact, ...
            'VariableNames', {'Case','Nx','dt_h','Time_s','Time_days', ...
            'AverageNumeric','AverageExact','AverageError'});

        caseResults(row).Case = caseId;
        caseResults(row).x = x;
        caseResults(row).t = t;
        caseResults(row).C_history = Cnum;
        caseResults(row).C_final = Cfinal;
        caseResults(row).C_exact_final = Cexact;
        caseResults(row).C_average_numeric = CavgNum;
        caseResults(row).C_average_exact = CavgExact;

        fprintf(['%-17s h=%8.3e m dt=%6.2f h steps=%4d ' ...
            'L2rel=%9.3e Linfrel=%9.3e AvgRel=%9.3e\n'], ...
            caseId, p.L/nx, dtHours, numel(t)-1, l2Rel, linfRel, avgRel);
    end
end

summary = struct2table(summaryRows);
profiles = vertcat(profileTables{:});
histories = vertcat(historyTables{:});
writetable(summary, fullfile(outDir, 'fick_summary.csv'));
writetable(profiles, fullfile(outDir, 'fick_profiles.csv'));
writetable(histories, fullfile(outDir, 'fick_average_history.csv'));

idxFine = summary.Nx == max(nxValues) & summary.dt_h == min(dtHoursValues);
idxGridCoarse = summary.Nx == min(nxValues) & summary.dt_h == min(dtHoursValues);
idxTimeCoarse = summary.Nx == max(nxValues) & summary.dt_h == max(dtHoursValues);

fineL2 = summary.L2Rel(idxFine);
fineLinf = summary.LinfRel(idxFine);
fineAvg = summary.AverageRelError(idxFine);
gridRatio = fineL2/summary.L2Rel(idxGridCoarse);
timeRatio = fineL2/summary.L2Rel(idxTimeCoarse);

testName = ["Fick fine-grid relative L2"; ...
    "Fick fine-grid relative Linf"; ...
    "Fick fine-grid mean relative error"; ...
    "Fick grid refinement error ratio"; ...
    "Fick time-step refinement error ratio"];
metric = ["L2Rel"; "LinfRel"; "AverageRelError"; ...
    "L2Rel_fine/L2Rel_coarse_grid"; "L2Rel_fine/L2Rel_coarse_dt"];
observed = [fineL2; fineLinf; fineAvg; gridRatio; timeRatio];
acceptance = ["<= 1.0e-2"; "<= 2.0e-2"; "<= 1.0e-2"; "< 1"; "< 1"];
passed = [fineL2 <= 1.0e-2; fineLinf <= 2.0e-2; fineAvg <= 1.0e-2; ...
    gridRatio < 1; timeRatio < 1];
status = table(testName, metric, observed, acceptance, passed, ...
    'VariableNames', {'Test','Metric','Observed','Acceptance','Passed'});

plot_fick_profiles(caseResults, p, outDir);
plot_fick_averages(caseResults, p, outDir);

result = struct('Parameters', p, 'Summary', summary, 'Profiles', profiles, ...
    'AverageHistories', histories, 'Cases', caseResults, 'Status', status);
save(fullfile(outDir, 'fick_verification_results.mat'), 'result');

fprintf('Fick verification status:\n');
disp(status);
end


function [x, t, Chistory, Cavg, C, M] = solve_fick_fe_be(L, D, Cs, tEnd, nx, dt)
% Consistent-mass linear FE and backward Euler, matching the production form.

nSteps = round(tEnd/dt);
assert(abs(nSteps*dt-tEnd) <= 100*eps(tEnd), ...
    'The selected time step must divide the final time exactly.');

x = linspace(0, L, nx+1)';
h = L/nx;
nNodes = nx+1;
M = spalloc(nNodes, nNodes, 3*nNodes);
K = spalloc(nNodes, nNodes, 3*nNodes);
Me = h/6*[2, 1; 1, 2];
Ke = D/h*[1, -1; -1, 1];
for e = 1:nx
    ids = [e, e+1];
    M(ids, ids) = M(ids, ids) + Me;
    K(ids, ids) = K(ids, ids) + Ke;
end

A = M/dt + K;
free = 2:nNodes;
Aff = A(free, free);
Afc = A(free, 1);

C = zeros(nNodes, 1);
C(1) = Cs;
t = (0:nSteps)'*dt;
Chistory = zeros(nNodes, nSteps+1);
Cavg = zeros(nSteps+1, 1);
Chistory(:,1) = C;
Cavg(1) = (ones(nNodes,1)'*(M*C))/L;

for k = 1:nSteps
    rhs = (M/dt)*C;
    Cnew = C;
    Cnew(1) = Cs;
    Cnew(free) = Aff \ (rhs(free)-Afc*Cs);
    C = Cnew;
    Chistory(:,k+1) = C;
    Cavg(k+1) = (ones(nNodes,1)'*(M*C))/L;
end
end


function C = fick_profile_exact(x, t, L, D, Cs, nTerms)
% Fourier-series solution for Dirichlet at x=0 and Neumann at x=L.
if t == 0
    C = zeros(size(x));
    C(1) = Cs;
    return;
end
n = 0:nTerms-1;
mu = (n+0.5)*pi;
decay = exp(-(D*t/L^2)*(mu.^2));
coeff = 2./mu;
C = Cs*(1-sin((x/L)*mu)*(coeff.*decay)');
C(1) = Cs;
end


function Cavg = fick_average_exact(t, L, D, Cs, nTerms)
% Exact spatial mean corresponding to fick_profile_exact.
n = 0:nTerms-1;
mu = (n+0.5)*pi;
weights = 2./(mu.^2);
decay = exp(-(t(:)*(D/L^2))*(mu.^2));
Cavg = Cs*(1-decay*weights');
Cavg(t(:) == 0) = 0;
end


function plot_fick_profiles(cases, p, outDir)
fig = figure('Visible','off','Color','w','Position',[100,100,1200,780]);
tiledlayout(fig, 2, 1, 'TileSpacing','compact', 'Padding','compact');

ax1 = nexttile;
xDense = linspace(0, p.L, 600)';
CexactDense = fick_profile_exact(xDense, p.tEnd, p.L, p.D, p.Cs, p.nSeries);
plot(ax1, xDense*1e3, CexactDense, 'k-', 'LineWidth', 2.2, ...
    'DisplayName', 'Analytical series');
hold(ax1, 'on');
styles = {'o-','s-','^-','d-','v-','>-'};
for i = 1:numel(cases)
    stride = max(1, floor(numel(cases(i).x)/18));
    plot(ax1, cases(i).x*1e3, cases(i).C_final, styles{i}, ...
        'LineWidth', 1.0, 'MarkerSize', 4, 'MarkerIndices', 1:stride:numel(cases(i).x), ...
        'DisplayName', strrep(cases(i).Case, '_', ', '));
end
xlabel(ax1, 'x (mm)');
ylabel(ax1, 'Uptake state, C_L (-)');
title(ax1, '1-D Fick diffusion at t = 20 days');
legend(ax1, 'Location','best');
grid(ax1, 'on'); box(ax1, 'on');

ax2 = nexttile;
hold(ax2, 'on');
for i = 1:numel(cases)
    err = cases(i).C_final-cases(i).C_exact_final;
    plot(ax2, cases(i).x*1e3, err, styles{i}, 'LineWidth', 1.0, ...
        'MarkerSize', 4, 'DisplayName', strrep(cases(i).Case, '_', ', '));
end
yline(ax2, 0, 'k:', 'HandleVisibility','off');
xlabel(ax2, 'x (mm)');
ylabel(ax2, 'C_{L,num} - C_{L,exact} (-)');
title(ax2, 'Pointwise uptake error');
legend(ax2, 'Location','best');
grid(ax2, 'on'); box(ax2, 'on');

exportgraphics(fig, fullfile(outDir, 'fick_profiles.png'), 'Resolution', 300);
close(fig);
end


function plot_fick_averages(cases, p, outDir)
fig = figure('Visible','off','Color','w','Position',[100,100,1200,780]);
tiledlayout(fig, 2, 1, 'TileSpacing','compact', 'Padding','compact');

ax1 = nexttile;
tDense = linspace(0, p.tEnd, 800)';
avgDense = fick_average_exact(tDense, p.L, p.D, p.Cs, p.nSeries);
plot(ax1, tDense/86400, avgDense, 'k-', 'LineWidth', 2.2, ...
    'DisplayName', 'Analytical series');
hold(ax1, 'on');
styles = {'o-','s-','^-','d-','v-','>-'};
for i = 1:numel(cases)
    stride = max(1, floor(numel(cases(i).t)/18));
    plot(ax1, cases(i).t/86400, cases(i).C_average_numeric, styles{i}, ...
        'LineWidth', 1.0, 'MarkerSize', 4, ...
        'MarkerIndices', 1:stride:numel(cases(i).t), ...
        'DisplayName', strrep(cases(i).Case, '_', ', '));
end
xlabel(ax1, 'Time (days)');
ylabel(ax1, 'Domain-mean uptake, C_L (-)');
title(ax1, 'Spatially averaged uptake');
legend(ax1, 'Location','best');
grid(ax1, 'on'); box(ax1, 'on');

ax2 = nexttile;
hold(ax2, 'on');
for i = 1:numel(cases)
    err = abs(cases(i).C_average_numeric-cases(i).C_average_exact);
    plot(ax2, cases(i).t(2:end)/86400, err(2:end), styles{i}, ...
        'LineWidth', 1.0, 'MarkerSize', 4, ...
        'DisplayName', strrep(cases(i).Case, '_', ', '));
end
set(ax2, 'YScale', 'log');
xlabel(ax2, 'Time (days)');
ylabel(ax2, 'Absolute domain-mean uptake error (-)');
title(ax2, 'Domain-mean uptake error');
legend(ax2, 'Location','best');
grid(ax2, 'on'); box(ax2, 'on');

exportgraphics(fig, fullfile(outDir, 'fick_average_history.png'), 'Resolution', 300);
close(fig);
end


function result = verify_fatigue_point(outDir)
% Reproduce the FatigueDamage constant-state update at one integration point.

p = struct();
p.W0 = 6.0e7;                 % J/m^3
p.beta = 2.0;
p.Pmin = 0.5e6;               % Pa
p.Pref = 21e6;                % Pa
p.pressureExponent = 1.2;
p.energyRateScale = 1.0;      % cycle^-1
p.contactFatigueCoeff = 0.25;
p.contactFatigueExponent = 1.0;
p.baseRate = 1.0e-6;          % cycle^-1
p.maxRate = 1.0e-5;           % cycle^-1
p.activationRatio = 0.05;
p.D0 = 0;
p.totalCycles = 120000;
p.cycleJumps = [1, 37, 1000];

caseDefs = struct([]);
caseDefs(1).Name = 'uncapped_base_below_threshold';
caseDefs(1).Wratio = 0.04;
caseDefs(1).Pmax = 21e6;
caseDefs(1).ExpectedBranch = 'uncapped';
caseDefs(1).ExpectedActivation = false;

caseDefs(2).Name = 'uncapped_active_low_amplitude';
caseDefs(2).Wratio = 0.10;
caseDefs(2).Pmax = p.Pmin + 1.0e-3*p.Pref;
caseDefs(2).ExpectedBranch = 'uncapped';
caseDefs(2).ExpectedActivation = true;

caseDefs(3).Name = 'capped_active_main_amplitude';
caseDefs(3).Wratio = 0.10;
caseDefs(3).Pmax = 21e6;
caseDefs(3).ExpectedBranch = 'capped';
caseDefs(3).ExpectedActivation = true;

summaryRows = struct([]);
historyTables = {};
caseRateRows = struct([]);
row = 0;
histRow = 0;

for iCase = 1:numel(caseDefs)
    c = caseDefs(iCase);
    Wplus = c.Wratio*p.W0;
    pressureAmplitude = max(0, c.Pmax-p.Pmin);
    pressureRatio = pressureAmplitude/p.Pref;
    activation = Wplus > p.activationRatio*p.W0;

    % Use pc/Pcref = 1 to exercise the same contact multiplier parameters.
    contactPressure = c.Pmax;
    contactPressureRef = c.Pmax;
    Mc = 1+p.contactFatigueCoeff* ...
        (contactPressure/contactPressureRef)^p.contactFatigueExponent;

    energeticRate = p.energyRateScale*double(activation)* ...
        pressureRatio^p.pressureExponent*(Wplus/p.W0)^p.beta;
    baseRateTerm = p.baseRate*pressureRatio^p.beta;
    rawRate = (energeticRate+baseRateTerm)*Mc;
    appliedRate = min(rawRate, p.maxRate);
    if rawRate > p.maxRate
        branch = 'capped';
    else
        branch = 'uncapped';
    end

    assert(strcmp(branch, c.ExpectedBranch), ...
        'Unexpected capped/uncapped branch for %s.', c.Name);
    assert(activation == c.ExpectedActivation, ...
        'Unexpected activation state for %s.', c.Name);

    caseRateRows(iCase).Case = string(c.Name);
    caseRateRows(iCase).Wplus_Jm3 = Wplus;
    caseRateRows(iCase).WplusOverW0 = c.Wratio;
    caseRateRows(iCase).Pmax_Pa = c.Pmax;
    caseRateRows(iCase).Pmin_Pa = p.Pmin;
    caseRateRows(iCase).PressureAmplitude_Pa = pressureAmplitude;
    caseRateRows(iCase).PressureRatio = pressureRatio;
    caseRateRows(iCase).Activation = activation;
    caseRateRows(iCase).ContactMultiplier = Mc;
    caseRateRows(iCase).EnergeticRate = energeticRate;
    caseRateRows(iCase).BaseRateTerm = baseRateTerm;
    caseRateRows(iCase).RawRatePerCycle = rawRate;
    caseRateRows(iCase).AppliedRatePerCycle = appliedRate;
    caseRateRows(iCase).Branch = string(branch);

    for iJump = 1:numel(p.cycleJumps)
        dN = p.cycleJumps(iJump);
        [N, Dnum, Dexact, maxHistoryError] = fatigue_discrete_history( ...
            p.D0, appliedRate, p.totalCycles, dN);

        row = row+1;
        summaryRows(row).Case = string(c.Name);
        summaryRows(row).Activation = activation;
        summaryRows(row).Branch = string(branch);
        summaryRows(row).Wplus_Jm3 = Wplus;
        summaryRows(row).WplusOverW0 = c.Wratio;
        summaryRows(row).PressureAmplitude_Pa = pressureAmplitude;
        summaryRows(row).PressureRatio = pressureRatio;
        summaryRows(row).ContactMultiplier = Mc;
        summaryRows(row).RawRatePerCycle = rawRate;
        summaryRows(row).AppliedRatePerCycle = appliedRate;
        summaryRows(row).CycleJump = dN;
        summaryRows(row).TotalCycles = p.totalCycles;
        summaryRows(row).NumericFinalD = Dnum(end);
        summaryRows(row).ExactFinalD = Dexact(end);
        summaryRows(row).FinalAbsError = abs(Dnum(end)-Dexact(end));
        summaryRows(row).MaxHistoryAbsError = maxHistoryError;

        if dN == 1000
            histRow = histRow+1;
            nHist = numel(N);
            historyTables{histRow,1} = table( ...
                repmat(string(c.Name), nHist, 1), ...
                repmat(string(branch), nHist, 1), ...
                repmat(activation, nHist, 1), N, Dnum, Dexact, Dnum-Dexact, ...
                'VariableNames', {'Case','Branch','Activation','Cycles', ...
                'D_numeric','D_exact','Error'});
        end
    end

    fprintf(['%-34s activation=%d raw=%10.3e applied=%10.3e ' ...
        'branch=%s\n'], c.Name, activation, rawRate, appliedRate, branch);
end

summary = struct2table(summaryRows);
histories = vertcat(historyTables{:});
rates = struct2table(caseRateRows);
writetable(summary, fullfile(outDir, 'fatigue_point_summary.csv'));
writetable(histories, fullfile(outDir, 'fatigue_point_histories.csv'));
writetable(rates, fullfile(outDir, 'fatigue_point_rates.csv'));

branchPass = any(rates.Branch == "capped") && any(rates.Branch == "uncapped") && ...
    any(rates.Activation & rates.Branch == "uncapped") && ...
    any(rates.Activation & rates.Branch == "capped");
maxError = max(summary.MaxHistoryAbsError);
jumpInvariant = max(groupsummary(summary, 'Case', 'range', 'NumericFinalD').range_NumericFinalD) ...
    <= 1.0e-12;

testName = ["Fatigue below-cap and rate-capped branch coverage"; ...
    "Fatigue discrete versus closed-form history"; ...
    "Fatigue cycle-jump invariance"];
metric = ["required branches present"; "maximum absolute history error"; ...
    "maximum final-D range across cycle jumps"];
observed = [double(branchPass); maxError; ...
    max(groupsummary(summary, 'Case', 'range', 'NumericFinalD').range_NumericFinalD)];
acceptance = ["= 1"; "<= 5.0e-12"; "<= 1.0e-12"];
passed = [branchPass; maxError <= 5.0e-12; jumpInvariant];
status = table(testName, metric, observed, acceptance, passed, ...
    'VariableNames', {'Test','Metric','Observed','Acceptance','Passed'});

plot_fatigue_results(histories, rates, outDir);

result = struct('Parameters', p, 'CaseDefinitions', caseDefs, ...
    'Rates', rates, 'Summary', summary, 'Histories', histories, 'Status', status);
save(fullfile(outDir, 'fatigue_verification_results.mat'), 'result');

fprintf('Fatigue verification status:\n');
disp(status);
end


function [N, Dnum, Dexact, maxError] = fatigue_discrete_history(D0, rate, totalCycles, dN)
% D_{n+1}=min(D_n+rate*DeltaN,1) versus D(N)=min(D0+rate*N,1).
nAlloc = ceil(totalCycles/dN)+1;
N = zeros(nAlloc,1);
Dnum = zeros(nAlloc,1);
Dexact = zeros(nAlloc,1);
Dnum(1) = D0;
Dexact(1) = D0;
k = 1;
while N(k) < totalCycles
    step = min(dN, totalCycles-N(k));
    k = k+1;
    N(k) = N(k-1)+step;
    Dnum(k) = min(Dnum(k-1)+rate*step, 1);
    Dexact(k) = min(D0+rate*N(k), 1);
end
N = N(1:k);
Dnum = Dnum(1:k);
Dexact = Dexact(1:k);
maxError = max(abs(Dnum-Dexact));
end


function plot_fatigue_results(histories, rates, outDir)
fig = figure('Visible','off','Color','w','Position',[100,100,1200,780]);
tiledlayout(fig, 2, 1, 'TileSpacing','compact', 'Padding','compact');

ax1 = nexttile;
hold(ax1, 'on');
names = unique(histories.Case, 'stable');
styles = {'o','s','^'};
colors = lines(numel(names));
shortLabels = ["base rate, below threshold", ...
    "activated, below rate cap", "activated, rate capped"];
for i = 1:numel(names)
    T = histories(histories.Case == names(i), :);
    plot(ax1, T.Cycles, T.D_exact, '-', 'Color', colors(i,:), ...
        'LineWidth', 1.8, 'DisplayName', shortLabels(i)+" exact");
    plot(ax1, T.Cycles, T.D_numeric, styles{i}, 'Color', colors(i,:), ...
        'MarkerSize', 3.5, 'LineStyle','none', ...
        'MarkerIndices', 1:max(1,floor(height(T)/20)):height(T), ...
        'DisplayName', shortLabels(i)+" discrete");
end
xlabel(ax1, 'Cycles N');
ylabel(ax1, 'Fatigue state, D_f (-)');
title(ax1, 'Constant-input fatigue update: exact versus discrete');
legend(ax1, 'Location','eastoutside', 'Interpreter','none');
grid(ax1, 'on'); box(ax1, 'on');

ax2 = nexttile;
barData = [rates.RawRatePerCycle, rates.AppliedRatePerCycle];
bar(ax2, barData, 'grouped');
set(ax2, 'YScale','log', 'XTick', 1:height(rates), ...
    'XTickLabel', shortLabels, 'XTickLabelRotation', 0, ...
    'TickLabelInterpreter', 'none');
ylabel(ax2, 'Rate (cycle^{-1})');
title(ax2, 'Raw and applied per-cycle rates');
legend(ax2, {'Raw rate','After per-cycle rate cap'}, 'Location','best');
grid(ax2, 'on'); box(ax2, 'on');

exportgraphics(fig, fullfile(outDir, 'fatigue_point_verification.png'), ...
    'Resolution', 300);
close(fig);
end


function write_parameter_provenance(outDir)
parameter = ["L"; "D_L"; "C_L_boundary"; "W0"; "beta"; ...
    "P_min"; "P_ref_fatigue"; "m_P"; "r0"; "r_base"; "r_max"; ...
    "A_c"; "m_c"];
value = [3.5e-3; 6.0e-13; 0.65; 6.0e7; 2.0; ...
    0.5e6; 21e6; 1.2; 1.0; 1.0e-6; 1.0e-5; 0.25; 1.0];
unit = ["m"; "m^2/s"; "-"; "J/m^3"; "-"; "Pa"; "Pa"; "-"; ...
    "cycle^-1"; "cycle^-1"; "cycle^-1"; "-"; "-"];
source = ["main_seal.m geometry"; "main_seal.m FluidDiffusion input"; ...
    "calibrated_degradation_options.m"; "calibrated_degradation_options.m"; ...
    "calibrated_degradation_options.m"; "main_seal.m FatigueDamage input"; ...
    "calibrated_degradation_options.m"; "calibrated_degradation_options.m"; ...
    "calibrated_degradation_options.m"; "calibrated_degradation_options.m"; ...
    "calibrated_degradation_options.m"; "main_seal.m FatigueDamage input"; ...
    "main_seal.m FatigueDamage input"];
T = table(parameter, value, unit, source, ...
    'VariableNames', {'Parameter','Value','Unit','Source'});
writetable(T, fullfile(outDir, 'parameter_provenance.csv'));
end


function print_acceptance_criteria()
fprintf('Predeclared acceptance criteria:\n');
fprintf('  Fick fine-grid relative L2             <= 1.0e-2\n');
fprintf('  Fick fine-grid relative Linf           <= 2.0e-2\n');
fprintf('  Fick fine-grid mean relative error     <= 1.0e-2\n');
fprintf('  Fick grid-refinement L2 error ratio    <  1 (common dt=0.25 h)\n');
fprintf('  Fick time-refinement L2 error ratio    <  1 (Nx=70)\n');
fprintf('  Fatigue required rate-branch coverage  =  1\n');
fprintf('  Fatigue max closed-form history error  <= 5.0e-12\n');
fprintf('  Fatigue final-D cycle-jump range       <= 1.0e-12\n\n');
end


function restore_graphics_defaults(fontName, fontSize, textFontName)
set(groot, 'defaultAxesFontName', fontName);
set(groot, 'defaultAxesFontSize', fontSize);
set(groot, 'defaultTextFontName', textFontName);
end
