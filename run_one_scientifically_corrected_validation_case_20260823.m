function endFile = run_one_scientifically_corrected_validation_case_20260823(caseKey)
%RUN_ONE_SCIENTIFICALLY_CORRECTED_VALIDATION_CASE_20260823 Run one isolated check.
%
% This worker uses the same parameters as
% run_scientifically_corrected_validation_20260823.  A small lock file lets
% the parent aggregation job wait for an already-running case instead of
% starting a duplicate process in the same output directory.

root = fileparts(mfilename('fullpath'));
validationRoot = fullfile(root, ...
    'Results_Paper_ScientificallyCorrected_20260823', 'Validation');
caseRoot = fullfile(validationRoot, 'cases');
lockRoot = fullfile(validationRoot, 'locks');
if ~isfolder(lockRoot), mkdir(lockRoot); end

caseKey = string(caseKey);
opts = common_options(caseRoot);
switch caseKey
    case "G50x70"
        opts.mesh = struct('Nx', 50, 'Ny', 70);
        opts.nameSuffix = "grid_50x70";
    case "dt20s"
        opts.initialDt = 20;
        opts.nMax = 215;
        opts.nameSuffix = "dt_20s";
    case "dt45s"
        opts.initialDt = 45;
        opts.nMax = 170;
        opts.nameSuffix = "dt_45s";
    otherwise
        error('ValidationWorker:UnknownCase', ...
            'Unknown validation case: %s', caseKey);
end

caseDir = fullfile(caseRoot, "Standard_21MPa_80C_" + opts.nameSuffix);
endFile = fullfile(caseDir, 'end.mat');
if isfile(endFile)
    fprintf('VALIDATION_CASE_ALREADY_COMPLETE %s\n', caseKey);
    return;
end

lockFile = fullfile(lockRoot, opts.nameSuffix + ".running");
fid = fopen(lockFile, 'w');
if fid < 0
    error('ValidationWorker:LockCreateFailed', ...
        'Could not create lock file: %s', lockFile);
end
fprintf(fid, 'PID=%d\nStarted=%s\n', feature('getpid'), ...
    char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
fclose(fid);
lockCleanup = onCleanup(@() delete_if_exists(lockFile)); %#ok<NASGU>

main_seal(1, opts);
if ~isfile(endFile)
    error('ValidationWorker:MissingEndFile', ...
        'Run returned without the expected end.mat: %s', endFile);
end
fprintf('VALIDATION_CASE_COMPLETE %s\n', caseKey);
end

function opts = common_options(caseRoot)
opts = struct();
opts.saveRoot = caseRoot;
opts.tmaxDays = 3;
opts.nMax = 180;
opts.plotEvery = inf;
opts.saveEvery = inf;
opts.verboseEvery = inf;
opts.useParpool = false;
opts.cleanResults = true;
opts.forceRestart = false;
opts.maxThreads = 4;
opts.alpha_coupling = 0.8;
opts.contactFatigueCoeff = 0.10;
opts.leakage_k_phi = 8e-18;
opts.leakagePhiExponent = 4.0;
opts.leakageAperturePhi = 1.2e-5;
opts.leakageApertureExponent = 3.0;
opts.leakageContactChi = 5.5;
opts.useContactPenalty = true;
end

function delete_if_exists(filePath)
if isfile(filePath)
    delete(filePath);
end
end
