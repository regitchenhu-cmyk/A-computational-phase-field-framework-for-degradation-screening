function validation = finalize_scientifically_corrected_validation_20260823()
%FINALIZE_SCIENTIFICALLY_CORRECTED_VALIDATION_20260823 Aggregate completed checks.
%
% This function is read-only with respect to solver-state MAT files.  The
% patched run_paper_validation_study loads a completed case whenever its
% end.mat exists, and writes only the two compact validation summaries and
% its diagnostic preview.

root = fileparts(mfilename('fullpath'));
addpath(root);
addpath(genpath(fullfile(root, 'Models')));
addpath(genpath(fullfile(root, 'Shapes')));
studyRoot = fullfile(root, 'Results_Paper_ScientificallyCorrected_20260823');
validationRoot = fullfile(studyRoot, 'Validation');
baseFile = fullfile(studyRoot, 'FixedVsContact', ...
    'Standard_21MPa_80C_contactAware', 'end.mat');

required = [ ...
    string(baseFile); ...
    string(fullfile(validationRoot, 'cases', ...
        'Standard_21MPa_80C_grid_30x42', 'end.mat')); ...
    string(fullfile(validationRoot, 'cases', ...
        'Standard_21MPa_80C_grid_50x70', 'end.mat')); ...
    string(fullfile(validationRoot, 'cases', ...
        'Standard_21MPa_80C_dt_20s', 'end.mat')); ...
    string(fullfile(validationRoot, 'cases', ...
        'Standard_21MPa_80C_dt_45s', 'end.mat'))];
missing = required(~isfile(required));
if ~isempty(missing)
    error('ValidationFinalize:Incomplete', ...
        'Cannot finalize; missing:\n  %s', strjoin(missing, '\n  '));
end

opts = struct();
opts.saveRoot = validationRoot;
opts.baseFile = baseFile;
opts.tmaxDays = 3;
opts.nMax = 180;
opts.plotEvery = inf;
opts.saveEvery = inf;
opts.verboseEvery = inf;
opts.useParpool = false;
opts.cleanResults = false;
opts.forceRestart = false;
opts.maxThreads = 4;
opts.alpha_coupling = 0.8;
opts.contactFatigueCoeff = 0.10;
opts.leakage_k_phi = 8e-18;
opts.leakagePhiExponent = 4.0;
opts.leakageAperturePhi = 1.2e-5;
opts.leakageApertureExponent = 3.0;
opts.leakageContactChi = 5.5;

validation = run_paper_validation_study(opts);
fprintf('SCIENTIFICALLY_CORRECTED_VALIDATION_FINALIZED\n');
end
