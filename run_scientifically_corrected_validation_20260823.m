function validation = run_scientifically_corrected_validation_20260823()
%RUN_SCIENTIFICALLY_CORRECTED_VALIDATION_20260823 Recompute paper checks.
%
% Runs the fixed-boundary/contact-aware pair first, then uses the newly
% computed contact-aware case as the baseline for the grid and time-step
% sensitivity checks.  All outputs are isolated from archived results.

root = fileparts(mfilename('fullpath'));
studyRoot = fullfile(root, 'Results_Paper_ScientificallyCorrected_20260823');

common = struct();
common.tmaxDays = 3;
common.nMax = 180;
common.plotEvery = inf;
common.saveEvery = inf;
common.verboseEvery = inf;
common.useParpool = false;
common.cleanResults = true;
common.forceRestart = false;
common.maxThreads = 4;
common.alpha_coupling = 0.8;
common.contactFatigueCoeff = 0.10;
common.leakage_k_phi = 8e-18;
common.leakagePhiExponent = 4.0;
common.leakageAperturePhi = 1.2e-5;
common.leakageApertureExponent = 3.0;
common.leakageContactChi = 5.5;

comparisonOptions = common;
comparisonOptions.saveRoot = fullfile(studyRoot, 'FixedVsContact');
comparison = run_contact_leakage_comparison(1, comparisonOptions);
writetable(comparison.summary, fullfile(studyRoot, 'fixed_vs_contact_summary.csv'));
close all;

baseFile = fullfile(comparisonOptions.saveRoot, ...
    'Standard_21MPa_80C_contactAware', 'end.mat');
if ~isfile(baseFile)
    error('Validation:MissingBaseline', ...
        'New contact-aware baseline was not written: %s', baseFile);
end

validationOptions = common;
validationOptions.saveRoot = fullfile(studyRoot, 'Validation');
validationOptions.baseFile = baseFile;
validation = run_paper_validation_study(validationOptions);
close all;

fprintf('SCIENTIFICALLY_CORRECTED_VALIDATION_COMPLETE\n');
end
