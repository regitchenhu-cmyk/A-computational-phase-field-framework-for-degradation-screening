function run_scientifically_corrected_case_20260823(caseName)
%RUN_SCIENTIFICALLY_CORRECTED_CASE_20260823 Run one non-overwriting study case.

root = fileparts(mfilename('fullpath'));
resultRoot = fullfile(root, 'Results_Parametric_ScientificallyCorrected_20260823');
run_one_corrected_case(caseName, resultRoot, 4);
end
