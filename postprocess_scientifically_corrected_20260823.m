function [T, report] = postprocess_scientifically_corrected_20260823()
%POSTPROCESS_SCIENTIFICALLY_CORRECTED_20260823 Final gated artifact pipeline.
%
% Run only after the 16-case campaign completes.  The preflight QA accepts
% no checkpoint substitute for end.mat.  Post-processing writes one canonical
% CSV, and the final QA re-derives all common values from MAT before accepting
% the table/figure data source.

root = fileparts(mfilename('fullpath'));
resultRoot = fullfile(root, 'Results_Parametric_ScientificallyCorrected_20260823');
outdir = fullfile(root, 'Figures_Parametric_ScientificallyCorrected_20260823');
summaryCsv = fullfile(outdir, 'parametric_summary.csv');

pre = struct('FailOnError',true,'RequireCsv',false,'SummaryCsv',summaryCsv);
qa_scientifically_corrected_campaign_20260823(resultRoot, pre);

T = postprocess_corrected_submission(resultRoot, outdir);

post = struct('FailOnError',true,'RequireCsv',true,'SummaryCsv',summaryCsv, ...
    'ReportDir',fullfile(outdir,'qa'));
report = qa_scientifically_corrected_campaign_20260823(resultRoot, post);
end
