function C = canonicalize_submission_summary(T)
%CANONICALIZE_SUBMISSION_SUMMARY Validate and order submission-facing fields.
%
% The corrected pipeline fails if a legacy intermediate table is supplied.
% This prevents archived column names from silently entering a new export.

C = T;
required = {'ConnectedDegradationExtentMm','ScreeningIndex','PhaseP99', ...
    'FatigueP99','AgingP99','TimeScreening07Days','TimeScreening1Days', ...
    'TimePhaseP99OneDays'};
missing = setdiff(required,C.Properties.VariableNames,'stable');
if ~isempty(missing)
    error('canonicalize_submission_summary:missingColumn', ...
        'Corrected summary lacks: %s',strjoin(missing,', '));
end

anchor = 'Cycles';
ordered = {'ConnectedDegradationExtentMm','ScreeningIndex','PhaseP99', ...
    'FatigueP99','AgingP99','TimeScreening07Days','TimeScreening1Days', ...
    'TimePhaseP99OneDays'};
ordered = intersect(ordered,C.Properties.VariableNames,'stable');
if ismember(anchor,C.Properties.VariableNames)
    C = movevars(C,ordered,'After',anchor);
end
end
