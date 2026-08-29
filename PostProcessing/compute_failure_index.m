function [screeningIndex, components] = compute_failure_index(data)
%COMPUTE_FAILURE_INDEX Deprecated compatibility wrapper.
%
% New submission-facing code calls compute_screening_index.  This wrapper is
% retained only so archived local analysis scripts do not break.

[screeningIndex, current] = compute_screening_index(data);
components = current;
components.crack = current.phase;
end
