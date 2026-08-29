function [screeningIndex, components] = compute_screening_index(data)
%COMPUTE_SCREENING_INDEX Percentile envelope used for within-study screening.
%
% The prescribed concentration boundary is excluded from the envelope.
% Fluid uptake already enters the degradation closures, so including its
% boundary-controlled maximum as another component would obscure comparison.

res = data.results;
targetSize = size(res.LFrac(:));

components.phase = bounded_field(res, 'phi_p99', targetSize, 0);
components.fatigue = bounded_field(res, 'Df_p99', targetSize, ...
    bounded_field(res, 'Df_max', targetSize, 0));
components.aging = bounded_field(res, 'alpha_p99', targetSize, ...
    bounded_field(res, 'alpha_max', targetSize, 0));
if isfield(res, 'CL_avg')
    components.fluid = min(1, max(0, res.CL_avg(:)));
else
    components.fluid = zeros(size(components.phase));
end

screeningIndex = max([components.phase, components.fatigue, ...
    components.aging], [], 2);
end

function value = bounded_field(res, name, targetSize, fallback)
if isfield(res, name) && ~isempty(res.(name))
    value = min(1, max(0, res.(name)(:)));
elseif isscalar(fallback)
    value = repmat(fallback, targetSize);
else
    value = fallback(:);
end
end
