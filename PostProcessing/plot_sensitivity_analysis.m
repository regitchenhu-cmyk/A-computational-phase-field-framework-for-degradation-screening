function plot_sensitivity_analysis(folders, labels, opts)
    %PLOT_SENSITIVITY_ANALYSIS Scenario-based relative response analysis.
    %This is a post hoc comparison, not a full design-of-experiments study.

    if nargin < 3 || isempty(opts)
        opts = pub_style();
    end

    T = collect_scenario_metrics(folders, labels);
    if height(T) < 2
        warning('At least two scenarios are required for sensitivity analysis.');
        return
    end

    base = T(1,:);
    responses = [T.CrackLengthMm, T.FluidAvg, T.AgingMax, T.FatigueMax, T.FailureIndex];
    baseValues = [base.CrackLengthMm, base.FluidAvg, base.AgingMax, ...
                  base.FatigueMax, base.FailureIndex];
    rel = responses ./ max(baseValues, eps) - 1;

    fig = figure('Visible', 'on');
    set(fig, 'Position', opts.fig_double, 'Color', 'w');

    b = bar(rel(2:end,:) * 100, 'grouped');
    for k = 1:numel(b)
        b(k).FaceColor = opts.colors(k,:);
    end

    yline(0, '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 0.6);
    set(gca, 'XTickLabel', T.Label(2:end), opts.axes_props{:});
    ylabel('Change from baseline [\%]', 'Interpreter', 'latex', ...
        'FontSize', opts.fontsize_label);
    lg = legend({'$L_{\mathrm{frac}}$', '$\bar{C}_L$', ...
        '$\alpha_{\max}$', '$D_{f,\max}$', 'FI'}, ...
        'Interpreter', 'latex', 'FontSize', opts.fontsize_legend, ...
        'Location', 'best');
    lg.Box = 'off';

    export_fig(fig, 'fig14_relative_sensitivity', opts, 'both');
end
