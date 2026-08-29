function plot_scenario_matrix(folders, labels, opts)
    %PLOT_SCENARIO_MATRIX Heatmap of final normalized degradation metrics.

    if nargin < 3 || isempty(opts)
        opts = pub_style();
    end

    T = collect_scenario_metrics(folders, labels);
    if isempty(T)
        warning('No scenario data available for matrix plot.');
        return
    end

    M = [T.CrackNorm, T.FatigueNorm, T.AgingNorm, T.FluidNorm, T.FailureIndex];
    metricLabels = {'Crack', 'Fatigue', 'Aging', 'Fluid', 'FI'};

    fig = figure('Visible', 'on');
    set(fig, 'Position', opts.fig_double, 'Color', 'w');

    imagesc(M);
    colormap(parula(256));
    caxis([0 max(1, max(M(:)))]);
    cb = colorbar;
    cb.Label.String = 'Normalized metric [-]';
    cb.Label.Interpreter = 'latex';
    cb.FontName = opts.fontname;
    cb.FontSize = opts.fontsize_tick;

    set(gca, 'XTick', 1:numel(metricLabels), 'XTickLabel', metricLabels, ...
        'YTick', 1:height(T), 'YTickLabel', T.Label, opts.axes_props{:});
    xtickangle(30);

    for i = 1:size(M,1)
        for j = 1:size(M,2)
            color = 'w';
            if M(i,j) < 0.55
                color = 'k';
            end
            text(j, i, sprintf('%.2f', M(i,j)), ...
                'HorizontalAlignment', 'center', 'Color', color, ...
                'FontName', opts.fontname, 'FontSize', opts.fontsize_tick);
        end
    end

    export_fig(fig, 'fig13_scenario_metric_matrix', opts, 'both');
end
