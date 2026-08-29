function plot_reliability(data, opts, outdir)
    %PLOT_RELIABILITY Comprehensive reliability assessment figure.
    %
    % Panel layout:
    %   (a) Degradation factor stacked area chart over time
    %   (b) Failure index evolution
    %   (c) Remaining life fraction estimate
    %   (d) Radar/spider chart of degradation modes

    tvec    = data.tvec;
    results = data.results;
    t_days  = tvec * opts.s2day;
    
    fig = figure('Visible', 'on');
    set(fig, 'Position', opts.fig_wide);
    set(fig, 'Color', 'w');
    
    tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    %% (a) Degradation factors stacked area
    ax1 = nexttile;
    
    % Normalize each degradation mode to [0,1]
    if isfield(results, 'phi_p99')
        crack_norm = min(1, results.phi_p99);
    else
        Lx = max(data.mesh.Nodes(:,1));
        crack_norm = min(1, results.LFrac / max(Lx, eps));
    end
    fluid_norm = min(1, results.CL_avg);
    if isfield(results, 'alpha_p99'), aging_norm = min(1, results.alpha_p99); else, aging_norm = min(1, results.alpha_max); end
    if isfield(results, 'Df_p99'), fatigue_norm = min(1, results.Df_p99); else, fatigue_norm = min(1, results.Df_max); end
    
    Y_area = [crack_norm(:), fatigue_norm(:), aging_norm(:), fluid_norm(:)];
    
    h = area(t_days, Y_area);
    h(1).FaceColor = opts.colors(7,:);  % Dark red - crack
    h(2).FaceColor = opts.colors(4,:);  % Purple - fatigue
    h(3).FaceColor = opts.colors(2,:);  % Orange - aging
    h(4).FaceColor = opts.colors(3,:);  % Green - fluid
    
    h(1).FaceAlpha = 0.7;
    h(2).FaceAlpha = 0.7;
    h(3).FaceAlpha = 0.7;
    h(4).FaceAlpha = 0.7;
    
    for i = 1:4
        h(i).EdgeColor = 'none';
    end
    
    lg = legend('Diffuse integrity P99', 'Fatigue P99', 'Aging P99', 'Fluid uptake (diagnostic)', ...
        'FontSize', opts.fontsize_legend, 'Location', 'northwest');
    lg.Box = 'off';
    xlabel('Time [days]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    ylabel('Normalized degradation [-]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    set(ax1, opts.axes_props{:});
    xlim([0 max(t_days)]);
    ylim([0 max(1, max(sum(Y_area,2))*1.1)]);
    add_subfig_label(ax1, 'a', opts);
    
    %% (b) Screening index over time
    ax2 = nexttile;
    
    FI = max([crack_norm(:), fatigue_norm(:), aging_norm(:)], [], 2);
    
    plot(t_days, FI, '-', 'Color', opts.colors(8,:), 'LineWidth', opts.linewidth_thick);
    hold on;
    
    xlabel('Time [days]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    ylabel('Screening index $I_S$ [-]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    set(ax2, opts.axes_props{:});
    xlim([0 max(t_days)]);
    ylim([0 1.1]);
    add_subfig_label(ax2, 'b', opts);
    
    %% (c) Individual degradation mode curves
    ax3 = nexttile;
    
    plot(t_days, crack_norm, '-', 'Color', opts.colors(7,:), ...
        'LineWidth', opts.linewidth); hold on;
    plot(t_days, fatigue_norm, '-', 'Color', opts.colors(4,:), ...
        'LineWidth', opts.linewidth);
    plot(t_days, aging_norm, '-', 'Color', opts.colors(2,:), ...
        'LineWidth', opts.linewidth);
    plot(t_days, fluid_norm, '-', 'Color', opts.colors(3,:), ...
        'LineWidth', opts.linewidth);
    
    yline(1.0, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
    
    lg = legend('$L_\mathrm{frac}/L_\mathrm{crit}$', ...
                '$D_{f,\max}$', ...
                '$\alpha_{\max}$', ...
                '$\bar{C}_L/C_\mathrm{sat}$', ...
        'Interpreter', 'latex', 'FontSize', opts.fontsize_legend, ...
        'Location', 'northwest');
    lg.Box = 'off';
    xlabel('Time [days]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    ylabel('Normalized degradation [-]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    set(ax3, opts.axes_props{:});
    xlim([0 max(t_days)]);
    max_component = max([crack_norm(:); fatigue_norm(:); aging_norm(:); fluid_norm(:)]);
    ylim([0 max(0.1, max_component*1.15)]);
    add_subfig_label(ax3, 'c', opts);
    
    %% (d) Radar/spider chart of final degradation
    ax4 = nexttile;
    
    categories = {'Integrity P99', 'Fatigue P99', 'Aging P99', 'Fluid diagnostic', 'Integrity P99'};
    values = [crack_norm(end), fatigue_norm(end), aging_norm(end), ...
              fluid_norm(end), crack_norm(end)];
    
    n_cat = 4;
    theta = linspace(0, 2*pi, n_cat+1);
    
    % Draw radar background
    for r = [0.25, 0.5, 0.75, 1.0]
        xr = r * cos(theta);
        yr = r * sin(theta);
        plot(xr, yr, '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.4);
        hold on;
    end
    
    % Draw axes
    for i = 1:n_cat
        plot([0 cos(theta(i))], [0 sin(theta(i))], '-', ...
            'Color', [0.8 0.8 0.8], 'LineWidth', 0.4);
    end
    
    % Plot values
    xv = values .* cos(theta);
    yv = values .* sin(theta);
    fill(xv, yv, opts.colors(1,:), 'FaceAlpha', 0.25, 'EdgeColor', opts.colors(1,:), ...
        'LineWidth', opts.linewidth);
    plot(xv, yv, 'o', 'MarkerFaceColor', opts.colors(1,:), ...
        'MarkerEdgeColor', 'w', 'MarkerSize', opts.markersize);
    
    % Category labels
    label_r = 1.15;
    for i = 1:n_cat
        text(label_r*cos(theta(i)), label_r*sin(theta(i)), categories{i}, ...
            'FontName', opts.fontname, 'FontSize', opts.fontsize_tick, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
    
    % Value labels at rings
    for r = [0.25, 0.5, 0.75, 1.0]
        text(r*cos(theta(1))+0.03, r*sin(theta(1))+0.05, sprintf('%.0f%%', r*100), ...
            'FontName', opts.fontname, 'FontSize', 8, 'Color', [0.6 0.6 0.6]);
    end
    
    axis equal;
    axis off;
    xlim([-1.4 1.4]);
    ylim([-1.4 1.4]);
    
    % Title with FI value
    title(sprintf('Final state (I_S = %.3f)', FI(end)), ...
        'FontName', opts.fontname, 'FontSize', opts.fontsize_title);
    add_subfig_label(ax4, 'd', opts, 'northwest');
    
    export_fig(fig, 'fig05_reliability', opts, 'both');
end
