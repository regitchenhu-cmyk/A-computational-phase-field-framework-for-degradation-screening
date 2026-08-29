function plot_time_histories(data, opts, outdir)
    %PLOT_TIME_HISTORIES Generate 2x2 panel of time-history curves.
    %
    % Panel layout:
    %   (a) Crack length vs time     (b) Average fluid uptake vs time
    %   (c) Aging (avg & max) vs time (d) Max fatigue damage vs time

    tvec    = data.tvec;
    results = data.results;
    
    % Time in days
    t_days = tvec * opts.s2day;
    
    fig = figure('Visible', 'on');
    set(fig, 'Position', opts.fig_double);
    set(fig, 'Color', 'w');
    
    tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    %% (a) Crack length
    ax1 = nexttile;
    plot(t_days, results.LFrac * opts.m2mm, '-', ...
        'Color', opts.colors(1,:), 'LineWidth', opts.linewidth);
    xlabel('Time [days]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    ylabel('Crack length $L_{\mathrm{frac}}$ [mm]', ...
        'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    set(ax1, opts.axes_props{:});
    xlim([0 max(t_days)]);
    ylim([0 max(0.1, max(results.LFrac*opts.m2mm)*1.1)]);
    add_subfig_label(ax1, 'a', opts);
    
    %% (b) Average fluid concentration
    ax2 = nexttile;
    plot(t_days, results.CL_avg, '-', ...
        'Color', opts.colors(3,:), 'LineWidth', opts.linewidth);
    hold on
    if isfield(results, 'CL_max')
        plot(t_days, results.CL_max, '--', ...
            'Color', opts.colors(3,:), 'LineWidth', opts.linewidth_thin);
        lg = legend('Average $\bar{C}_L$', 'Maximum $C_{L,\max}$', ...
            'Interpreter', 'latex', 'FontSize', opts.fontsize_legend, ...
            'Location', 'southeast');
        lg.Box = 'off';
    end
    xlabel('Time [days]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    ylabel('Fluid concentration $C_L$ [-]', ...
        'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    set(ax2, opts.axes_props{:});
    xlim([0 max(t_days)]);
    add_subfig_label(ax2, 'b', opts);
    
    %% (c) Aging
    ax3 = nexttile;
    plot(t_days, results.alpha_avg, '-', ...
        'Color', opts.colors(2,:), 'LineWidth', opts.linewidth);
    hold on
    plot(t_days, results.alpha_max, '--', ...
        'Color', opts.colors(2,:), 'LineWidth', opts.linewidth_thin);
    xlabel('Time [days]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    ylabel('Aging variable $\alpha_{\mathrm{age}}$ [-]', ...
        'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    lg = legend('Average $\bar{\alpha}$', 'Maximum $\alpha_{\max}$', ...
        'Interpreter', 'latex', 'FontSize', opts.fontsize_legend, ...
        'Location', 'southeast');
    lg.Box = 'off';
    set(ax3, opts.axes_props{:});
    xlim([0 max(t_days)]);
    ylim([0 max(0.01, max(results.alpha_max)*1.1)]);
    add_subfig_label(ax3, 'c', opts);
    
    %% (d) Fatigue damage
    ax4 = nexttile;
    plot(t_days, results.Df_max, '-', ...
        'Color', opts.colors(4,:), 'LineWidth', opts.linewidth);
    xlabel('Time [days]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    ylabel('Max fatigue damage $D_{f,\max}$ [-]', ...
        'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    set(ax4, opts.axes_props{:});
    xlim([0 max(t_days)]);
    ylim([0 max(0.01, max(results.Df_max)*1.1)]);
    add_subfig_label(ax4, 'd', opts);
    
    export_fig(fig, 'fig02_time_histories', opts, 'both');
    
    
    %% ===== Additional: Crack length vs cycle count =====
    if isfield(results, 'totalCycles')
        fig2 = figure('Visible', 'on');
        set(fig2, 'Position', opts.fig_single);
        set(fig2, 'Color', 'w');
        
        plot(results.totalCycles, results.LFrac * opts.m2mm, '-', ...
            'Color', opts.colors(1,:), 'LineWidth', opts.linewidth);
        xlabel('Number of pressure cycles $N$ [-]', ...
            'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
        ylabel('Crack length $L_{\mathrm{frac}}$ [mm]', ...
            'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
        set(gca, opts.axes_props{:});
        xlim([0 max(results.totalCycles)]);
        ylim([0 max(0.1, max(results.LFrac*opts.m2mm)*1.1)]);
        
        export_fig(fig2, 'fig02b_crack_vs_cycles', opts, 'both');
    end
end
