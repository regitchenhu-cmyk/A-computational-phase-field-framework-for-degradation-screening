function compare_scenarios(folders, labels, opts)
    %COMPARE_SCENARIOS Compare results across multiple simulation scenarios.
    %
    % Usage:
    %   folders = {'./Results_Seal/Standard_21MPa_80C/', ...
    %              './Results_Seal/HighTemp_21MPa_120C/', ...
    %              './Results_Seal/HighPress_35MPa_80C/'};
    %   labels = {'Standard', 'High-Temp', 'High-Press'};
    %   opts = pub_style();
    %   compare_scenarios(folders, labels, opts);
    %
    % Generates:
    %   fig07: Crack length comparison
    %   fig08: Multi-mechanism degradation comparison
    %   fig09: Failure index comparison
    
    if nargin < 3
        opts = pub_style();
    end
    
    addpath(genpath('./PostProcessing'))
    
    outdir = './Figures';
    if ~exist(outdir, 'dir'); mkdir(outdir); end
    
    nscen = length(folders);
    
    % Load final results from each scenario
    all_data = cell(1, nscen);
    for s = 1:nscen
        % Find the last step or end.mat
        files = dir(fullfile(folders{s}, '*.mat'));
        steps = [];
        has_end = false;
        for i = 1:length(files)
            [~, fname, ~] = fileparts(files(i).name);
            num = str2double(fname);
            if ~isnan(num)
                steps = [steps, num];
            elseif strcmp(fname, 'end')
                has_end = true;
            end
        end
        
        if has_end
            all_data{s} = load(fullfile(folders{s}, 'end.mat'));
        elseif ~isempty(steps)
            last = max(steps);
            all_data{s} = load(fullfile(folders{s}, [num2str(last) '.mat']));
        else
            fprintf('WARNING: No data found in %s\n', folders{s});
            continue;
        end
        fprintf('Loaded scenario "%s": t_max = %.1f days\n', ...
            labels{s}, all_data{s}.tvec(end) * opts.s2day);
    end
    
    % Line styles for each scenario
    lstyles = {'-', '--', '-.', ':'};
    
    %% ===== Figure 7: Crack length comparison =====
    fig1 = figure('Visible', 'on');
    set(fig1, 'Position', opts.fig_single);
    set(fig1, 'Color', 'w');
    
    for s = 1:nscen
        if isempty(all_data{s}); continue; end
        t_d = all_data{s}.tvec * opts.s2day;
        Lf = all_data{s}.results.LFrac * opts.m2mm;
        plot(t_d, Lf, lstyles{mod(s-1,4)+1}, ...
            'Color', opts.colors(s,:), ...
            'LineWidth', opts.linewidth);
        hold on;
    end
    
    lg = legend(labels, 'Interpreter', 'latex', ...
        'FontSize', opts.fontsize_legend, 'Location', 'northwest');
    lg.Box = 'off';
    xlabel('Time [days]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    ylabel('Crack length $L_{\mathrm{frac}}$ [mm]', ...
        'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    set(gca, opts.axes_props{:});
    
    export_fig(fig1, 'fig07_crack_comparison', opts, 'both');
    
    
    %% ===== Figure 8: Multi-mechanism degradation comparison (2x2) =====
    fig2 = figure('Visible', 'on');
    set(fig2, 'Position', opts.fig_double);
    set(fig2, 'Color', 'w');
    
    tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % (a) Crack length
    ax1 = nexttile;
    for s = 1:nscen
        if isempty(all_data{s}); continue; end
        t_d = all_data{s}.tvec * opts.s2day;
        plot(t_d, all_data{s}.results.LFrac * opts.m2mm, lstyles{mod(s-1,4)+1}, ...
            'Color', opts.colors(s,:), 'LineWidth', opts.linewidth);
        hold on;
    end
    ylabel('$L_{\mathrm{frac}}$ [mm]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    xlabel('Time [days]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    set(ax1, opts.axes_props{:});
    add_subfig_label(ax1, 'a', opts);
    
    % (b) Fluid uptake
    ax2 = nexttile;
    for s = 1:nscen
        if isempty(all_data{s}); continue; end
        t_d = all_data{s}.tvec * opts.s2day;
        plot(t_d, all_data{s}.results.CL_avg, lstyles{mod(s-1,4)+1}, ...
            'Color', opts.colors(s,:), 'LineWidth', opts.linewidth);
        hold on;
    end
    ylabel('$\bar{C}_L$ [-]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    xlabel('Time [days]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    lg = legend(labels, 'Interpreter', 'latex', ...
        'FontSize', opts.fontsize_legend, 'Location', 'southeast');
    lg.Box = 'off';
    set(ax2, opts.axes_props{:});
    add_subfig_label(ax2, 'b', opts);
    
    % (c) Aging
    ax3 = nexttile;
    for s = 1:nscen
        if isempty(all_data{s}); continue; end
        t_d = all_data{s}.tvec * opts.s2day;
        plot(t_d, all_data{s}.results.alpha_max, lstyles{mod(s-1,4)+1}, ...
            'Color', opts.colors(s,:), 'LineWidth', opts.linewidth);
        hold on;
    end
    ylabel('$\alpha_{\max}$ [-]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    xlabel('Time [days]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    set(ax3, opts.axes_props{:});
    add_subfig_label(ax3, 'c', opts);
    
    % (d) Fatigue
    ax4 = nexttile;
    for s = 1:nscen
        if isempty(all_data{s}); continue; end
        t_d = all_data{s}.tvec * opts.s2day;
        plot(t_d, all_data{s}.results.Df_max, lstyles{mod(s-1,4)+1}, ...
            'Color', opts.colors(s,:), 'LineWidth', opts.linewidth);
        hold on;
    end
    ylabel('$D_{f,\max}$ [-]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    xlabel('Time [days]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    set(ax4, opts.axes_props{:});
    add_subfig_label(ax4, 'd', opts);
    
    export_fig(fig2, 'fig08_degradation_comparison', opts, 'both');
    
    
    %% ===== Figure 9: Screening-index comparison =====
    fig3 = figure('Visible', 'on');
    set(fig3, 'Position', opts.fig_single);
    set(fig3, 'Color', 'w');
    
    max_t = 0;
    for s = 1:nscen
        if ~isempty(all_data{s})
            max_t = max(max_t, all_data{s}.tvec(end) * opts.s2day);
        end
    end
    
    hold on;
    for s = 1:nscen
        if isempty(all_data{s}); continue; end
        t_d = all_data{s}.tvec * opts.s2day;
        FI = compute_failure_index(all_data{s});
        
        plot(t_d, FI, lstyles{mod(s-1,4)+1}, ...
            'Color', opts.colors(s,:), 'LineWidth', opts.linewidth_thick);
    end
    
    lg = legend(labels, 'Interpreter', 'latex', ...
        'FontSize', opts.fontsize_legend, 'Location', 'northwest');
    lg.Box = 'off';
    xlabel('Time [days]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    ylabel('Screening index $I_S$ [-]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    set(gca, opts.axes_props{:});
    xlim([0 max_t]);
    ylim([0 1.1]);
    
    export_fig(fig3, 'fig09_FI_comparison', opts, 'both');
    
    
    %% ===== Figure 10: Final-state bar chart comparison =====
    fig4 = figure('Visible', 'on');
    set(fig4, 'Position', opts.fig_single);
    set(fig4, 'Color', 'w');
    
    bar_data = zeros(nscen, 4);
    for s = 1:nscen
        if isempty(all_data{s}); continue; end
        [~, comp] = compute_failure_index(all_data{s});
        bar_data(s,:) = [comp.crack(end)*100, ...
                         comp.fatigue(end)*100, ...
                         comp.aging(end)*100, ...
                         comp.fluid(end)*100];
    end
    
    b = bar(bar_data);
    b(1).FaceColor = opts.colors(7,:);
    b(2).FaceColor = opts.colors(4,:);
    b(3).FaceColor = opts.colors(2,:);
    b(4).FaceColor = opts.colors(3,:);
    
    set(gca, 'XTickLabel', labels);
    ylabel('Degradation [\%]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    lg = legend('Integrity P99', 'Fatigue P99', 'Aging P99', 'Fluid diagnostic', ...
        'FontSize', opts.fontsize_legend, 'Location', 'northwest');
    lg.Box = 'off';
    set(gca, opts.axes_props{:});
    ylim([0 max(100, max(bar_data(:))*1.15)]);
    
    export_fig(fig4, 'fig10_scenario_bars', opts, 'both');
end
