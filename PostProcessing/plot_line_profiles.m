function plot_line_profiles(data, opts, outdir)
    %PLOT_LINE_PROFILES Extract and plot 1D profiles across the seal.
    %
    % Generates:
    %   fig11: Horizontal profiles at mid-height (phi, CL, alpha, Df vs x)
    %   fig12: Vertical profiles at key x-locations (phi, sigma_H vs y)
    
    physics = data.physics;
    mesh = data.mesh;
    
    scale = opts.m2mm;
    
    Lx = max(mesh.Nodes(:,1));
    Ly = max(mesh.Nodes(:,2));
    
    allNodes = mesh.GetAllNodesForGroup(1);
    x_all = mesh.Nodes(allNodes, 1);
    y_all = mesh.Nodes(allNodes, 2);
    
    %% Extract horizontal profiles at mid-height
    y_mid = Ly / 2;
    tol_y = Ly / 50;  % tolerance for selecting nodes near mid-height
    
    mid_nodes = allNodes(abs(y_all - y_mid) < tol_y);
    [x_sorted, sort_idx] = sort(mesh.Nodes(mid_nodes, 1));
    mid_nodes_sorted = mid_nodes(sort_idx);
    
    % Get field values along the line
    fields = {"phi", "CL", "alpha_a", "D_f"};
    field_labels = {'$\phi$ [-]', '$C_L$ [-]', '$\alpha_{\mathrm{age}}$ [-]', '$D_f$ [-]'};
    field_colors = [opts.colors(7,:); opts.colors(3,:); opts.colors(2,:); opts.colors(4,:)];
    
    profiles = cell(1, length(fields));
    for f = 1:length(fields)
        [dofT, dofS] = physics.dofSpace.getDofType(fields(f));
        dofs = physics.dofSpace.getDofIndices(dofT, mid_nodes_sorted);
        profiles{f} = physics.StateVec{dofS}(dofs);
    end
    
    %% Figure 11: Horizontal profiles (2x2)
    fig1 = figure('Visible', 'on');
    set(fig1, 'Position', opts.fig_double);
    set(fig1, 'Color', 'w');
    
    tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    panel_labels = {'a', 'b', 'c', 'd'};
    
    for f = 1:length(fields)
        ax = nexttile;
        plot(x_sorted * scale, profiles{f}, '-', ...
            'Color', field_colors(f,:), 'LineWidth', opts.linewidth);
        xlabel('$x$ [mm]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
        ylabel(field_labels{f}, 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
        set(ax, opts.axes_props{:});
        xlim([0 Lx*scale]);
        add_subfig_label(ax, panel_labels{f}, opts);
    end
    
    export_fig(fig1, 'fig11_horizontal_profiles', opts, 'both');
    
    
    %% Extract vertical profiles at key x-locations
    x_locs = [0.25*Lx, 0.5*Lx, 0.75*Lx];  % 25%, 50%, 75% of seal width
    x_labels_str = {'$x/W = 0.25$', '$x/W = 0.50$', '$x/W = 0.75$'};
    tol_x = Lx / 50;
    
    fig2 = figure('Visible', 'on');
    set(fig2, 'Position', opts.fig_double);
    set(fig2, 'Color', 'w');
    
    tl = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    for loc = 1:length(x_locs)
        ax = nexttile;
        
        loc_nodes = allNodes(abs(x_all - x_locs(loc)) < tol_x);
        [y_sorted, sort_idx] = sort(mesh.Nodes(loc_nodes, 2));
        loc_nodes_sorted = loc_nodes(sort_idx);
        
        for f = 1:length(fields)
            [dofT, dofS] = physics.dofSpace.getDofType(fields(f));
            dofs = physics.dofSpace.getDofIndices(dofT, loc_nodes_sorted);
            vals = physics.StateVec{dofS}(dofs);
            
            plot(y_sorted * scale, vals, '-', ...
                'Color', field_colors(f,:), 'LineWidth', opts.linewidth);
            hold on;
        end
        
        xlabel('$y$ [mm]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
        ylabel('Field value [-]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
        title(x_labels_str{loc}, 'Interpreter', 'latex', 'FontSize', opts.fontsize_title);
        set(ax, opts.axes_props{:});
        xlim([0 Ly*scale]);
        
        if loc == length(x_locs)
            lg = legend('$\phi$', '$C_L$', '$\alpha$', '$D_f$', ...
                'Interpreter', 'latex', 'FontSize', opts.fontsize_legend, ...
                'Location', 'northeast');
            lg.Box = 'off';
        end
        
        add_subfig_label(ax, panel_labels{loc}, opts);
    end
    
    export_fig(fig2, 'fig12_vertical_profiles', opts, 'both');
end
