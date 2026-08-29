function plot_gc_degradation(data, opts, outdir)
    %PLOT_GC_DEGRADATION Visualize effective fracture energy Gc_eff.
    %
    % Shows how the three degradation mechanisms combine to reduce Gc
    % spatially across the seal cross-section.
    %
    % Panel layout:
    %   (a) f_fluid = 1 - eta_f * CL/CL_sat
    %   (b) f_aging = 1 - eta_a * alpha
    %   (c) f_fatigue = 1 - eta_fat * D_f
    %   (d) Gc_eff / Gc_0 = f_fluid * f_aging * f_fatigue
    
    physics = data.physics;
    mesh = data.mesh;
    scale = opts.m2mm;
    
    % Get degradation parameters from the phase-field model
    pf_model = physics.models{1};
    eta_fluid   = pf_model.fluidDegrade;
    eta_aging   = pf_model.agingDegrade;
    if isprop(pf_model, 'fatigueDegrade') && ~isempty(pf_model.fatigueDegrade)
        eta_fatigue = pf_model.fatigueDegrade;
    else
        eta_fatigue = 0.0;
    end
    
    % Get nodal values
    [~, stp_CL] = physics.dofSpace.getDofType({"CL"});
    [~, stp_ag] = physics.dofSpace.getDofType({"alpha_a"});
    [~, stp_df] = physics.dofSpace.getDofType({"D_f"});
    
    allNodes = mesh.GetAllNodesForGroup(1);
    
    CL_dofs = physics.dofSpace.getDofIndices(physics.dofSpace.getDofType({"CL"}), allNodes);
    AG_dofs = physics.dofSpace.getDofIndices(physics.dofSpace.getDofType({"alpha_a"}), allNodes);
    DF_dofs = physics.dofSpace.getDofIndices(physics.dofSpace.getDofType({"D_f"}), allNodes);
    
    CL_vals = max(0, physics.StateVec{stp_CL}(CL_dofs));
    AG_vals = max(0, min(1, physics.StateVec{stp_ag}(AG_dofs)));
    DF_vals = max(0, min(1, physics.StateVec{stp_df}(DF_dofs)));
    
    % Compute degradation factors
    f_fluid = max(0.02, 1 - eta_fluid * min(1, CL_vals / 1.0));
    f_aging = max(0.02, 1 - eta_aging * AG_vals);
    f_fatigue = max(0.02, 1 - eta_fatigue * DF_vals);
    f_total = f_fluid .* f_aging .* f_fatigue;
    f_total = max(f_total, 0.01);
    
    % Coordinates
    x_all = mesh.Nodes(allNodes, 1) * scale;
    y_all = mesh.Nodes(allNodes, 2) * scale;
    
    % Grid for interpolation
    nx = 400; ny = 400;
    [Xg, Yg] = meshgrid(linspace(min(x_all), max(x_all), nx), ...
                         linspace(min(y_all), max(y_all), ny));
    
    fig = figure('Visible', 'on');
    set(fig, 'Position', opts.fig_wide);
    set(fig, 'Color', 'w');
    
    tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Custom Gc colormap: green (good) -> yellow -> red (bad)
    cmap_gc = [linspace(0.2,0.9,128)' linspace(0.7,0.9,128)' linspace(0.2,0.2,128)';
               linspace(0.9,0.8,128)' linspace(0.9,0.15,128)' linspace(0.2,0.15,128)'];
    cmap_gc = flipud(cmap_gc);  % Red = low, Green = high
    
    %% (a) Fluid degradation factor
    ax1 = nexttile;
    Zg = griddata(x_all, y_all, f_fluid, Xg, Yg, 'linear');
    pcolor(Xg, Yg, Zg);
    shading interp;
    colormap(ax1, cmap_gc);
    caxis([0 1]);
    pub_colorbar('$f_{\rm fluid}$ [-]', opts);
    xlabel('$x$ [mm]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    ylabel('$y$ [mm]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    set(ax1, opts.axes_props{:});
    axis image;
    add_subfig_label(ax1, 'a', opts);
    
    %% (b) Aging degradation factor
    ax2 = nexttile;
    Zg = griddata(x_all, y_all, f_aging, Xg, Yg, 'linear');
    pcolor(Xg, Yg, Zg);
    shading interp;
    colormap(ax2, cmap_gc);
    caxis([0 1]);
    pub_colorbar('$f_{\rm aging}$ [-]', opts);
    xlabel('$x$ [mm]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    ylabel('$y$ [mm]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    set(ax2, opts.axes_props{:});
    axis image;
    add_subfig_label(ax2, 'b', opts);
    
    %% (c) Fatigue degradation factor
    ax3 = nexttile;
    Zg = griddata(x_all, y_all, f_fatigue, Xg, Yg, 'linear');
    pcolor(Xg, Yg, Zg);
    shading interp;
    colormap(ax3, cmap_gc);
    caxis([0 1]);
    pub_colorbar('$f_{\rm fatigue}$ [-]', opts);
    xlabel('$x$ [mm]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    ylabel('$y$ [mm]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    set(ax3, opts.axes_props{:});
    axis image;
    add_subfig_label(ax3, 'c', opts);
    
    %% (d) Combined: Gc_eff / Gc_0
    ax4 = nexttile;
    Zg = griddata(x_all, y_all, f_total, Xg, Yg, 'linear');
    pcolor(Xg, Yg, Zg);
    shading interp;
    colormap(ax4, cmap_gc);
    caxis([0 1]);
    pub_colorbar('$G_{c,\rm eff}/G_{c,0}$ [-]', opts);
    xlabel('$x$ [mm]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    ylabel('$y$ [mm]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    set(ax4, opts.axes_props{:});
    axis image;
    add_subfig_label(ax4, 'd', opts);
    
    % Add annotation for minimum Gc_eff
    [min_gc, min_idx] = min(f_total);
    title(sprintf('$\\min(G_{c,\\mathrm{eff}}/G_{c,0}) = %.3f$', min_gc), ...
        'Interpreter', 'latex', 'FontSize', opts.fontsize_title);
    
    export_fig(fig, 'fig06_gc_degradation', opts, 'both');
end
