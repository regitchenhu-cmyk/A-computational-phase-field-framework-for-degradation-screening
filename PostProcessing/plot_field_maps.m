function plot_field_maps(data, opts, outdir)
    %PLOT_FIELD_MAPS Generate 2x3 panel of field contour maps.
    %
    % Panel layout:
    %   (a) Phase-field phi      (b) Hydrostatic stress    (c) Fluid concentration
    %   (d) Aging variable       (e) Fatigue damage         (f) Von Mises stress
    
    physics = data.physics;
    mesh = data.mesh;
    
    fig = figure('Visible', 'on');
    set(fig, 'Position', opts.fig_wide);
    set(fig, 'Color', 'w');
    
    t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Coordinate scaling: m -> mm
    scale = opts.m2mm;
    
    %% (a) Phase-field damage
    ax1 = nexttile;
    plot_nodal_field(physics, mesh, "phi", 0, "Internal", scale);
    cmap_phi = [ones(128,1) linspace(1,0,128)' linspace(1,0,128)';  % white->red
                linspace(1,0,128)' zeros(128,1) zeros(128,1)];       % red->black
    colormap(ax1, flipud(cmap_phi));
    caxis([0 1]);
    cb = pub_colorbar('$\phi$ [-]', opts);
    xlabel('$x$ [mm]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    ylabel('$y$ [mm]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    set(ax1, opts.axes_props{:});
    axis image
    add_subfig_label(ax1, 'a', opts);
    
    %% (b) Hydrostatic stress
    ax2 = nexttile;
    plot_ip_field(physics, mesh, "sh", "Internal", scale, opts.Pa2MPa);
    colormap(ax2, opts.cmap_diverge);
    cb = pub_colorbar('$\sigma_H$ [MPa]', opts);
    xlabel('$x$ [mm]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    ylabel('$y$ [mm]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    set(ax2, opts.axes_props{:});
    axis image
    add_subfig_label(ax2, 'b', opts);
    
    %% (c) Fluid concentration
    ax3 = nexttile;
    plot_nodal_field(physics, mesh, "CL", 0, "Internal", scale);
    colormap(ax3, opts.cmap_fluid);
    CL_vals = get_nodal_values(physics, "CL");
    CL_max = max(0.01, max(CL_vals));
    CL_min = min(0, min(CL_vals));
    caxis([CL_min CL_max]);
    cb = pub_colorbar('$C_L$ [-]', opts);
    xlabel('$x$ [mm]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    ylabel('$y$ [mm]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    set(ax3, opts.axes_props{:});
    axis image
    add_subfig_label(ax3, 'c', opts);
    
    %% (d) Aging variable
    ax4 = nexttile;
    plot_nodal_field(physics, mesh, "alpha_a", 0, "Internal", scale);
    colormap(ax4, flipud(autumn(256)));  % yellow->red, visible on white
    AG_vals = get_nodal_values(physics, "alpha_a");
    caxis([0 max(0.01, max(AG_vals))]);
    cb = pub_colorbar('$\alpha_{age}$ [-]', opts);
    xlabel('$x$ [mm]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    ylabel('$y$ [mm]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    set(ax4, opts.axes_props{:});
    axis image
    add_subfig_label(ax4, 'd', opts);
    
    %% (e) Fatigue damage
    ax5 = nexttile;
    plot_nodal_field(physics, mesh, "D_f", 0, "Internal", scale);
    colormap(ax5, opts.cmap_seq_blue);  % white->blue
    DF_vals = get_nodal_values(physics, "D_f");
    caxis([0 max(0.01, max(DF_vals))]);
    cb = pub_colorbar('$D_f$ [-]', opts);
    xlabel('$x$ [mm]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    ylabel('$y$ [mm]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    set(ax5, opts.axes_props{:});
    axis image
    add_subfig_label(ax5, 'e', opts);
    
    %% (f) Deformed shape with phi overlay (auto-scaled)
    ax6 = nexttile;
    % Auto-determine displacement scale for visible deformation
    [~, dxSteps_u] = physics.dofSpace.getDofType({"dx","dy"});
    max_disp = max([max(abs(physics.StateVec{dxSteps_u(1)})), ...
                    max(abs(physics.StateVec{dxSteps_u(2)}))]);
    Lchar = max(mesh.Nodes(:,1));
    if max_disp > 0
        auto_scale = 0.1 * Lchar / max_disp;  % 10% of domain size
    else
        auto_scale = 1;
    end
    plot_nodal_field(physics, mesh, "phi", auto_scale, "Internal", scale);
    colormap(ax6, flipud(cmap_phi));
    caxis([0 1]);
    cb = pub_colorbar('$\phi$ [-]', opts);
    xlabel('$x$ [mm]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    ylabel('$y$ [mm]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    title(sprintf('Deformed ($\\times %.0f$)', auto_scale), ...
        'Interpreter', 'latex', 'FontSize', opts.fontsize_title);
    set(ax6, opts.axes_props{:});
    axis image
    add_subfig_label(ax6, 'f', opts);
    
    export_fig(fig, 'fig01_field_maps', opts, 'both');
end


%% ===== Helper: plot nodal field with coordinate scaling =====
function plot_nodal_field(physics, mesh, dofName, dispscale, plotloc, coord_scale)
    [dxTypes, dxSteps] = physics.dofSpace.getDofType({"dx";"dy";dofName});
    
    for g = 1:length(mesh.Elementgroups)
        if mesh.Elementgroups{g}.name == plotloc && mesh.Elementgroups{g}.type == "Q9"
            X = []; Y = []; Z = [];
            for el = 1:size(mesh.Elementgroups{g}.Elems, 1)
                elnodes = mesh.Elementgroups{g}.Elems(el,:);
                order = [1 3 9 7];
                
                zdofs = physics.dofSpace.getDofIndices(dxTypes(3), elnodes);
                if dispscale >= 0
                    xdofs = physics.dofSpace.getDofIndices(dxTypes(1), elnodes);
                    ydofs = physics.dofSpace.getDofIndices(dxTypes(2), elnodes);
                    X(el,:) = (mesh.Nodes(elnodes(order),1) + ...
                        dispscale*physics.StateVec{dxSteps(1)}(xdofs(order))) * coord_scale;
                    Y(el,:) = (mesh.Nodes(elnodes(order),2) + ...
                        dispscale*physics.StateVec{dxSteps(2)}(ydofs(order))) * coord_scale;
                else
                    X(el,:) = mesh.Nodes(elnodes(order),1) * coord_scale;
                    Y(el,:) = mesh.Nodes(elnodes(order),2) * coord_scale;
                end
                Z(el,:) = physics.StateVec{dxSteps(3)}(zdofs(order));
            end
            patch(X', Y', Z', Z', 'EdgeColor', 'None', 'FaceColor', 'interp');
            hold on;
        end
    end
end


%% ===== Helper: plot IP field with coordinate and value scaling =====
function plot_ip_field(physics, mesh, varName, plotloc, coord_scale, val_scale)
    for g = 1:length(mesh.Elementgroups)
        if mesh.Elementgroups{g}.name == plotloc
            x = []; y = []; z = [];
            for e = 1:size(mesh.Elementgroups{g}.Elems, 1)
                IPvar = physics.Request_Info(varName, e, "Interior");
                xy = mesh.getIPCoords(g, e);
                
                meanx = mean(xy(1,:));
                meany = mean(xy(2,:));
                for i = 1:size(xy,2)
                    if xy(1,i) < meanx
                        xy(1,i) = xy(1,i) + 1e-9;
                    else
                        xy(1,i) = xy(1,i) - 1e-9;
                    end
                    if xy(2,i) < meany
                        xy(2,i) = xy(2,i) + 1e-9;
                    else
                        xy(2,i) = xy(2,i) - 1e-9;
                    end
                end
                
                x(e,:) = xy(1,:) * coord_scale;
                y(e,:) = xy(2,:) * coord_scale;
                z(e,:) = IPvar * val_scale;
            end
            
            xlims = [min(x(:)), max(x(:))];
            ylims = [min(y(:)), max(y(:))];
            [xi,yi] = meshgrid(linspace(xlims(1),xlims(2),500), ...
                               linspace(ylims(1),ylims(2),500));
            zi = griddata(x(:), y(:), z(:), xi, yi, 'linear');
            s = surf(xi, yi, zi);
            s.EdgeColor = 'none';
            view(2);
            hold on;
        end
    end
end


%% ===== Helper: get nodal values for a given dof name =====
function vals = get_nodal_values(physics, dofName)
    [dxTypes, dxSteps] = physics.dofSpace.getDofType({dofName});
    vals = physics.StateVec{dxSteps(1)};
end
