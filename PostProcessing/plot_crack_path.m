function plot_crack_path(data, opts, outdir)
    %PLOT_CRACK_PATH Detailed crack path visualization with deformed mesh.
    %
    % Generates:
    %   (a) Undeformed reference with crack overlay
    %   (b) Deformed configuration with phi contour
    %   (c) Crack-tip detail (zoomed view)
    
    physics = data.physics;
    mesh = data.mesh;
    scale = opts.m2mm;
    
    fig = figure('Visible', 'on');
    set(fig, 'Position', [50 50 1200 500]);
    set(fig, 'Color', 'w');
    
    tl = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    %% (a) Undeformed with phi > 0.9 contour (crack path)
    ax1 = nexttile;
    
    % Plot full phi field
    plot_nodal_field_local(physics, mesh, "phi", 0, "Internal", scale);
    colormap(ax1, flipud(hot(256)));
    caxis([0 1]);
    hold on;
    
    % Overlay crack path contour (phi = 0.95)
    [X_grid, Y_grid, Phi_grid] = interpolate_nodal_field(physics, mesh, "phi", scale);
    if ~isempty(Phi_grid)
        [~, hc] = contour(X_grid, Y_grid, Phi_grid, [0.95 0.95], ...
            'LineColor', 'w', 'LineWidth', 1.5, 'LineStyle', '-');
    end
    
    % Draw boundary outline
    Lx = max(mesh.Nodes(:,1)) * scale;
    Ly = max(mesh.Nodes(:,2)) * scale;
    rectangle('Position', [0 0 Lx Ly], 'EdgeColor', 'k', 'LineWidth', 0.8);
    
    axis image;
    xlabel('$x$ [mm]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    ylabel('$y$ [mm]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    pub_colorbar('$\phi$ [-]', opts);
    set(ax1, opts.axes_props{:});
    add_subfig_label(ax1, 'a', opts);
    
    %% (b) Deformed configuration with magnified displacement
    ax2 = nexttile;
    disp_scale = 30;  % Magnification factor for visualization
    
    plot_nodal_field_local(physics, mesh, "phi", disp_scale, "Internal", scale);
    colormap(ax2, flipud(hot(256)));
    caxis([0 1]);
    
    axis image;
    xlabel('$x$ [mm]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    ylabel('$y$ [mm]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    title(sprintf('Deformed ($\\times %d$)', disp_scale), ...
        'Interpreter', 'latex', 'FontSize', opts.fontsize_title);
    pub_colorbar('$\phi$ [-]', opts);
    set(ax2, opts.axes_props{:});
    add_subfig_label(ax2, 'b', opts);
    
    %% (c) Crack-tip detail: stress near crack tip
    ax3 = nexttile;
    
    plot_ip_field_local(physics, mesh, "sh", "Internal", scale, opts.Pa2MPa);
    colormap(ax3, opts.cmap_diverge);
    hold on;
    
    % Overlay crack contour
    if ~isempty(Phi_grid)
        contour(X_grid, Y_grid, Phi_grid, [0.5 0.5], ...
            'LineColor', 'k', 'LineWidth', 1.2, 'LineStyle', '--');
    end
    
    % Zoom to crack tip region
    % Estimate crack tip location from phi field
    [~, idx] = max(Phi_grid(:));
    [iy, ix] = ind2sub(size(Phi_grid), idx);
    
    % Find rightmost point where phi > 0.5 along the midline
    mid_row = round(size(Phi_grid, 1) / 2);
    phi_mid = Phi_grid(mid_row, :);
    tip_cols = find(phi_mid > 0.5);
    if ~isempty(tip_cols)
        tip_x = X_grid(1, max(tip_cols));
        tip_y = Y_grid(mid_row, 1);
    else
        tip_x = Lx/2;
        tip_y = Ly/2;
    end
    
    % Zoom window around crack tip
    zoom_r = max(Lx, Ly) * 0.3;
    xlim([max(0, tip_x - zoom_r), min(Lx, tip_x + zoom_r)]);
    ylim([max(0, tip_y - zoom_r), min(Ly, tip_y + zoom_r)]);
    
    axis image;
    xlabel('$x$ [mm]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    ylabel('$y$ [mm]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    title('Crack-tip $\sigma_H$', ...
        'Interpreter', 'latex', 'FontSize', opts.fontsize_title);
    pub_colorbar('$\sigma_H$ [MPa]', opts);
    set(ax3, opts.axes_props{:});
    add_subfig_label(ax3, 'c', opts);
    
    export_fig(fig, 'fig04_crack_path', opts, 'both');
end


%% ===== Local helpers =====
function [Xg, Yg, Zg] = interpolate_nodal_field(physics, mesh, dofName, coord_scale)
    [dxTypes, dxSteps] = physics.dofSpace.getDofType({dofName});
    
    allNodes = mesh.GetAllNodesForGroup(1);
    dofs = physics.dofSpace.getDofIndices(dxTypes(1), allNodes);
    
    x_all = mesh.Nodes(allNodes, 1) * coord_scale;
    y_all = mesh.Nodes(allNodes, 2) * coord_scale;
    z_all = physics.StateVec{dxSteps(1)}(dofs);
    
    xlims = [min(x_all) max(x_all)];
    ylims = [min(y_all) max(y_all)];
    
    nx = 300; ny = 300;
    [Xg, Yg] = meshgrid(linspace(xlims(1), xlims(2), nx), ...
                         linspace(ylims(1), ylims(2), ny));
    try
        Zg = griddata(x_all, y_all, z_all, Xg, Yg, 'linear');
    catch
        Zg = [];
    end
end

function plot_nodal_field_local(physics, mesh, dofName, dispscale, plotloc, coord_scale)
    [dxTypes, dxSteps] = physics.dofSpace.getDofType({"dx";"dy";dofName});
    for g = 1:length(mesh.Elementgroups)
        if mesh.Elementgroups{g}.name == plotloc && mesh.Elementgroups{g}.type == "Q9"
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

function plot_ip_field_local(physics, mesh, varName, plotloc, coord_scale, val_scale)
    for g = 1:length(mesh.Elementgroups)
        if mesh.Elementgroups{g}.name == plotloc
            for e = 1:size(mesh.Elementgroups{g}.Elems, 1)
                IPvar = physics.Request_Info(varName, e, "Interior");
                xy = mesh.getIPCoords(g, e);
                meanx = mean(xy(1,:)); meany = mean(xy(2,:));
                for i = 1:size(xy,2)
                    xy(1,i) = xy(1,i) + sign(meanx - xy(1,i) + 1e-20)*1e-9;
                    xy(2,i) = xy(2,i) + sign(meany - xy(2,i) + 1e-20)*1e-9;
                end
                x(e,:) = xy(1,:) * coord_scale;
                y(e,:) = xy(2,:) * coord_scale;
                z(e,:) = IPvar * val_scale;
            end
            [xi,yi] = meshgrid(linspace(min(x(:)),max(x(:)),500), ...
                               linspace(min(y(:)),max(y(:)),500));
            zi = griddata(x(:),y(:),z(:),xi,yi,'linear');
            s = surf(xi,yi,zi); s.EdgeColor = 'none'; view(2);
            hold on;
        end
    end
end
