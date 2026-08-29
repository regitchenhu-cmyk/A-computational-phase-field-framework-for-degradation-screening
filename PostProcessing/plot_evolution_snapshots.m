function plot_evolution_snapshots(result_folder, steps, opts, outdir)
    %PLOT_EVOLUTION_SNAPSHOTS Multi-panel temporal evolution of damage fields.
    %
    % Generates a figure with N columns × 3 rows:
    %   Row 1: Phase-field phi at each time step
    %   Row 2: Fluid concentration CL at each time step
    %   Row 3: Aging variable alpha_a at each time step

    nsteps = length(steps);
    scale = opts.m2mm;
    
    %% ===== Figure A: Phase-field evolution =====
    fig1 = figure('Visible', 'on');
    set(fig1, 'Position', [50 50 280*nsteps 750]);
    set(fig1, 'Color', 'w');
    
    tl = tiledlayout(3, nsteps, 'TileSpacing', 'tight', 'Padding', 'compact');
    
    for i = 1:nsteps
        fname = fullfile(result_folder, [num2str(steps(i)) '.mat']);
        if ~exist(fname, 'file')
            warning('File not found: %s, skipping.', fname);
            continue;
        end
        data = load(fname);
        physics = data.physics;
        mesh = data.mesh;
        
        t_val = data.tvec(end);
        if t_val > 86400
            t_str = sprintf('$t = %.1f$ d', t_val * opts.s2day);
        elseif t_val > 3600
            t_str = sprintf('$t = %.1f$ h', t_val * opts.s2hr);
        else
            t_str = sprintf('$t = %.0f$ s', t_val);
        end
        
        %% Row 1: Phase-field
        ax = nexttile(i);
        plot_nodal_pub(physics, mesh, "phi", 0, "Internal", scale);
        colormap(ax, flipud(hot(256)));
        caxis([0 1]);
        axis image; axis off;
        title(t_str, 'Interpreter', 'latex', 'FontSize', opts.fontsize_title);
        if i == 1
            ylabel('$\phi$', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label, ...
                'Rotation', 0, 'HorizontalAlignment', 'right');
        end
        if i == nsteps
            pub_colorbar('', opts);
        end
        
        %% Row 2: Fluid concentration
        ax = nexttile(nsteps + i);
        plot_nodal_pub(physics, mesh, "CL", 0, "Internal", scale);
        colormap(ax, parula(256));
        c_max = max(0.01, max(get_nodal_vals(physics, "CL")));
        caxis([0 c_max]);
        axis image; axis off;
        if i == 1
            ylabel('$C_L$', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label, ...
                'Rotation', 0, 'HorizontalAlignment', 'right');
        end
        if i == nsteps
            pub_colorbar('', opts);
        end
        
        %% Row 3: Aging
        ax = nexttile(2*nsteps + i);
        plot_nodal_pub(physics, mesh, "alpha_a", 0, "Internal", scale);
        colormap(ax, copper(256));
        a_max = max(0.01, max(get_nodal_vals(physics, "alpha_a")));
        caxis([0 a_max]);
        axis image; axis off;
        if i == 1
            ylabel('$\alpha$', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label, ...
                'Rotation', 0, 'HorizontalAlignment', 'right');
        end
        if i == nsteps
            pub_colorbar('', opts);
        end
    end
    
    export_fig(fig1, 'fig03_evolution_snapshots', opts, 'both');
    
    
    %% ===== Figure B: Fatigue + aging evolution (2 rows) =====
    fig2 = figure('Visible', 'on');
    set(fig2, 'Position', [50 50 280*nsteps 520]);
    set(fig2, 'Color', 'w');
    
    tl = tiledlayout(2, nsteps, 'TileSpacing', 'tight', 'Padding', 'compact');
    
    for i = 1:nsteps
        fname = fullfile(result_folder, [num2str(steps(i)) '.mat']);
        if ~exist(fname, 'file'); continue; end
        data = load(fname);
        physics = data.physics;
        mesh = data.mesh;
        
        t_val = data.tvec(end);
        if t_val > 86400
            t_str = sprintf('$t = %.1f$ d', t_val * opts.s2day);
        else
            t_str = sprintf('$t = %.1f$ h', t_val * opts.s2hr);
        end
        
        %% Row 1: Fatigue damage
        ax = nexttile(i);
        plot_nodal_pub(physics, mesh, "D_f", 0, "Internal", scale);
        colormap(ax, cool(256));
        df_max = max(0.01, max(get_nodal_vals(physics, "D_f")));
        caxis([0 df_max]);
        axis image; axis off;
        title(t_str, 'Interpreter', 'latex', 'FontSize', opts.fontsize_title);
        if i == 1
            ylabel('$D_f$', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label, ...
                'Rotation', 0, 'HorizontalAlignment', 'right');
        end
        if i == nsteps
            pub_colorbar('', opts);
        end
        
        %% Row 2: Deformed phi
        ax = nexttile(nsteps + i);
        plot_nodal_pub(physics, mesh, "phi", 50, "Internal", scale);
        colormap(ax, flipud(hot(256)));
        caxis([0 1]);
        axis image; axis off;
        if i == 1
            ylabel('$\phi$(def.)', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label, ...
                'Rotation', 0, 'HorizontalAlignment', 'right');
        end
        if i == nsteps
            pub_colorbar('', opts);
        end
    end
    
    export_fig(fig2, 'fig03b_fatigue_evolution', opts, 'both');
end


%% ===== Local helpers =====
function plot_nodal_pub(physics, mesh, dofName, dispscale, plotloc, coord_scale)
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

function vals = get_nodal_vals(physics, dofName)
    [dxTypes, dxSteps] = physics.dofSpace.getDofType({dofName});
    vals = physics.StateVec{dxSteps(1)};
end
