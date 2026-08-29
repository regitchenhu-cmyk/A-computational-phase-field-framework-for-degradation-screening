function generate_animation(result_folder, field_name, varargin)
    %GENERATE_ANIMATION Create MP4 animation of field evolution.
    %
    % Usage:
    %   generate_animation('./Results_Seal/Standard_21MPa_80C/', 'phi')
    %   generate_animation('./Results_Seal/Standard_21MPa_80C/', 'CL', ...
    %       'FrameRate', 15, 'StepInterval', 5)
    %
    % Optional parameters:
    %   'FrameRate'     - Video frame rate (default: 10)
    %   'StepInterval'  - Plot every N-th step (default: 2)
    %   'DispScale'     - Displacement magnification (default: 0)
    %   'CLimits'       - Color axis limits [min max] (default: auto)
    %   'Colormap'      - Colormap name or matrix (default: auto)

    addpath(genpath('./Models'))
    addpath(genpath('./Shapes'))
    addpath(genpath('./PostProcessing'))
    
    p = inputParser;
    addParameter(p, 'FrameRate', 10);
    addParameter(p, 'StepInterval', 2);
    addParameter(p, 'DispScale', 0);
    addParameter(p, 'CLimits', []);
    addParameter(p, 'Colormap', []);
    parse(p, varargin{:});
    
    opts = pub_style();
    scale = opts.m2mm;
    
    % Detect available steps (numeric-named .mat files only)
    files = dir(fullfile(result_folder, '*.mat'));
    steps = [];
    for i = 1:length(files)
        [~, fname, ~] = fileparts(files(i).name);
        num = str2double(fname);
        if ~isnan(num); steps = [steps, num]; end
    end
    if isempty(steps)
        error('No numeric step .mat files found for animation in: %s\nAnimation requires multiple saved steps (e.g. 10.mat, 20.mat, ...).', result_folder);
    end
    steps = sort(steps);
    steps = steps(1:p.Results.StepInterval:end);
    
    fprintf('Generating animation: %d frames\n', length(steps));
    
    % Determine colormap
    switch field_name
        case 'phi'
            cmap = flipud(hot(256));
            clabel = '\phi';
            if isempty(p.Results.CLimits); clims = [0 1]; else; clims = p.Results.CLimits; end
        case 'CL'
            cmap = parula(256);
            clabel = 'C_L';
            if isempty(p.Results.CLimits); clims = []; else; clims = p.Results.CLimits; end
        case 'alpha_a'
            cmap = copper(256);
            clabel = '\alpha_{age}';
            if isempty(p.Results.CLimits); clims = [0 1]; else; clims = p.Results.CLimits; end
        case 'D_f'
            cmap = cool(256);
            clabel = 'D_f';
            if isempty(p.Results.CLimits); clims = [0 1]; else; clims = p.Results.CLimits; end
        otherwise
            cmap = jet(256);
            clabel = field_name;
            clims = p.Results.CLimits;
    end
    
    if ~isempty(p.Results.Colormap)
        if ischar(p.Results.Colormap) || isstring(p.Results.Colormap)
            cmap = feval(p.Results.Colormap, 256);
        else
            cmap = p.Results.Colormap;
        end
    end
    
    % Setup video
    outdir = './Animations';
    if ~exist(outdir, 'dir'); mkdir(outdir); end
    
    vidname = fullfile(outdir, [field_name '_evolution.mp4']);
    vid = VideoWriter(vidname, 'MPEG-4');
    vid.FrameRate = p.Results.FrameRate;
    vid.Quality = 95;
    open(vid);
    
    fig = figure('Visible', 'off');
    set(fig, 'Position', [100 100 800 600]);
    set(fig, 'Color', 'w');
    
    for i = 1:length(steps)
        fname = fullfile(result_folder, [num2str(steps(i)) '.mat']);
        data = load(fname);
        
        clf(fig);
        
        % Main field plot
        ax1 = axes('Position', [0.08 0.25 0.60 0.65]);
        plot_nodal_anim(data.physics, data.mesh, field_name, ...
            p.Results.DispScale, "Internal", scale);
        colormap(ax1, cmap);
        if ~isempty(clims); caxis(clims); end
        axis image;
        
        cb = colorbar;
        cb.FontName = opts.fontname;
        cb.FontSize = opts.fontsize_tick;
        
        xlabel('$x$ [mm]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
        ylabel('$y$ [mm]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
        set(ax1, opts.axes_props{:});
        
        % Time annotation
        t_val = data.tvec(end);
        if t_val > 86400*30
            t_str = sprintf('t = %.1f days', t_val * opts.s2day);
        elseif t_val > 3600
            t_str = sprintf('t = %.1f hours', t_val * opts.s2hr);
        else
            t_str = sprintf('t = %.0f s', t_val);
        end
        title(sprintf('%s    (Step %d,  %s)', clabel, steps(i), t_str), ...
            'FontName', opts.fontname, 'FontSize', opts.fontsize_title);
        
        % Inset: crack length over time
        ax2 = axes('Position', [0.74 0.55 0.24 0.35]);
        plot(data.tvec * opts.s2day, data.results.LFrac * opts.m2mm, '-', ...
            'Color', opts.colors(1,:), 'LineWidth', 1.2);
        hold on;
        plot(t_val * opts.s2day, data.results.LFrac(end) * opts.m2mm, 'o', ...
            'MarkerFaceColor', opts.colors(1,:), 'MarkerEdgeColor', 'w', ...
            'MarkerSize', 6);
        xlabel('Days', 'FontSize', 9);
        ylabel('L_{frac} [mm]', 'FontSize', 9);
        set(ax2, 'FontSize', 8, 'Box', 'on', 'LineWidth', 0.4);
        
        % Inset: percentile screening index over time
        ax3 = axes('Position', [0.74 0.10 0.24 0.35]);
        FI = compute_failure_index(data);
        plot(data.tvec * opts.s2day, FI, '-', ...
            'Color', [0.1 0.1 0.1], 'LineWidth', 1.2);
        hold on;
        plot(t_val * opts.s2day, FI(end), 'o', ...
            'MarkerFaceColor', [0.1 0.1 0.1], 'MarkerEdgeColor', 'w', ...
            'MarkerSize', 6);
        xlabel('Days', 'FontSize', 9);
        ylabel('I_S [-]', 'FontSize', 9);
        set(ax3, 'FontSize', 8, 'Box', 'on', 'LineWidth', 0.4);
        ylim([0 max(1, max(FI)*1.1)]);
        
        drawnow;
        frame = getframe(fig);
        writeVideo(vid, frame);
        
        if mod(i, 10) == 0
            fprintf('  Frame %d / %d\n', i, length(steps));
        end
    end
    
    close(vid);
    close(fig);
    fprintf('Animation saved: %s\n', vidname);
end


function plot_nodal_anim(physics, mesh, dofName, dispscale, plotloc, coord_scale)
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
