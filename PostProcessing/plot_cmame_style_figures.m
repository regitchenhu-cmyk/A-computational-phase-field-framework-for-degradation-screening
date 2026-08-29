function plot_cmame_style_figures(resultRoot, outdir, summaryCsv)
    %PLOT_CMAME_STYLE_FIGURES CMAME-inspired figure set.
    %
    % Design choices follow the supplied CMAME reference: compact multi-panel
    % figures, serif text, top colorbars, white backgrounds, thin axes, light
    % mesh overlays and local color scaling for nearly saturated fields.

    if nargin < 1 || isempty(resultRoot)
        error('plot_cmame_style_figures:resultRootRequired', ...
            'An explicit resultRoot is required; no legacy campaign fallback is allowed.');
    end
    if nargin < 2 || isempty(outdir)
        error('plot_cmame_style_figures:outdirRequired', ...
            'An explicit output directory is required.');
    end
    if nargin < 3 || isempty(summaryCsv)
        summaryCsv = fullfile(outdir, 'parametric_summary.csv');
    end
    if ~isfile(summaryCsv)
        error('plot_cmame_style_figures:missingSummary', ...
            'Canonical summary CSV not found: %s', summaryCsv);
    end
    if ~exist(outdir, 'dir'), mkdir(outdir); end

    addpath(genpath('./Models'))
    addpath(genpath('./Shapes'))

    set(0, 'DefaultFigureColor', 'w');
    set(0, 'DefaultAxesFontName', 'Times New Roman');
    set(0, 'DefaultTextFontName', 'Times New Roman');
    set(0, 'DefaultAxesFontSize', 8);

    plot_field_matrix(resultRoot, outdir);
    plot_response_curves(resultRoot, outdir);
    plot_mesh_lh_l9_summary(summaryCsv, outdir);
end

function plot_field_matrix(resultRoot, outdir)
    cases = { ...
        'Grid_G28x40', 'baseline'; ...
        'Orthogonal_O03_P21_T120_f100', 'thermal/frequency'; ...
        'Orthogonal_O07_P35_T080_f100', 'pressure/frequency'};
    fields = {'phi', 'CL', 'alpha', 'Df'};
    fieldTitles = {'$\phi$', '$C_L$', '$\alpha_a$', '$D_f$'};

    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 17.6 13.2]);
    tl = tiledlayout(size(cases,1), numel(fields), 'TileSpacing', 'compact', 'Padding', 'compact');

    for r = 1:size(cases,1)
        data = load(fullfile(resultRoot, cases{r,1}, 'end.mat'), 'mesh', 'physics');
        F = nodal_fields(data);
        for c = 1:numel(fields)
            ax = nexttile(tl);
            values = F.(fields{c});
            draw_nodal_field(ax, data.mesh, values, cmap_blue_white_red(256));
            add_structured_mesh_overlay(ax, data.mesh, [0.35 0.35 0.35], 0.16);
            axis(ax, 'equal');
            axis(ax, 'tight');
            ax.XTick = [];
            ax.YTick = [];
            ax.Box = 'on';
            ax.LineWidth = 0.45;
            if r == 1
                title(ax, fieldTitles{c}, 'Interpreter', 'latex', 'FontSize', 8.5);
            end
            if c == 1
                text(ax, -0.15, 0.5, sprintf('(%c) %s', char('a'+r-1), cases{r,2}), ...
                    'Units', 'normalized', 'Rotation', 90, 'HorizontalAlignment', 'center', ...
                    'FontSize', 7.4);
            end
            cb = colorbar(ax, 'northoutside');
            cb.Box = 'off';
            cb.FontSize = 6.5;
            cb.TickDirection = 'out';
            set_local_clim(ax, values);
        end
    end

    export_cmame(fig, outdir, 'cmame_field_matrix');
end

function plot_response_curves(resultRoot, outdir)
    cases = { ...
        'Grid_G28x40', 'baseline', [0.00 0.26 0.62]; ...
        'Orthogonal_O03_P21_T120_f100', 'T/f severe', [0.82 0.28 0.06]; ...
        'Orthogonal_O07_P35_T080_f100', 'P/f severe', [0.55 0.20 0.64]};

    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 17.6 7.2]);
    tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    ax1 = nexttile(tl); hold(ax1, 'on');
    ax2 = nexttile(tl); hold(ax2, 'on');
    for k = 1:size(cases,1)
        data = load(fullfile(resultRoot, cases{k,1}, 'end.mat'), 'results', 'tvec', 'mesh');
        t = data.tvec(:) / 86400;
        c = cases{k,3};
        [screeningIndex, ~] = compute_screening_index(data);
        fatigue = data.results.Df_p99(:);
        fluid = data.results.CL_avg(:);

        markerIdx = unique(round(linspace(1, numel(t), 7)));
        plot(ax1, t, fluid, '-', 'Color', c, 'LineWidth', 1.25, ...
            'Marker', marker_for(k), 'MarkerIndices', markerIdx, ...
            'MarkerSize', 4.0, 'MarkerFaceColor', 'w');
        plot(ax1, t, fatigue, linestyle_for(k), 'Color', c, 'LineWidth', 1.10, ...
            'Marker', marker_for(k+3), 'MarkerIndices', markerIdx, ...
            'MarkerSize', 3.8, 'MarkerFaceColor', c);
        plot(ax2, t, screeningIndex, linestyle_for(k), 'Color', c, 'LineWidth', 1.45, ...
            'Marker', marker_for(k), 'MarkerIndices', markerIdx, ...
            'MarkerSize', 4.2, 'MarkerFaceColor', 'w', 'DisplayName', cases{k,2});
    end
    format_curve_axes(ax1);
    xlabel(ax1, 'Time [days]');
    ylabel(ax1, 'Component [-]');
    title(ax1, '(a) Fluid (solid) and fatigue (dashed)');
    legend(ax1, {'baseline $C_L$', 'baseline $D_f$', 'T/f $C_L$', 'T/f $D_f$', ...
        'P/f $C_L$', 'P/f $D_f$'}, 'Interpreter', 'latex', 'Box', 'off', ...
        'Location', 'northwest', 'FontSize', 7);

    format_curve_axes(ax2);
    xlabel(ax2, 'Time [days]');
    ylabel(ax2, 'Screening index I_S [-]');
    title(ax2, '(b) Percentile screening envelope');
    legend(ax2, 'Box', 'off', 'Location', 'northwest', 'FontSize', 7);
    ylim(ax2, [0 1.03]);

    export_cmame(fig, outdir, 'cmame_response_curves');
end

function plot_mesh_lh_l9_summary(summaryCsv, outdir)
    % Endpoint panels consume only the explicitly supplied canonical CSV.
    T = readtable(summaryCsv, 'TextType', 'string');

    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 17.6 10.5]);
    tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    S = T(startsWith(T.Label, "Grid_"), :);
    [nx, ny] = parse_grid(S.Label);
    ne = nx .* ny;
    ax = nexttile(tl); hold(ax, 'on');
    plot(ax, ne, S.ScreeningIndex, '-o', 'Color', [0 0.25 0.60], ...
        'MarkerFaceColor', 'w', 'MarkerEdgeColor', [0 0.25 0.60], ...
        'LineWidth', 1.25, 'MarkerSize', 4.2);
    plot(ax, ne, S.FatigueP99, '--s', 'Color', [0.80 0.32 0.00], ...
        'MarkerFaceColor', 'w', 'MarkerEdgeColor', [0.80 0.32 0.00], ...
        'LineWidth', 1.20, 'MarkerSize', 4.2);
    format_curve_axes(ax);
    xlabel(ax, 'Number of elements');
    ylabel(ax, 'Index [-]');
    title(ax, '(a) Mesh sensitivity');
    legend(ax, {'$I_S$', '$D_{f,P99}$'}, 'Interpreter', 'latex', 'Box', 'off', 'Location', 'northwest', 'FontSize', 7);

    S = T(startsWith(T.Label, "LengthScale_"), :);
    lh = parse_lh(S.Label);
    [lh, order] = sort(lh); S = S(order,:);
    ax = nexttile(tl); hold(ax, 'on');
    plot(ax, lh, 100*(S.ScreeningIndex/S.ScreeningIndex(2)-1), ...
        '-.^', 'Color', [0.45 0.22 0.62], 'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', [0.45 0.22 0.62], 'LineWidth', 1.25, 'MarkerSize', 4.2);
    yline(ax, 0, ':', 'Color', [0.25 0.25 0.25], 'LineWidth', 0.7);
    format_curve_axes(ax);
    xlabel(ax, '$l/h$', 'Interpreter', 'latex');
    ylabel(ax, 'Change from $l/h=3$ [\%]', 'Interpreter', 'latex');
    title(ax, '(b) Length-scale sensitivity');

    O = T(startsWith(T.Label, "Orthogonal_"), :);
    ax = nexttile(tl); hold(ax, 'on');
    bar(ax, 1:height(O), O.ScreeningIndex, 0.72, ...
        'FaceColor', [0.18 0.44 0.70], 'EdgeColor', 'none');
    format_curve_axes(ax);
    xlabel(ax, 'L9 case');
    ylabel(ax, 'Screening index I_S [-]');
    title(ax, '(c) Orthogonal cases');
    xlim(ax, [0.4 height(O)+0.6]);
    ylim(ax, [0 1.05]);

    [P, Temp, Freq] = parse_orthogonal(O.Label);
    ax = nexttile(tl); hold(ax, 'on');
    main_effect(ax, P, O.FatigueP99, [0 0.25 0.60], '-', 'o');
    main_effect(ax, Temp, O.FatigueP99, [0.00 0.48 0.22], '--', 's');
    main_effect(ax, Freq, O.FatigueP99, [0.80 0.32 0.00], '-.', '^');
    format_curve_axes(ax);
    xlabel(ax, 'Factor level');
    ylabel(ax, 'Mean $D_{f,P99}$ [-]', 'Interpreter', 'latex');
    title(ax, '(d) Descriptive fatigue level means');
    legend(ax, {'Pressure', 'Temperature', 'Frequency'}, 'Box', 'off', 'Location', 'northwest', 'FontSize', 7);

    export_cmame(fig, outdir, 'cmame_parametric_summary');
end

function F = nodal_fields(data)
    mesh = data.mesh;
    physics = data.physics;
    allNodes = mesh.GetAllNodesForGroup(mesh.getGroupIndex("Internal"));
    mdl = physics.models{1};

    phiDofs = mdl.dofSpace.getDofIndices(mdl.dofTypeIndices(3), allNodes);
    clDofs = mdl.dofSpace.getDofIndices(mdl.dofTypeIndices(4), allNodes);
    aDofs = mdl.dofSpace.getDofIndices(mdl.dofTypeIndices(5), allNodes);
    dfDofs = mdl.dofSpace.getDofIndices(mdl.dofTypeIndices(6), allNodes);

    F.nodes = allNodes;
    F.phi = physics.StateVec{mdl.phi_step}(phiDofs);
    F.CL = physics.StateVec{mdl.CL_step}(clDofs);
    F.alpha = physics.StateVec{mdl.aging_step}(aDofs);
    F.Df = physics.StateVec{mdl.fatigue_step}(dfDofs);
end

function draw_nodal_field(ax, mesh, nodalValues, cmap)
    x = mesh.Nodes(:,1) * 1e3;
    y = mesh.Nodes(:,2) * 1e3;
    xq = linspace(min(x), max(x), 260);
    yq = linspace(min(y), max(y), 360);
    [Xq, Yq] = meshgrid(xq, yq);
    Vq = griddata(x, y, nodalValues(:), Xq, Yq, 'natural');
    imagesc(ax, xq, yq, Vq);
    set(ax, 'YDir', 'normal');
    axis(ax, 'image');
    colormap(ax, cmap);
end

function add_structured_mesh_overlay(ax, mesh, color, alphaVal)
    x = unique(round(mesh.Nodes(:,1) * 1e3, 10));
    y = unique(round(mesh.Nodes(:,2) * 1e3, 10));
    x = x(1:2:end);
    y = y(1:2:end);
    lineColor = (1-alphaVal)*[1 1 1] + alphaVal*color;
    hold(ax, 'on');
    for i = 1:numel(x)
        plot(ax, [x(i) x(i)], [min(y) max(y)], '-', 'Color', lineColor, 'LineWidth', 0.18);
    end
    for j = 1:numel(y)
        plot(ax, [min(x) max(x)], [y(j) y(j)], '-', 'Color', lineColor, 'LineWidth', 0.18);
    end
end

function set_local_clim(ax, values)
    lo = min(values(:));
    hi = max(values(:));
    if abs(hi - lo) < 1e-10
        lo = lo - 1e-3;
        hi = hi + 1e-3;
    end
    pad = 0.03 * max(abs(hi-lo), eps);
    clim(ax, [lo-pad, hi+pad]);
end

function format_curve_axes(ax)
    ax.Box = 'on';
    ax.LineWidth = 0.55;
    ax.TickDir = 'out';
    ax.FontSize = 8;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.GridAlpha = 0.12;
    ax.Layer = 'top';
end

function export_cmame(fig, outdir, name)
    exportgraphics(fig, fullfile(outdir, name + ".png"), 'Resolution', 600);
    exportgraphics(fig, fullfile(outdir, name + ".eps"), 'ContentType', 'vector');
    close(fig);
end

function cmap = cmap_blue_white_red(n)
    if nargin < 1, n = 256; end
    n1 = floor(n/2);
    n2 = n - n1;
    blue = [0.12 0.24 0.62];
    white = [0.96 0.96 0.96];
    red = [0.72 0.05 0.05];
    cmap = [interp1([1 n1], [blue; white], 1:n1); ...
            interp1([1 n2], [white; red], 1:n2)];
end

function [nx, ny] = parse_grid(labels)
    nx = zeros(numel(labels), 1); ny = nx;
    for i = 1:numel(labels)
        tok = regexp(labels(i), 'G(\d+)x(\d+)', 'tokens', 'once');
        nx(i) = str2double(tok{1});
        ny(i) = str2double(tok{2});
    end
end

function lh = parse_lh(labels)
    lh = zeros(numel(labels), 1);
    for i = 1:numel(labels)
        tok = regexp(labels(i), 'lh([\d.]+)', 'tokens', 'once');
        lh(i) = str2double(tok{1});
    end
end

function [P, T, f] = parse_orthogonal(labels)
    P = zeros(numel(labels), 1); T = P; f = P;
    for i = 1:numel(labels)
        tok = regexp(labels(i), 'P(\d+)_T(\d+)_f(\d+)', 'tokens', 'once');
        P(i) = str2double(tok{1});
        T(i) = str2double(tok{2});
        f(i) = str2double(tok{3})/1000;
    end
end

function main_effect(ax, x, y, color, lineStyle, marker)
    levels = unique(x);
    mu = zeros(size(levels));
    for i = 1:numel(levels)
        mu(i) = mean(y(x == levels(i)), 'omitnan');
    end
    plot(ax, 1:numel(levels), mu, lineStyle, 'Color', color, ...
        'Marker', marker, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', color, ...
        'LineWidth', 1.25, 'MarkerSize', 4.2);
    ax.XTick = 1:numel(levels);
    ax.XTickLabel = ["low", "mid", "high"];
end

function ls = linestyle_for(k)
    styles = {'-', '--', '-.', ':'};
    ls = styles{mod(k-1, numel(styles)) + 1};
end

function mk = marker_for(k)
    markers = {'o', 's', '^', 'd', 'v', '>'};
    mk = markers{mod(k-1, numel(markers)) + 1};
end
