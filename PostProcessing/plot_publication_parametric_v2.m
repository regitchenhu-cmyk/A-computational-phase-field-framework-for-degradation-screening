function plot_publication_parametric_v2(summaryCsv, resultRoot, outdir)
    %PLOT_PUBLICATION_PARAMETRIC_V2 Publication-grade parametric figures.
    %
    % These figures report a percentile screening index rather than an
    % experimentally calibrated failure threshold.

    if nargin < 1 || isempty(summaryCsv)
        error('plot_publication_parametric_v2:summaryRequired', ...
            'An explicit canonical summary CSV is required; legacy fallback is disabled.');
    end
    if nargin < 2 || isempty(resultRoot)
        error('plot_publication_parametric_v2:resultRootRequired', ...
            'An explicit resultRoot is required; legacy fallback is disabled.');
    end
    if nargin < 3 || isempty(outdir)
        outdir = fileparts(summaryCsv);
    end
    if ~isfile(summaryCsv)
        error('plot_publication_parametric_v2:missingSummary', ...
            'Canonical summary CSV not found: %s', summaryCsv);
    end
    if ~exist(outdir, 'dir'), mkdir(outdir); end

    addpath(genpath('./Models'))
    addpath(genpath('./Shapes'))

    T = readtable(summaryCsv, 'TextType', 'string');
    set(0, 'DefaultFigureColor', 'w');
    set(0, 'DefaultAxesFontName', 'Times New Roman');
    set(0, 'DefaultTextFontName', 'Times New Roman');

    plot_baseline_history(resultRoot, outdir);
    plot_mesh(T, outdir);
    plot_length_scale(T, outdir);
    plot_l9(T, outdir);
    plot_fatigue_localization(resultRoot, outdir);
end

function plot_baseline_history(resultRoot, outdir)
    C = palette();
    data = load(fullfile(resultRoot, 'Grid_G28x40', 'end.mat'), 'results', 'tvec');
    t = data.tvec(:) / 86400;
    res = data.results;
    fi = res.screening_index(:);

    fig = figure('Visible', 'off', 'Position', [80 80 1000 380]);
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    plot(t, res.CL_avg, '-', 'Color', C.blue, 'LineWidth', 1.8);
    hold on
    plot(t, res.alpha_p99, '-', 'Color', C.green, 'LineWidth', 1.8);
    plot(t, res.Df_p99, '-', 'Color', C.orange, 'LineWidth', 1.8);
    plot(t, res.phi_p99, '-', 'Color', C.red, 'LineWidth', 1.8);
    xlabel('Time [days]');
    ylabel('Normalized component [-]');
    legend({'Fluid average (diagnostic)','Aging P99','Fatigue P99','Phase field P99'}, ...
        'Box', 'off', 'Location', 'northwest', 'FontSize', 7.8);
    title('(a) Baseline degradation history');
    style_axes(gca);

    nexttile;
    plot(t, fi, '-', 'Color', C.purple, 'LineWidth', 2.0);
    xlabel('Time [days]');
    ylabel('Screening index I_S [-]');
    ylim([0 max(0.25, 1.08*max(fi))]);
    title('(b) Percentile screening envelope');
    style_axes(gca);

    export_pub(fig, outdir, "publication_baseline_history");
    close(fig);
end

function C = palette()
    C.blue   = [0.000, 0.278, 0.671];
    C.orange = [0.835, 0.369, 0.000];
    C.green  = [0.000, 0.502, 0.235];
    C.purple = [0.459, 0.306, 0.694];
    C.gray   = [0.320, 0.320, 0.320];
    C.light  = [0.910, 0.930, 0.950];
    C.red    = [0.760, 0.110, 0.160];
end

function style_axes(ax)
    ax.Box = 'on';
    ax.LineWidth = 0.8;
    ax.TickDir = 'out';
    ax.TickLength = [0.012 0.012];
    ax.FontSize = 8.5;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.GridAlpha = 0.18;
    ax.MinorGridAlpha = 0.08;
    ax.Layer = 'top';
end

function export_pub(fig, outdir, name)
    set(fig, 'Units', 'centimeters');
    exportgraphics(fig, fullfile(outdir, name + ".png"), 'Resolution', 600);
    exportgraphics(fig, fullfile(outdir, name + ".eps"), 'ContentType', 'vector');
end

function plot_mesh(T, outdir)
    C = palette();
    S = T(startsWith(T.Label, "Grid_"), :);
    [nx, ny] = parse_grid_labels(S.Label);
    ne = nx .* ny;

    fig = figure('Visible', 'off', 'Position', [80 80 1100 360]);
    tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    plot(ne, S.ScreeningIndex, '-o', 'Color', C.blue, ...
        'MarkerFaceColor', C.blue, 'LineWidth', 1.8, 'MarkerSize', 5);
    xlabel('Number of elements');
    ylabel('Screening index I_S [-]');
    legend({'I_S'}, ...
        'Location', 'northwest', 'Box', 'off', 'FontSize', 7.8);
    title('(a) Mesh response');
    style_axes(gca);

    nexttile;
    plot(ne, S.FluidAvg, '-o', 'Color', C.blue, 'MarkerFaceColor', C.blue, 'LineWidth', 1.6);
    hold on
    plot(ne, S.AgingP99, '-s', 'Color', C.green, 'MarkerFaceColor', C.green, 'LineWidth', 1.6);
    plot(ne, S.FatigueP99, '-^', 'Color', C.orange, 'MarkerFaceColor', C.orange, 'LineWidth', 1.6);
    xlabel('Number of elements');
    ylabel('Component value [-]');
    ylim([0.06 0.25]);
    legend({'Fluid average (diagnostic)','Aging P99','Fatigue P99'}, ...
        'Location', 'northwest', 'Box', 'off', 'FontSize', 7.8);
    title('(b) Converged components');
    style_axes(gca);

    nexttile;
    yyaxis left
    plot(ne, S.FatigueMax, '--o', 'Color', C.red, 'MarkerFaceColor', 'w', 'LineWidth', 1.3);
    hold on
    plot(ne, S.FatigueP99, '-o', 'Color', C.orange, 'MarkerFaceColor', C.orange, 'LineWidth', 1.8);
    ylabel('Fatigue damage [-]');
    ylim([0 1.05]);
    yyaxis right
    bar(ne, S.FatigueFracGt09IP, 0.22, 'FaceColor', C.gray, 'EdgeColor', 'none', 'FaceAlpha', 0.45);
    ylabel('Fraction D_f > 0.9 [-]');
    ylim([0 2.5e-4]);
    xlabel('Number of elements');
    title('(c) Corner-localization check');
    legend({'D_{f,max}','D_{f,P99}','area fraction'}, ...
        'Location', 'northwest', 'Box', 'off', 'FontSize', 7.6);
    style_axes(gca);

    export_pub(fig, outdir, "publication_mesh_independence");
    close(fig);
end

function plot_length_scale(T, outdir)
    C = palette();
    S = T(startsWith(T.Label, "LengthScale_"), :);
    lh = zeros(height(S), 1);
    for i = 1:height(S)
        tok = regexp(S.Label(i), 'lh([\d.]+)', 'tokens', 'once');
        lh(i) = str2double(tok{1});
    end
    [lh, order] = sort(lh);
    S = S(order, :);

    fig = figure('Visible', 'off', 'Position', [80 80 1000 360]);
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    plot(lh, S.ScreeningIndex, '-o', 'Color', C.blue, ...
        'MarkerFaceColor', C.blue, 'LineWidth', 1.8, 'MarkerSize', 5);
    hold on
    plot(lh, S.FluidAvg, '-s', 'Color', C.green, ...
        'MarkerFaceColor', C.green, 'LineWidth', 1.4, 'MarkerSize', 5);
    xlabel('Regularization ratio l/h [-]');
    ylabel('Final index [-]');
    legend({'I_S','Fluid average (diagnostic)'}, 'Box', 'off', 'Location', 'northeast', 'FontSize', 7.8);
    title('(a) Length-scale response');
    style_axes(gca);

    nexttile;
    base = S.ScreeningIndex(2);
    rel = 100 * (S.ScreeningIndex - base) ./ base;
    plot(lh, rel, '-o', 'Color', C.purple, 'MarkerFaceColor', C.purple, ...
        'LineWidth', 1.8, 'MarkerSize', 5);
    yline(0, '-', 'Color', C.gray, 'LineWidth', 0.8);
    xlabel('Regularization ratio l/h [-]');
    ylabel('Change from l/h = 3 [%]');
    title('(b) Relative variation');
    style_axes(gca);

    export_pub(fig, outdir, "publication_length_scale_sensitivity");
    close(fig);
end

function plot_l9(T, outdir)
    C = palette();
    O = T(startsWith(T.Label, "Orthogonal_"), :);
    [P, Temp, Freq] = parse_orthogonal_labels(O.Label);

    fig = figure('Visible', 'off', 'Position', [80 80 1150 720]);
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    x = 1:height(O);
    b = bar(x, O.ScreeningIndex, 0.72, 'FaceColor', 'flat', 'EdgeColor', 'none');
    for i = 1:numel(x)
        if P(i) == 21, b.CData(i,:) = C.blue; end
        if P(i) == 28, b.CData(i,:) = C.green; end
        if P(i) == 35, b.CData(i,:) = C.orange; end
    end
    ylim([0 1.05]);
    xlim([0.4 height(O)+0.6]);
    xlabel('L9 case');
    ylabel('Screening index I_S [-]');
    title('(a) Case-wise response');
    style_axes(gca);

    nexttile;
    plot_main_effect(P, O.ScreeningIndex, C.blue, 'Pressure [MPa]');
    hold on
    plot_main_effect(Temp, O.ScreeningIndex, C.green, 'Temperature [degC]');
    plot_main_effect(Freq*1000, O.ScreeningIndex, C.orange, 'Frequency [mHz]');
    ylabel('Mean I_S [-]');
    title('(b) Main effects on I_S');
    legend({'Pressure','Temperature','Frequency'}, 'Box', 'off', 'Location', 'southeast', 'FontSize', 7.8);
    style_axes(gca);

    nexttile;
    plot_main_effect(P, O.FatigueP99, C.blue, 'Pressure [MPa]');
    hold on
    plot_main_effect(Temp, O.FatigueP99, C.green, 'Temperature [degC]');
    plot_main_effect(Freq*1000, O.FatigueP99, C.orange, 'Frequency [mHz]');
    ylabel('Mean D_{f,P99} [-]');
    title('(c) Main effects on fatigue');
    style_axes(gca);

    nexttile;
    plot_main_effect(P, O.FluidAvg, C.blue, 'Pressure [MPa]');
    hold on
    plot_main_effect(Temp, O.FluidAvg, C.green, 'Temperature [degC]');
    plot_main_effect(Freq*1000, O.FluidAvg, C.orange, 'Frequency [mHz]');
    ylabel('Mean fluid uptake [-]');
    title('(d) Main effects on fluid uptake');
    style_axes(gca);

    export_pub(fig, outdir, "publication_L9_orthogonal_effects");
    close(fig);
end

function plot_main_effect(x, y, color, ~)
    levels = unique(x);
    mu = zeros(size(levels));
    for k = 1:numel(levels)
        mu(k) = mean(y(x == levels(k)), 'omitnan');
    end
    xn = 1:numel(levels);
    plot(xn, mu, '-o', 'Color', color, 'MarkerFaceColor', color, ...
        'LineWidth', 1.6, 'MarkerSize', 4.8);
    set(gca, 'XTick', xn, 'XTickLabel', ["low","mid","high"]);
    xlabel('Factor level');
end

function plot_fatigue_localization(resultRoot, outdir)
    C = palette();
    data = load(fullfile(resultRoot, 'Grid_G40x56', 'end.mat'), 'physics', 'mesh');
    Df = data.physics.models{8}.Df_ip;
    vals = zeros(numel(Df), 1);
    g = data.mesh.getGroupIndex("Internal");
    nElem = size(Df, 1);
    nIp = size(Df, 2);
    xy = zeros(nElem*nIp, 2);
    row = 0;
    for e = 1:nElem
        eip = data.mesh.getIPCoords(g, e);
        if size(eip, 2) ~= 2 && size(eip, 1) == 2
            eip = eip.';
        end
        for ip = 1:nIp
            row = row + 1;
            xy(row,:) = eip(ip,:);
            vals(row) = Df(e, ip);
        end
    end
    p99 = prctile(vals, 99);

    fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 17.6 5.6], 'Color', 'w');
    tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    histogram(vals, 70, 'FaceColor', C.blue, 'EdgeColor', 'none', 'FaceAlpha', 0.88);
    xline(p99, '-', 'Color', C.orange, 'LineWidth', 1.5);
    xline(max(vals), '--', 'Color', C.red, 'LineWidth', 1.2);
    xlabel('Integration-point D_f [-]');
    ylabel('Count');
    title('');
    text(0.03, 0.94, '(a) Distribution', 'Units', 'normalized', ...
        'FontSize', 8.5, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'Margin', 1.5);
    ylabel('Count ($\times 10^4$)', 'Interpreter', 'latex');
    yticks([0 1e4 2e4]);
    yticklabels({'0','1','2'});
    legend({'$D_f$','P99','max'}, 'Interpreter', 'latex', ...
        'Box', 'off', 'Location', 'northwest', 'FontSize', 7.4);
    style_axes(gca);

    nexttile;
    scatter(xy(:,1)*1e3, xy(:,2)*1e3, 8, vals, 'filled');
    hold on
    hot = vals >= p99;
    scatter(xy(hot,1)*1e3, xy(hot,2)*1e3, 28, vals(hot), 'o', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.45);
    axis equal tight
    xlabel('x [mm]');
    ylabel('y [mm]');
    title('');
    text(0.03, 0.92, '(b) Full-field $D_f$', 'Units', 'normalized', ...
        'Interpreter', 'latex', 'FontSize', 7.8, 'FontWeight', 'bold', ...
        'BackgroundColor', 'w', 'Margin', 1.5);
    colormap(gca, hot_colormap());
    cb = colorbar('eastoutside');
    ylabel(cb, 'D_f [-]');
    clim([min(vals) max(vals)+eps]);
    style_axes(gca);

    nexttile;
    zoom = xy(:,1) <= 0.20e-3 & xy(:,2) <= 0.20e-3;
    scatter(xy(zoom,1)*1e3, xy(zoom,2)*1e3, 20, [0.78 0.80 0.82], 'filled');
    hold on
    hotz = hot & zoom;
    scatter(xy(hotz,1)*1e3, xy(hotz,2)*1e3, 70, vals(hotz), 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.4);
    axis equal
    xlim([0 0.20]); ylim([0 0.20]);
    xlabel('x [mm]');
    ylabel('y [mm]');
    title('');
    text(0.03, 0.91, '(c) Corner zoom', 'Units', 'normalized', ...
        'FontSize', 7.8, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'Margin', 1.5);
    colormap(gca, hot_colormap());
    clim([0.5 1.0]);
    text(0.105, 0.065, sprintf('P99 = %.3f\nf(D_f>0.9) = %.1e', ...
        p99, nnz(vals>0.9)/numel(vals)), 'FontSize', 6.5, ...
        'HorizontalAlignment', 'center', 'Interpreter', 'none', ...
        'BackgroundColor', 'w', 'Margin', 1.5);
    style_axes(gca);

    export_pub(fig, outdir, "publication_fatigue_localization");
    close(fig);
end

function cmap = hot_colormap()
    cmap = [linspace(1, 0.55, 128)', linspace(0.86, 0.05, 128)', linspace(0.35, 0.05, 128)'];
end

function [nx, ny] = parse_grid_labels(labels)
    nx = zeros(numel(labels), 1);
    ny = zeros(numel(labels), 1);
    for i = 1:numel(labels)
        tok = regexp(labels(i), 'G(\d+)x(\d+)', 'tokens', 'once');
        nx(i) = str2double(tok{1});
        ny(i) = str2double(tok{2});
    end
end

function [P, Temp, Freq] = parse_orthogonal_labels(labels)
    P = zeros(numel(labels), 1);
    Temp = zeros(numel(labels), 1);
    Freq = zeros(numel(labels), 1);
    for i = 1:numel(labels)
        tok = regexp(labels(i), 'P(\d+)_T(\d+)_f(\d+)', 'tokens', 'once');
        P(i) = str2double(tok{1});
        Temp(i) = str2double(tok{2});
        Freq(i) = str2double(tok{3}) / 1000;
    end
end
