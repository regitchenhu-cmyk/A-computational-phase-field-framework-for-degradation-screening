function opts = pub_style()
    %PUB_STYLE Returns a struct of publication-quality style settings
    % for aerospace hydraulic seal simulation figures.
    %
    % Usage:
    %   opts = pub_style();
    %   set(gcf, 'Position', opts.fig_single);
    %   set(gca, opts.axes_props{:});
    %
    % Targets: AIAA Journal / Tribology International / IJMS style
    
    %% Figure sizes (pixels at 300 DPI)
    % Single column: 90mm = 3.54in
    % Double column: 190mm = 7.48in
    opts.fig_single = [100 100 540 420];    % ~3.54" x 2.76"
    opts.fig_double = [100 100 1120 420];   % ~7.48" x 2.76"
    opts.fig_tall   = [100 100 540 680];    % ~3.54" x 4.5"
    opts.fig_wide   = [100 100 1120 840];   % ~7.48" x 5.5"
    opts.fig_square = [100 100 540 540];    % ~3.54" x 3.54"
    
    %% Font settings
    opts.fontname   = 'Times New Roman';
    opts.fontsize   = 11;          % Body text
    opts.fontsize_label = 12;      % Axis labels
    opts.fontsize_title = 13;      % Subplot titles
    opts.fontsize_legend = 10;     % Legend
    opts.fontsize_tick  = 10;      % Tick labels
    opts.fontsize_annotation = 9;  % Annotations / subfigure labels
    
    %% Line styles
    opts.linewidth = 1.5;
    opts.linewidth_thin = 0.8;
    opts.linewidth_thick = 2.0;
    opts.markersize = 6;
    opts.markersize_small = 4;
    
    %% Color palettes
    % Primary palette (colorblind-safe, print-friendly)
    opts.colors = [
        0.000 0.447 0.741;   % Blue
        0.850 0.325 0.098;   % Red-Orange
        0.466 0.674 0.188;   % Green
        0.494 0.184 0.556;   % Purple
        0.929 0.694 0.125;   % Gold
        0.301 0.745 0.933;   % Light Blue
        0.635 0.078 0.184;   % Dark Red
        0.100 0.100 0.100;   % Near Black
    ];
    
    % Sequential colormaps for contour plots
    opts.cmap_damage  = hot(256);           % Phase-field damage
    opts.cmap_stress  = jet(256);           % Stress fields
    opts.cmap_fluid   = parula(256);        % Fluid concentration
    opts.cmap_aging   = copper(256);        % Aging
    opts.cmap_fatigue = cool(256);          % Fatigue damage
    
    % Custom sequential colormap (white -> blue -> dark blue)
    r = linspace(1, 0.05, 256)';
    g = linspace(1, 0.20, 256)';
    b = linspace(1, 0.55, 256)';
    opts.cmap_seq_blue = [r g b];
    
    % Custom diverging (blue-white-red)
    n = 128;
    r1 = [linspace(0.2, 1, n)'; linspace(1, 0.7, n)'];
    g1 = [linspace(0.3, 1, n)'; linspace(1, 0.15, n)'];
    b1 = [linspace(0.7, 1, n)'; linspace(1, 0.15, n)'];
    opts.cmap_diverge = [r1 g1 b1];
    
    %% Standard axes properties (apply via set(gca, ...))
    opts.axes_props = { ...
        'FontName',     opts.fontname, ...
        'FontSize',     opts.fontsize_tick, ...
        'LineWidth',    0.6, ...
        'Box',          'on', ...
        'TickDir',      'in', ...
        'TickLength',   [0.015 0.015], ...
        'XMinorTick',   'on', ...
        'YMinorTick',   'on', ...
        'XGrid',        'off', ...
        'YGrid',        'off', ...
        'Layer',        'top' ...
    };
    
    %% Export settings
    opts.dpi = 600;
    opts.format = '-dpng';         % '-dpng', '-depsc2', '-dtiff', '-dpdf'
    opts.format_vec = '-depsc2';   % Vector format for line plots
    opts.format_ras = '-dpng';     % Raster format for contour plots
    
    %% Unit conversion helpers
    opts.Pa2MPa = 1e-6;
    opts.m2mm   = 1e3;
    opts.s2hr   = 1/3600;
    opts.s2day  = 1/86400;
    opts.s2yr   = 1/(86400*365);
end
