function plot_crack_growth_rate(data, opts, outdir)
    %PLOT_CRACK_GROWTH_RATE Plot crack growth rate da/dN versus crack length.

    if nargin < 2 || isempty(opts)
        opts = pub_style();
    end
    if nargin < 3 || isempty(outdir)
        outdir = './Figures';
    end
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end

    res = data.results;
    if ~isfield(res, 'totalCycles') || numel(res.totalCycles) < 3
        warning('Crack growth rate requires totalCycles history.');
        return
    end

    N = res.totalCycles(:);
    a = res.LFrac(:) * opts.m2mm;
    valid = isfinite(N) & isfinite(a);
    N = N(valid);
    a = a(valid);
    [N, idx] = unique(N, 'stable');
    a = a(idx);

    if numel(N) < 3
        warning('Not enough unique cycle points for crack growth rate.');
        return
    end

    dadN = gradient(a) ./ max(gradient(N), eps);
    dadN = max(dadN, eps);

    fig = figure('Visible', 'on');
    set(fig, 'Position', opts.fig_single, 'Color', 'w');

    semilogy(a, dadN, '-', 'Color', opts.colors(1,:), ...
        'LineWidth', opts.linewidth);
    xlabel('Crack length $a$ [mm]', 'Interpreter', 'latex', ...
        'FontSize', opts.fontsize_label);
    ylabel('$da/dN$ [mm/cycle]', 'Interpreter', 'latex', ...
        'FontSize', opts.fontsize_label);
    set(gca, opts.axes_props{:});
    xlim([0 max(0.1, max(a) * 1.05)]);

    export_fig(fig, 'fig15_crack_growth_rate', opts, 'both');
end
