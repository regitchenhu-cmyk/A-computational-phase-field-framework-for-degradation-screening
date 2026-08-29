function add_subfig_label(ax, label_char, opts, position)
    %ADD_SUBFIG_LABEL Add (a), (b), (c) style label to subplot.
    %
    % Usage:
    %   add_subfig_label(gca, 'a', opts)
    %   add_subfig_label(gca, 'b', opts, 'northwest')  % default: northwest
    
    if nargin < 4
        position = 'northwest';
    end
    
    axes(ax);
    xl = xlim; yl = ylim;
    dx = xl(2) - xl(1);
    dy = yl(2) - yl(1);
    
    switch position
        case 'northwest'
            xp = xl(1) + 0.02*dx;
            yp = yl(2) - 0.06*dy;
            ha = 'left'; va = 'top';
        case 'northeast'
            xp = xl(2) - 0.02*dx;
            yp = yl(2) - 0.06*dy;
            ha = 'right'; va = 'top';
        case 'southwest'
            xp = xl(1) + 0.02*dx;
            yp = yl(1) + 0.06*dy;
            ha = 'left'; va = 'bottom';
        case 'southeast'
            xp = xl(2) - 0.02*dx;
            yp = yl(1) + 0.06*dy;
            ha = 'right'; va = 'bottom';
    end
    
    text(xp, yp, ['(' label_char ')'], ...
        'FontName', opts.fontname, ...
        'FontSize', opts.fontsize_title, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', ha, ...
        'VerticalAlignment', va, ...
        'BackgroundColor', 'w', ...
        'EdgeColor', 'none', ...
        'Margin', 2);
end
