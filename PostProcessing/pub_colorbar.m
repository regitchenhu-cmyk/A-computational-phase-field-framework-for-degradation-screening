function cb = pub_colorbar(label_text, opts, varargin)
    %PUB_COLORBAR Create publication-quality colorbar.
    %
    % Usage:
    %   cb = pub_colorbar('\phi [-]', opts);
    %   cb = pub_colorbar('$\sigma_H$ [MPa]', opts, 'Location', 'southoutside');
    
    cb = colorbar(varargin{:});
    cb.FontName = opts.fontname;
    cb.FontSize = opts.fontsize_tick;
    cb.LineWidth = 0.5;
    cb.TickDirection = 'out';
    
    if ~isempty(label_text)
        cb.Label.String = label_text;
        cb.Label.FontSize = opts.fontsize_label;
        cb.Label.FontName = opts.fontname;
        cb.Label.Interpreter = 'latex';
    end
end
