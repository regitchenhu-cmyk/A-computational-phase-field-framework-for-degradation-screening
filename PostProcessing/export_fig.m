function export_fig(fig, filename, opts, format_type)
    %EXPORT_FIG Export figure with publication-quality settings.
    %
    % Usage:
    %   export_fig(gcf, 'fig01_damage', opts)           % Default PNG 600dpi
    %   export_fig(gcf, 'fig01_damage', opts, 'vector')  % EPS vector
    %   export_fig(gcf, 'fig01_damage', opts, 'both')    % PNG + EPS
    %   export_fig(gcf, 'fig01_damage', opts, 'pdf')     % PDF
    %
    % Creates output in opts.outputDir when provided, otherwise ./Figures/.
    
    if nargin < 4
        format_type = 'raster';
    end
    
    if isfield(opts, 'outputDir') && ~isempty(opts.outputDir)
        outdir = opts.outputDir;
    else
        outdir = './Figures';
    end
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end
    
    % Apply final figure-level styling
    set(fig, 'Color', 'w');
    set(fig, 'PaperPositionMode', 'auto');
    set(fig, 'InvertHardcopy', 'off');
    
    % Tighten layout
    try
        set(fig, 'Units', 'pixels');
    catch
    end
    
    filepath = fullfile(outdir, filename);
    
    switch format_type
        case 'raster'
            print(fig, filepath, opts.format_ras, ['-r' num2str(opts.dpi)]);
            fprintf('Exported: %s.png (%d dpi)\n', filepath, opts.dpi);
            
        case 'vector'
            print(fig, filepath, opts.format_vec);
            fprintf('Exported: %s.eps\n', filepath);
            
        case 'pdf'
            print(fig, filepath, '-dpdf', '-bestfit');
            fprintf('Exported: %s.pdf\n', filepath);
            
        case 'both'
            print(fig, filepath, opts.format_ras, ['-r' num2str(opts.dpi)]);
            print(fig, filepath, opts.format_vec);
            fprintf('Exported: %s.png + .eps\n', filepath);
            
        case 'tiff'
            print(fig, filepath, '-dtiff', ['-r' num2str(opts.dpi)]);
            fprintf('Exported: %s.tiff (%d dpi)\n', filepath, opts.dpi);
    end
end
