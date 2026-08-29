function T = diagnose_fatigue_localization(resultRoot, outputDir)
    %DIAGNOSE_FATIGUE_LOCALIZATION Quantify integration-point fatigue spikes.
    %
    % The fatigue model stores D_f at integration points. A plain maximum can
    % be dominated by one or two pathological points, so this report compares
    % max values against percentiles and damaged integration-point fractions.

    if nargin < 1 || isempty(resultRoot)
        resultRoot = './Results_Parametric_Calibrated';
    end
    if nargin < 2 || isempty(outputDir)
        outputDir = './Figures_Parametric_Calibrated';
    end
    if ~exist(outputDir, 'dir'), mkdir(outputDir); end

    addpath(genpath('./Models'))
    addpath(genpath('./Shapes'))

    folders = dir(resultRoot);
    folders = folders([folders.isdir]);
    folders = folders(~ismember({folders.name}, {'.','..'}));

    rows = {};
    for i = 1:numel(folders)
        caseName = folders(i).name;
        matfile = fullfile(resultRoot, caseName, 'end.mat');
        if ~exist(matfile, 'file')
            continue
        end

        data = load(matfile, 'physics', 'mesh', 'results', 'tvec');
        [Df, modelIndex] = get_fatigue_ip(data.physics);
        if isempty(Df)
            warning('No fatigue integration-point field found in %s', caseName);
            continue
        end

        flat = Df(:);
        nTotal = numel(flat);
        p = prctile(flat, [50 90 95 99 99.9]);
        rows(end+1,:) = { ...
            caseName, matfile, modelIndex, nTotal, ...
            max(flat), mean(flat), p(1), p(2), p(3), p(4), p(5), ...
            nnz(flat > 0.5), nnz(flat > 0.7), nnz(flat > 0.9), ...
            nnz(flat > 0.5)/nTotal, nnz(flat > 0.7)/nTotal, nnz(flat > 0.9)/nTotal}; %#ok<AGROW>

        if contains(caseName, 'Grid_')
            plot_fatigue_ip_map(data, Df, caseName, outputDir);
        end
    end

    T = cell2table(rows, 'VariableNames', { ...
        'Case','File','FatigueModelIndex','NumIntegrationPoints', ...
        'DfMax','DfMean','DfP50','DfP90','DfP95','DfP99','DfP999', ...
        'CountDfGt05','CountDfGt07','CountDfGt09', ...
        'FracDfGt05','FracDfGt07','FracDfGt09'});

    csvFile = fullfile(outputDir, 'fatigue_localization_diagnostics.csv');
    writetable(T, csvFile);
    fprintf('Wrote fatigue localization diagnostics: %s\n', csvFile);
end

function [Df, modelIndex] = get_fatigue_ip(physics)
    Df = [];
    modelIndex = NaN;
    for m = 1:numel(physics.models)
        mdl = physics.models{m};
        if isprop(mdl, 'Df_ip') && ~isempty(mdl.Df_ip)
            Df = mdl.Df_ip;
            modelIndex = m;
            return
        end
    end
end

function plot_fatigue_ip_map(data, Df, caseName, outputDir)
    mesh = data.mesh;
    g = mesh.getGroupIndex("Internal");
    nElem = size(Df, 1);
    nIp = size(Df, 2);
    xy = zeros(nElem*nIp, 2);
    vals = Df(:);

    row = 0;
    for e = 1:nElem
        eip = mesh.getIPCoords(g, e);
        if size(eip, 2) ~= 2 && size(eip, 1) == 2
            eip = eip.';
        end
        for ip = 1:nIp
            row = row + 1;
            xy(row,:) = eip(ip,:);
        end
    end

    fig = figure('Visible', 'off', 'Color', 'w');
    scatter(xy(:,1)*1e3, xy(:,2)*1e3, 12, vals, 'filled');
    axis equal tight
    box on
    xlabel('x [mm]');
    ylabel('y [mm]');
    title(sprintf('%s: integration-point D_f', strrep(caseName, '_', '\_')));
    cb = colorbar;
    ylabel(cb, 'D_f [-]');
    clim([0 1]);
    colormap(turbo);

    pngFile = fullfile(outputDir, sprintf('fatigue_ip_map_%s.png', caseName));
    exportgraphics(fig, pngFile, 'Resolution', 300);
    close(fig);
end
