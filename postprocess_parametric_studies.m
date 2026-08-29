function T = postprocess_parametric_studies(resultRoot, outdir, manifest)
    %POSTPROCESS_PARAMETRIC_STUDIES Summarize parametric study results.

    if nargin < 1 || isempty(resultRoot)
        resultRoot = './Results_Parametric';
    end
    if nargin < 2 || isempty(outdir)
        outdir = './Figures_Parametric';
    end
    if nargin < 3
        manifest = table();
    end
    if ~exist(outdir, 'dir'), mkdir(outdir); end

    addpath(genpath('./PostProcessing'))
    addpath(genpath('./Models'))
    addpath(genpath('./Shapes'))

    opts = pub_style();
    opts.outputDir = outdir;

    folders = {};
    labels = {};
    studies = {};
    if ~isempty(manifest)
        requiredVars = {'Case','Study'};
        if ~all(ismember(requiredVars, manifest.Properties.VariableNames))
            error('postprocess_parametric_studies:invalidManifest', ...
                'Manifest must contain Case and Study columns.');
        end
        for i = 1:height(manifest)
            folder = fullfile(resultRoot, manifest.Case(i));
            endFile = fullfile(folder, 'end.mat');
            if ~isfile(endFile)
                error('postprocess_parametric_studies:incompleteCampaign', ...
                    'Required final result is missing: %s', endFile);
            end
            folders{end+1} = char(folder); %#ok<AGROW>
            labels{end+1} = char(manifest.Case(i)); %#ok<AGROW>
            studies{end+1} = char(manifest.Study(i)); %#ok<AGROW>
        end
        dirs = dir(resultRoot);
        actual = string({dirs([dirs.isdir]).name});
        actual = actual(~startsWith(actual,'.'));
        extras = setdiff(actual(:), string(manifest.Case));
        if ~isempty(extras)
            error('postprocess_parametric_studies:extraCases', ...
                'Unexpected case directories in canonical campaign: %s', ...
                strjoin(extras, ', '));
        end
    else
        dirs = dir(resultRoot);
        for i = 1:numel(dirs)
            if ~dirs(i).isdir || startsWith(dirs(i).name, '.')
                continue
            end
            folder = fullfile(resultRoot, dirs(i).name);
            if isempty(dir(fullfile(folder, '*.mat')))
                continue
            end
            folders{end+1} = folder; %#ok<AGROW>
            labels{end+1} = dirs(i).name; %#ok<AGROW>
            studies{end+1} = classify_study(dirs(i).name); %#ok<AGROW>
        end
    end

    T = collect_scenario_metrics(folders, labels);
    if isempty(T)
        warning('No parametric results found in %s.', resultRoot);
        return
    end
    if height(T) ~= numel(labels)
        error('postprocess_parametric_studies:skippedCase', ...
            'Expected %d summary rows but collected %d.', numel(labels), height(T));
    end
    T.Study = studies(:);
    T = movevars(T, 'Study', 'Before', 'Label');
    if ~isempty(manifest)
        T = canonicalize_submission_summary(T);
    end

    % Write the canonical CSV exactly once.  All endpoint tables and figures
    % consume this same validated table rather than independently reloading MAT.
    summaryCsv = fullfile(outdir, 'parametric_summary.csv');
    writetable(T, summaryCsv);
    T = readtable(summaryCsv, 'TextType', 'string');
    generate_results_table(T, [], fullfile(outdir, 'parametric_summary.tex'), false);

    plot_grouped_metric(T, 'grid', 'Grid independence', outdir, opts);
    plot_grouped_metric(T, 'length_scale', 'Phase-field length-scale sensitivity', outdir, opts);
    plot_orthogonal_main_effects(T, outdir, opts);
end

function study = classify_study(name)
    if startsWith(name, 'Grid_')
        study = 'grid';
    elseif startsWith(name, 'LengthScale_')
        study = 'length_scale';
    elseif startsWith(name, 'Orthogonal_')
        study = 'orthogonal';
    else
        study = 'other';
    end
end

function plot_grouped_metric(T, studyName, plotTitle, outdir, opts)
    idx = strcmp(T.Study, studyName);
    if ~any(idx), return; end
    S = T(idx,:);

    fig = figure('Visible', 'on');
    set(fig, 'Position', opts.fig_double, 'Color', 'w');

    x = 1:height(S);
    yyaxis left
    require_columns(S, {'TimeScreening1Days','TimeScreening07Days'});
    t1 = S.TimeScreening1Days;
    t07 = S.TimeScreening07Days;
    plot(x, t1, '-o', 'Color', opts.colors(1,:), ...
        'LineWidth', opts.linewidth, 'MarkerFaceColor', opts.colors(1,:));
    ylabel('$t_{I_S=1}$ [days]', 'Interpreter', 'latex', ...
        'FontSize', opts.fontsize_label);

    yyaxis right
    plot(x, t07, '-s', 'Color', opts.colors(7,:), ...
        'LineWidth', opts.linewidth, 'MarkerFaceColor', opts.colors(7,:));
    ylabel('$t_{I_S=0.7}$ [days]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);

    set(gca, 'XTick', x, 'XTickLabel', S.Label, opts.axes_props{:});
    xtickangle(30);
    title(plotTitle, 'FontName', opts.fontname, 'FontSize', opts.fontsize_title);
    export_fig(fig, ['parametric_' studyName], opts, 'both');
end

function plot_orthogonal_main_effects(T, outdir, opts)
    O = T(strcmp(T.Study, 'orthogonal'),:);
    if isempty(O), return; end

    P = zeros(height(O),1);
    Temp = zeros(height(O),1);
    Freq = zeros(height(O),1);
    for i = 1:height(O)
        tok = regexp(O.Label{i}, 'P(\d+)_T(\d+)_f(\d+)', 'tokens', 'once');
        P(i) = str2double(tok{1});
        Temp(i) = str2double(tok{2});
        Freq(i) = str2double(tok{3}) / 1000;
    end

    fig = figure('Visible', 'on');
    set(fig, 'Position', opts.fig_double, 'Color', 'w');
    tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

    require_columns(O, {'TimeScreening1Days'});
    t1 = O.TimeScreening1Days;
    nexttile; plot_effect(P, t1, opts, 'Pressure [MPa]');
    nexttile; plot_effect(Temp, t1, opts, 'Temperature [$^\circ$C]');
    nexttile; plot_effect(Freq, t1, opts, 'Frequency [Hz]');

    export_fig(fig, 'parametric_orthogonal_main_effects', opts, 'both');
end

function plot_effect(x, y, opts, xlab)
    levels = unique(x);
    meanY = zeros(size(levels));
    for i = 1:numel(levels)
        meanY(i) = mean(y(x == levels(i)));
    end
    plot(levels, meanY, '-o', 'Color', opts.colors(1,:), ...
        'LineWidth', opts.linewidth, 'MarkerFaceColor', opts.colors(1,:));
    xlabel(xlab, 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    ylabel('Mean $t_{I_S=1}$ [days]', 'Interpreter', 'latex', 'FontSize', opts.fontsize_label);
    set(gca, opts.axes_props{:});
end

function require_columns(T, names)
missing = setdiff(names,T.Properties.VariableNames,'stable');
if ~isempty(missing)
    error('postprocess_parametric_studies:missingColumn', ...
        'Corrected table lacks: %s.',strjoin(missing,', '));
end
end
