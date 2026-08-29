function T = postprocess_corrected_submission(resultRoot, outdir)
    %POSTPROCESS_CORRECTED_SUBMISSION Build corrected tables and figures.

    if nargin < 1 || isempty(resultRoot)
        resultRoot = './Results_Parametric_Corrected';
    end
    if nargin < 2 || isempty(outdir)
        outdir = './Figures_Parametric_Corrected';
    end
    if ~exist(outdir, 'dir'), mkdir(outdir); end

    manifest = scientifically_corrected_manifest_20260823();
    T = postprocess_parametric_studies(resultRoot, outdir, manifest);
    summaryCsv = fullfile(outdir, 'parametric_summary.csv');
    % Re-read the one authoritative CSV before every downstream endpoint
    % artifact.  This also catches serialization/order mistakes immediately.
    C = readtable(summaryCsv, 'TextType', 'string');
    if height(C) ~= height(manifest) || ~isequal(string(C.Label), manifest.Case)
        error('postprocess_corrected_submission:canonicalOrder', ...
            'Canonical summary CSV does not follow the 16-case manifest.');
    end
    T = C;
    plot_cmame_style_figures(resultRoot, outdir, summaryCsv);
    plot_publication_parametric_v2(summaryCsv, resultRoot, outdir);
    write_l9_range_table(T, outdir);
end

function write_l9_range_table(T, outdir)
    O = T(strcmp(T.Study, 'orthogonal'), :);
    P = zeros(height(O),1); Temp = P; Freq = P;
    for i = 1:height(O)
        tok = regexp(char(string(O.Label(i))), ...
            'P(\d+)_T(\d+)_f(\d+)', 'tokens', 'once');
        P(i) = str2double(tok{1});
        Temp(i) = str2double(tok{2});
        Freq(i) = str2double(tok{3})/1000;
    end

    responseNames = ["ScreeningIndex"; "FatigueP99"; "AgingP99"];
    responses = {O.ScreeningIndex, O.FatigueP99, O.AgingP99};
    rows = cell(numel(responses), 14);
    for r = 1:numel(responses)
        y = responses{r};
        pMeans = level_means(P, y);
        tMeans = level_means(Temp, y);
        fMeans = level_means(Freq, y);
        ranges = [range(pMeans), range(tMeans), range(fMeans)];
        labels = ["Pressure", "Temperature", "Frequency"];
        ranking = rank_effects(labels, ranges, y);
        rows(r,:) = {responseNames(r), pMeans(1), pMeans(2), pMeans(3), ...
            tMeans(1), tMeans(2), tMeans(3), fMeans(1), fMeans(2), fMeans(3), ...
            ranges(1), ranges(2), ranges(3), ranking};
    end
    R = cell2table(rows, 'VariableNames', {'Response','PressureLevel1Mean', ...
        'PressureLevel2Mean','PressureLevel3Mean','TemperatureLevel1Mean', ...
        'TemperatureLevel2Mean','TemperatureLevel3Mean','FrequencyLevel1Mean', ...
        'FrequencyLevel2Mean','FrequencyLevel3Mean','PressureRange', ...
        'TemperatureRange','FrequencyRange','RangeRanking'});
    writetable(R, fullfile(outdir, 'l9_main_effect_ranges.csv'));
end

function ranking = rank_effects(labels, ranges, response)
    % Do not rank effects that differ only by floating-point round-off.
    tol = max(1e-12, 1e-9*max(abs(response), [], 'omitnan'));
    clean = ranges;
    clean(abs(clean) <= tol) = 0;
    [sorted, order] = sort(clean, 'descend');
    ranking = labels(order(1));
    for i = 2:numel(order)
        if abs(sorted(i) - sorted(i-1)) <= tol
            separator = " = ";
        else
            separator = " > ";
        end
        ranking = ranking + separator + labels(order(i));
    end
end

function mu = level_means(x, y)
    levels = unique(x);
    mu = zeros(numel(levels), 1);
    for i = 1:numel(levels)
        mu(i) = mean(y(x == levels(i)), 'omitnan');
    end
end
