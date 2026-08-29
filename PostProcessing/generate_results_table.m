function T = generate_results_table(folders, labels, outfile, writeCsv)
    %GENERATE_RESULTS_TABLE Export final scenario metrics as LaTeX and CSV.

    addpath(genpath('./Models'))
    addpath(genpath('./Shapes'))

    if nargin < 3 || isempty(outfile)
        outfile = './Figures/table_results.tex';
    end
    if nargin < 4 || isempty(writeCsv)
        writeCsv = true;
    end

    outdir = fileparts(outfile);
    if ~isempty(outdir) && ~exist(outdir, 'dir')
        mkdir(outdir);
    end

    if istable(folders)
        % Canonical pipelines pass the already validated summary table so the
        % LaTeX table and figures cannot silently reload different MAT files.
        T = folders;
    else
        if nargin < 2 || isempty(labels)
            labels = folders;
        end
        T = collect_scenario_metrics(folders, labels);
    end

    fid = fopen(outfile, 'w');
    if fid < 0
        error('Cannot open table output file: %s', outfile);
    end
    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, '\\begin{table}[htbp]\n');
    fprintf(fid, '\\centering\n');
    fprintf(fid, '\\caption{Summary of seal degradation simulation results.}\n');
    fprintf(fid, '\\label{tab:seal_results}\n');
    fprintf(fid, '\\begin{tabular}{lcccccccc}\n');
    fprintf(fid, '\\hline\\hline\n');
    fprintf(fid, 'Scenario & $t_{\\max}$ & $t_{I_S=1}$ & $L_{\\mathrm{deg}}$ & $\\bar{C}_L$ & $\\alpha_{\\max}$ & $D_{f,\\max}$ & Cycles & $I_S$ \\\\\n');
    fprintf(fid, ' & [days] & [days] & [mm] & [-] & [-] & [-] & [-] & [-] \\\\\n');
    fprintf(fid, '\\hline\n');

    require_columns(T, {'TimeScreening1Days','ConnectedDegradationExtentMm','ScreeningIndex'});
    timeScreening1 = T.TimeScreening1Days;
    degradationExtent = T.ConnectedDegradationExtentMm;
    screeningIndex = T.ScreeningIndex;
    for s = 1:height(T)
        fprintf(fid, '%s & %.1f & %.3f & %.3f & %.4f & %.4f & %.4f & %.0f & %.4f \\\\\n', ...
            char(string(T.Label(s))), T.TimeDays(s), timeScreening1(s), degradationExtent(s), T.FluidAvg(s), ...
            T.AgingMax(s), T.FatigueMax(s), T.Cycles(s), screeningIndex(s));
    end

    fprintf(fid, '\\hline\\hline\n');
    fprintf(fid, '\\end{tabular}\n');
    fprintf(fid, '\\end{table}\n');

    if writeCsv
        csvfile = replace(outfile, '.tex', '.csv');
        writetable(T, csvfile);
        fprintf('CSV summary written: %s\n', csvfile);
    end

    fprintf('LaTeX table written: %s\n', outfile);
    disp(T);
end

function require_columns(T,names)
missing=setdiff(names,T.Properties.VariableNames,'stable');
if ~isempty(missing)
    error('generate_results_table:missingColumn', ...
        'Corrected table lacks: %s.',strjoin(missing,', '));
end
end
