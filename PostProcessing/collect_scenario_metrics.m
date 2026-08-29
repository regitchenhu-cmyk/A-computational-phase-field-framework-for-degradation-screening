function T = collect_scenario_metrics(folders, labels)
    %COLLECT_SCENARIO_METRICS Build a table of final metrics for scenarios.

    if nargin < 2 || isempty(labels)
        labels = folders;
    end

    rows = {};
    for s = 1:numel(folders)
        try
            [data, filepath] = load_latest_result(folders{s});
        catch ME
            warning('Skipping %s: %s', folders{s}, ME.message);
            continue
        end

        res = data.results;
        [screeningIndex, comp] = compute_screening_index(data);
        fatigueStats = fatigue_ip_stats(data);
        tdays = data.tvec(:) / 86400;
        tScreening07 = first_crossing_time(tdays, screeningIndex, 0.7);
        tScreening1 = first_crossing_time(tdays, screeningIndex, 1.0);
        tPhaseOne = first_crossing_time(tdays, comp.phase, 1.0);
        if isfield(res, 'degradation_connected')
            spanningPath = logical(res.degradation_connected(end));
        else
            spanningPath = false;
        end
        rows(end+1,:) = { ...
            char(labels{s}), char(folders{s}), char(filepath), ...
            data.tvec(end) / 86400, ...
            res.LFrac(end) * 1e3, spanningPath, res.CL_avg(end), res.CL_max(end), ...
            res.alpha_avg(end), res.alpha_max(end), res.Df_max(end), ...
            res.totalCycles(end), screeningIndex(end), ...
            comp.phase(end), comp.fatigue(end), comp.aging(end), ...
            tScreening07, tScreening1, tPhaseOne, ...
            fatigueStats.mean, fatigueStats.p95, fatigueStats.p999, ...
            fatigueStats.fracGt05, fatigueStats.fracGt07, fatigueStats.fracGt09}; %#ok<AGROW>
    end

    if isempty(rows)
        T = table();
        return
    end

    T = cell2table(rows, 'VariableNames', { ...
        'Label','Folder','File','TimeDays','ConnectedDegradationExtentMm','SpanningDegradationPath','FluidAvg', ...
        'FluidMax','AgingAvg','AgingMax','FatigueMax','Cycles', ...
        'ScreeningIndex','PhaseP99','FatigueP99','AgingP99', ...
        'TimeScreening07Days','TimeScreening1Days','TimePhaseP99OneDays', ...
        'FatigueMeanIP','FatigueP95IP','FatigueP999IP', ...
        'FatigueFracGt05IP','FatigueFracGt07IP','FatigueFracGt09IP'});
end

function t = first_crossing_time(tvec, y, threshold)
    idx = find(y(:) >= threshold, 1, 'first');
    if isempty(idx)
        t = NaN;
    else
        t = tvec(idx);
    end
end

function stats = fatigue_ip_stats(data)
    stats = struct('mean', NaN, 'p95', NaN, 'p99', NaN, 'p999', NaN, ...
        'fracGt05', NaN, 'fracGt07', NaN, 'fracGt09', NaN);
    if ~isfield(data, 'physics') || isempty(data.physics)
        return
    end
    Df = [];
    for m = 1:numel(data.physics.models)
        mdl = data.physics.models{m};
        if isprop(mdl, 'Df_ip') && ~isempty(mdl.Df_ip)
            Df = mdl.Df_ip(:);
            break
        end
    end
    if isempty(Df)
        return
    end
    stats.mean = mean(Df);
    p = prctile(Df, [95 99 99.9]);
    stats.p95 = p(1);
    stats.p99 = p(2);
    stats.p999 = p(3);
    stats.fracGt05 = nnz(Df > 0.5) / numel(Df);
    stats.fracGt07 = nnz(Df > 0.7) / numel(Df);
    stats.fracGt09 = nnz(Df > 0.9) / numel(Df);
end
