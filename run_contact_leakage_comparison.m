function comparison = run_contact_leakage_comparison(scenario, overrides)
    %RUN_CONTACT_LEAKAGE_COMPARISON Compare fixed-BC and contact-aware seal runs.
    %
    % Example:
    %   comparison = run_contact_leakage_comparison(1, struct('tmaxDays', 2, 'nMax', 80));
    
    if nargin < 1 || isempty(scenario)
        scenario = 1;
    end
    if nargin < 2 || isempty(overrides)
        overrides = struct();
    end
    
    base = overrides;
    base.plotEvery = get_option(base, 'plotEvery', inf);
    base.saveEvery = get_option(base, 'saveEvery', inf);
    base.verboseEvery = get_option(base, 'verboseEvery', inf);
    base.useParpool = get_option(base, 'useParpool', false);
    base.cleanResults = get_option(base, 'cleanResults', true);
    base.saveRoot = get_option(base, 'saveRoot', './Results_ContactLeakageComparison');
    
    fixed = base;
    fixed.useContactPenalty = false;
    fixed.nameSuffix = "fixedBC";
    
    contact = base;
    contact.useContactPenalty = true;
    contact.nameSuffix = "contactAware";
    
    fprintf('\n=== Running fixed displacement boundary case ===\n');
    [physicsFixed, tFixed, resultsFixed] = main_seal(scenario, fixed);
    
    fprintf('\n=== Running contact-aware penalty boundary case ===\n');
    [physicsContact, tContact, resultsContact] = main_seal(scenario, contact);
    
    comparison.fixed.physics = physicsFixed;
    comparison.fixed.tvec = tFixed;
    comparison.fixed.results = resultsFixed;
    comparison.contact.physics = physicsContact;
    comparison.contact.tvec = tContact;
    comparison.contact.results = resultsContact;
    
    comparison.summary = summarize_pair(tFixed, resultsFixed, tContact, resultsContact);
    
    fig = plot_contact_leakage_comparison(tFixed, resultsFixed, tContact, resultsContact);
    outDir = fullfile(string(base.saveRoot), "figures");
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end
    outFile = fullfile(outDir, "contact_leakage_comparison.png");
    exportgraphics(fig, outFile, 'Resolution', 300);
    comparison.figure = outFile;
    
    disp(comparison.summary)
    fprintf('Saved comparison figure: %s\n', outFile);
end

function summary = summarize_pair(tFixed, rFixed, tContact, rContact)
    summary = table;
    summary.Case = ["Fixed BC"; "Contact-aware"];
    summary.FinalTime_days = [tFixed(end); tContact(end)] / 86400;
    summary.FinalDfMax = [rFixed.Df_max(end); rContact.Df_max(end)];
    summary.FinalLFrac_mm = [rFixed.LFrac(end); rContact.LFrac(end)] * 1e3;
    summary.FinalLeakage_Q_m3s = [safe_last(rFixed, 'leakage_Q'); safe_last(rContact, 'leakage_Q')];
    summary.FinalPoiseuille_Q_m3s = [safe_last(rFixed, 'leakage_Q_poiseuille'); ...
                                     safe_last(rContact, 'leakage_Q_poiseuille')];
    summary.DarcyLeakTime_days = [safe_last(rFixed, 'leakage_time_darcy'); ...
                                  safe_last(rContact, 'leakage_time_darcy')];
    summary.PoiseuilleLeakTime_days = [safe_last(rFixed, 'leakage_time_poiseuille'); ...
                                       safe_last(rContact, 'leakage_time_poiseuille')];
    summary.LeakagePath = [logical(safe_last(rFixed, 'leakage_connected')); ...
                           logical(safe_last(rContact, 'leakage_connected'))];
    summary.FinalTopPc_MPa = [safe_last(rFixed, 'contact_top_mean'); ...
                              safe_last(rContact, 'contact_top_mean')] / 1e6;
end

function fig = plot_contact_leakage_comparison(tFixed, rFixed, tContact, rContact)
    fig = figure('Color', 'w', 'Position', [80 80 1100 760]);
    tf = tFixed / 86400;
    tc = tContact / 86400;
    
    subplot(2,2,1)
    plot(tf, rFixed.Df_max, 'k--', 'LineWidth', 1.4); hold on
    plot(tc, rContact.Df_max, 'r-', 'LineWidth', 1.5);
    xlabel('Time [days]'); ylabel('D_f max [-]');
    title('Fatigue Damage'); grid on
    legend('Fixed BC','Contact-aware','Location','best')
    
    subplot(2,2,2)
    plot(tf, rFixed.LFrac*1e3, 'k--', 'LineWidth', 1.4); hold on
    plot(tc, rContact.LFrac*1e3, 'r-', 'LineWidth', 1.5);
    xlabel('Time [days]'); ylabel('L_{frac} [mm]');
    title('Damage Front'); grid on
    
    subplot(2,2,3)
    semilogy(tf, max(safe_series(rFixed, 'leakage_Q'), realmin), 'k--', 'LineWidth', 1.4); hold on
    semilogy(tc, max(safe_series(rContact, 'leakage_Q'), realmin), 'r-', 'LineWidth', 1.5);
    xlabel('Time [days]'); ylabel('Q [m^3/s]');
    title('Leakage Rate'); grid on
    
    subplot(2,2,4)
    plot(tc, safe_series(rContact, 'contact_top_mean')/1e6, 'Color', [0.85 0.25 0.15], 'LineWidth', 1.5); hold on
    plot(tc, safe_series(rContact, 'contact_bottom_mean')/1e6, 'Color', [0.10 0.45 0.75], 'LineWidth', 1.3);
    plot(tc, safe_series(rContact, 'contact_right_mean')/1e6, 'Color', [0.20 0.60 0.30], 'LineWidth', 1.3);
    xlabel('Time [days]'); ylabel('Mean p_c [MPa]');
    title('Contact-aware Case'); grid on
    legend('Top','Bottom','Right','Location','best')
    
    sgtitle('Seal Boundary Model Comparison')
end

function y = safe_series(results, fieldName)
    if isfield(results, fieldName)
        y = results.(fieldName);
    else
        firstField = fieldnames(results);
        y = zeros(size(results.(firstField{1})));
    end
end

function y = safe_last(results, fieldName)
    s = safe_series(results, fieldName);
    y = s(end);
end

function value = get_option(s, name, defaultValue)
    if isfield(s, name) && ~isempty(s.(name))
        value = s.(name);
    else
        value = defaultValue;
    end
end
