function postprocess_seal(result_folder, steps_to_plot)
    %POSTPROCESS_SEAL Master post-processing for seal failure simulation.
    %
    % Generates publication-quality figures from saved simulation data.
    %
    % Usage:
    %   postprocess_seal('./Results_Seal/Standard_21MPa_80C/')
    %   postprocess_seal('./Results_Seal/Standard_21MPa_80C/', [50 100 200])
    %
    % Generates:
    %   Fig 1: Field contour maps (phi, sigma_H, CL, alpha, D_f)
    %   Fig 2: Time history curves (crack length, fluid, aging, fatigue)
    %   Fig 3: Degradation evolution snapshots
    %   Fig 4: Crack path and deformed configuration
    %   Fig 5: Reliability assessment dashboard
    %   Fig 6: Energy and convergence analysis
    
    addpath(genpath('./Models'))
    addpath(genpath('./Shapes'))
    addpath(genpath('./PostProcessing'))
    
    opts = pub_style();
    
    outdir = './Figures';
    if ~exist(outdir, 'dir'); mkdir(outdir); end
    
    %% Detect available result files
    files = dir(fullfile(result_folder, '*.mat'));
    
    if isempty(files)
        error('No .mat files found in: %s', result_folder);
    end
    
    % Separate numeric-named files (time steps) from named files (e.g. "end")
    available_steps = [];
    named_files = {};
    for i = 1:length(files)
        [~, name, ~] = fileparts(files(i).name);
        num = str2double(name);
        if ~isnan(num)
            available_steps = [available_steps, num];
        else
            named_files{end+1} = name;
        end
    end
    available_steps = sort(available_steps);
    
    % Diagnostic output
    fprintf('Directory: %s\n', result_folder);
    fprintf('  Found %d numeric step files', length(available_steps));
    if ~isempty(available_steps)
        fprintf(': [%d ... %d]', available_steps(1), available_steps(end));
    end
    fprintf('\n');
    if ~isempty(named_files)
        fprintf('  Found named files: %s\n', strjoin(named_files, ', '));
    end
    
    % Determine which file to load as the "final state"
    % Priority: "end.mat" > last numeric step > any other .mat
    has_end = any(strcmp(named_files, 'end'));
    
    if has_end
        fprintf('Loading final state from "end.mat"...\n');
        data_end = load(fullfile(result_folder, 'end.mat'));
    elseif ~isempty(available_steps)
        last_step = available_steps(end);
        fprintf('Loading final state (step %d)...\n', last_step);
        data_end = load(fullfile(result_folder, [num2str(last_step) '.mat']));
    else
        % Fall back: load the first .mat file found
        fprintf('No numeric steps or end.mat found. Loading: %s\n', files(1).name);
        data_end = load(fullfile(result_folder, files(1).name));
    end
    
    % Validate loaded data has required fields
    required_fields = {'physics', 'mesh', 'tvec', 'results'};
    for f = 1:length(required_fields)
        if ~isfield(data_end, required_fields{f})
            error('Loaded file is missing required field "%s". Check file contents.', ...
                required_fields{f});
        end
    end
    
    if nargin < 2 || isempty(steps_to_plot)
        if ~isempty(available_steps)
            n = length(available_steps);
            if n >= 4
                idx = round(linspace(1, n, 4));
            else
                idx = 1:n;
            end
            steps_to_plot = available_steps(idx);
        else
            % Only end.mat exists, no multi-step snapshots possible
            steps_to_plot = [];
            fprintf('Note: No numeric step files found, skipping evolution snapshots.\n');
        end
    end
    
    %% ===== Figure 1: Final-state field maps =====
    fprintf('Generating Fig 1: Field contour maps...\n');
    plot_field_maps(data_end, opts, outdir);
    
    %% ===== Figure 2: Time history curves =====
    fprintf('Generating Fig 2: Time histories...\n');
    plot_time_histories(data_end, opts, outdir);
    
    %% ===== Figure 3: Degradation evolution snapshots =====
    if ~isempty(steps_to_plot)
        fprintf('Generating Fig 3: Evolution snapshots...\n');
        plot_evolution_snapshots(result_folder, steps_to_plot, opts, outdir);
    else
        fprintf('Skipping Fig 3: No multi-step data available.\n');
    end
    
    %% ===== Figure 4: Crack path with deformation =====
    fprintf('Generating Fig 4: Crack path...\n');
    plot_crack_path(data_end, opts, outdir);
    
    %% ===== Figure 5: Reliability dashboard =====
    fprintf('Generating Fig 5: Reliability assessment...\n');
    plot_reliability(data_end, opts, outdir);
    
    %% ===== Figure 6: Gc degradation map =====
    fprintf('Generating Fig 6: Gc degradation map...\n');
    plot_gc_degradation(data_end, opts, outdir);
    
    fprintf('\n=== All figures saved to %s ===\n', outdir);
end
