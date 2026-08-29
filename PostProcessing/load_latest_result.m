function [data, filepath] = load_latest_result(result_folder)
    %LOAD_LATEST_RESULT Load end.mat or the highest numeric step result.

    if nargin < 1 || ~isfolder(result_folder)
        error('Result folder not found: %s', result_folder);
    end

    endfile = fullfile(result_folder, 'end.mat');
    if isfile(endfile)
        filepath = endfile;
        data = load(filepath);
        return
    end

    files = dir(fullfile(result_folder, '*.mat'));
    steps = [];
    names = {};
    for i = 1:numel(files)
        [~, name] = fileparts(files(i).name);
        step = str2double(name);
        if ~isnan(step)
            steps(end+1) = step; %#ok<AGROW>
            names{end+1} = files(i).name; %#ok<AGROW>
        end
    end

    if isempty(steps)
        error('No end.mat or numeric step .mat files found in: %s', result_folder);
    end

    [~, idx] = max(steps);
    filepath = fullfile(result_folder, names{idx});
    data = load(filepath);
end
