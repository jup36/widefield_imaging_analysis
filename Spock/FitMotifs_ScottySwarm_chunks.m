function [swarm_id, save_fn] = FitMotifs_ScottySwarm_chunks( ...
    fn, dependency_id, s_conn, ...
    parameter_class, nChunks, varargin)
% FitMotifs_ScottySwarm_chunks
%
% Fit motifs on Scotty by submitting one sbatch job per data chunk.
%
% Backward compatible:
%   Old call still works:
%
%       FitMotifs_ScottySwarm_chunks(fn, dependency_id, s_conn, parameter_class, nChunks)
%
%   New call can specify output directory/name:
%
%       FitMotifs_ScottySwarm_chunks(fn, dependency_id, s_conn, parameter_class, nChunks, ...
%           'save_dir', motif_output_dir, ...
%           'save_header', save_header, ...
%           'L', lagVal, ...
%           'K', kVal)
%
% Inputs
%   fn              - path to master *_processed.mat file
%   dependency_id   - numeric string; if empty, jobs submit immediately
%   s_conn          - struct from ssh2_command_scotty('connect',...)
%   parameter_class - e.g. 'general_params_dual_L10_K5'
%   nChunks         - number of train/test chunks
%
% Optional name-value inputs
%   save_dir        - output folder for fitted chunk files
%   save_header     - output file header/base name
%   L               - lag value for filename/job name
%   K               - K value for filename/job name
%   sbatch_time     - sbatch time in minutes
%   sbatch_memory   - sbatch memory in GB
%   script_prefix   - prefix for dynamic script/job names, default 'motif'
%
% Returns
%   swarm_id        - 1 x nChunks cell of Slurm job ID strings
%   save_fn         - 1 x nChunks cell of output paths

%% Parse optional inputs

p = inputParser;

addParameter(p, 'save_dir', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'save_header', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'L', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
addParameter(p, 'K', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
addParameter(p, 'sbatch_time', 600, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'sbatch_memory', 16, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'script_prefix', 'motif', @(x) ischar(x) || isstring(x));

parse(p, varargin{:});

save_dir      = char(p.Results.save_dir);
save_header   = char(p.Results.save_header);
L_save        = p.Results.L;
K_save        = p.Results.K;
sbatch_time   = p.Results.sbatch_time;
sbatch_memory = p.Results.sbatch_memory;
script_prefix = char(p.Results.script_prefix);

%% Basic checks

fn = char(fn);
dependency_id = char(string(dependency_id));
dependency_id = strtrim(dependency_id);
parameter_class = char(parameter_class);

if exist(fn, 'file') ~= 2
    error('Processed file does not exist:\n%s', fn);
end

if isempty(nChunks) || ~isnumeric(nChunks) || ~isscalar(nChunks) || nChunks < 1
    error('nChunks must be a positive scalar.');
end

if nChunks ~= round(nChunks)
    error('nChunks must be an integer scalar.');
end

nChunks = double(nChunks);

%% Load params once to verify parameter class exists

gp = loadobj(feval(parameter_class)); %#ok<NASGU>

%% Build output filenames once

[fn_dir, fn_temp, fn_ext] = fileparts(fn);

% Old behavior: if save_dir is not supplied, save beside processed file.
if isempty(save_dir)
    save_dir = fn_dir;
end

if exist(save_dir, 'dir') ~= 7
    mkdir(save_dir);
end

% Old behavior: if save_header is not supplied, use processed-file base name.
if isempty(save_header)
    save_header = fn_temp;
else
    % If user passed image header only, append processed/channel suffix from fn.
    %
    % Example:
    %   save_header = m1045_122424_base_day4-8_img
    %   fn_temp     = m1045_122424_task_day4-8_img_processed_red_dff_combined
    %
    % Desired:
    %   m1045_122424_base_day4-8_img_processed_red_dff_combined
    processedSuffix = regexprep(fn_temp, '^m\d{4}_\d{6}.*?_img', '');
    
    if ~contains(save_header, 'processed')
        save_header = [save_header processedSuffix];
    end
end

swarm_id = cell(1, nChunks);
save_fn  = cell(1, nChunks);

for i = 1:nChunks
    
    if ~isempty(L_save) && ~isempty(K_save)
        save_fn{i} = fullfile(save_dir, ...
            sprintf('%s_fit_L%d_K%d_chunk%d%s', ...
            save_header, L_save, K_save, i, fn_ext));
    else
        % Backward-compatible old filename convention
        save_fn{i} = fullfile(save_dir, ...
            sprintf('%s_fit_chunk%d%s', ...
            save_header, i, fn_ext));
    end
end

%% Console report

fprintf('\nFitMotifs_ScottySwarm_chunks output plan:\n');
fprintf('  Processed file: %s\n', fn);
fprintf('  Output dir    : %s\n', save_dir);
fprintf('  Parameter     : %s\n', parameter_class);
fprintf('  nChunks       : %d\n', nChunks);

if isempty(dependency_id)
    fprintf('  Dependency    : none\n');
else
    fprintf('  Dependency    : afterok:%s\n', dependency_id);
end

fprintf('\nExample output chunk file:\n%s\n', save_fn{1});

%% Submit one job per chunk

for i = 1:nChunks
    
    %% 1. Build informative dynamic script / job stem
    
    [~, saveBase] = fileparts(save_fn{i});
    
    % Example saveBase:
    %   m1045_122424_base_day4-8_img_processed_red_dff_combined_fit_L1_K1_chunk1
    %
    % Desired session stem:
    %   m1045_122424_base_day4-8_img
    sessionStem = regexprep(saveBase, '_processed.*$', '');
    
    if isempty(sessionStem) || strcmp(sessionStem, saveBase)
        [~, fnBase] = fileparts(fn);
        sessionStem = regexprep(fnBase, '_processed.*$', '');
    end
    
    if ~isempty(L_save) && ~isempty(K_save)
        scriptStem = sprintf('%s_%s_L%d_K%d_chunk%d', ...
            script_prefix, sessionStem, L_save, K_save, i);
    else
        scriptStem = sprintf('%s_%s_chunk%d', ...
            script_prefix, sessionStem, i);
    end
    
    scriptStem = sanitize_for_slurm(scriptStem);
    
    %% 2. Write a Scotty-style Slurm script
    
    script_name = WriteBashScriptWinScotty(scriptStem, ...
        'FitMotifs_Spock', ...
        {ConvertWinToBucketPath(fn), ...
        ConvertWinToBucketPath(save_fn{i}), ...
        i, parameter_class}, ...
        {"'%s'", "'%s'", "%d", "'%s'"}, ...
        'sbatch_time', sbatch_time, ...
        'sbatch_memory', sbatch_memory, ...
        'sbatch_path', ...
        "/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/Spock/");
    
    %% 3. Submit, with dependency only if provided
    
    dynamicScriptDir = '/jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts';
    
    if isempty(dependency_id)
        follow = sprintf('cd %s ; sbatch %s', ...
            dynamicScriptDir, script_name);
    else
        follow = sprintf('cd %s ; sbatch --dependency=afterok:%s %s', ...
            dynamicScriptDir, dependency_id, script_name);
    end
    
    resp = ssh2_command_scotty(s_conn, follow);
    
    %% 4. Extract clean numeric job ID
    
    if ~isfield(resp, 'command_result') || isempty(resp.command_result)
        error('FitMotifs_ScottySwarm:EmptyResponse', ...
            'Empty sbatch response for chunk %d.', i);
    end
    
    cmdResult = resp.command_result;
    
    if ischar(cmdResult)
        cmdResultC = cellstr(cmdResult);
    elseif isstring(cmdResult)
        cmdResultC = cellstr(cmdResult);
    elseif iscell(cmdResult)
        cmdResultC = cmdResult;
    else
        error('FitMotifs_ScottySwarm:BadResponseType', ...
            'Unexpected sbatch response type for chunk %d.', i);
    end
    
    lineI = contains(cmdResultC, 'Submitted batch job');
    line = cmdResultC(lineI);
    
    if isempty(line)
        fprintf(2, '\nSBATCH response for chunk %d:\n', i);
        disp(cmdResultC(:));
        
        error('FitMotifs_ScottySwarm:NoJobID', ...
            'Failed to capture job ID for chunk %d.', i);
    end
    
    thisID = regexprep(line{end}, '[^\d]', '');
    
    if isempty(thisID)
        fprintf(2, '\nSBATCH response for chunk %d:\n', i);
        disp(cmdResultC(:));
        
        error('FitMotifs_ScottySwarm:NoNumericJobID', ...
            'Failed to parse numeric job ID for chunk %d.', i);
    end
    
    swarm_id{i} = thisID;
    
    fprintf('Chunk %d submitted as job %s\n', i, thisID);
    fprintf('  Script: %s\n', script_name);
    fprintf('  Output: %s\n', save_fn{i});
end

end

%% Helper: sanitize Slurm/script names
function s = sanitize_for_slurm(s)

s = char(s);

% Slurm-safe: letters, numbers, underscore, dash only
s = regexprep(s, '[^\w\-]', '_');

% Avoid very long dynamic script/job names.
% The session name plus L/K/chunk should still remain readable.
maxLen = 120;

if numel(s) > maxLen
    s = s(1:maxLen);
end

end