function script_name = WriteBashScriptWinScotty( ...
    uniqID, func_name, input_val, input_type, varargin)
% WriteBashScriptWinScotty
% ----------------------------------------------------------
% Create a Slurm batch script for the Scotty cluster that runs
% a MATLAB function in headless mode.
%
% Backward compatible with old usages.
%
% New optional fields:
%   'sbatch_array'  : e.g. '1-15%2'. If empty, no array line is written.
%   'sbatch_output' : custom output path, e.g. 'out/job_%A_%a.out'.
%                    If empty, defaults to out/<jobName>_%j.out.
%
% Important behavior:
%   If gp.sbatch_name is empty, use uniqID as the Slurm job name.
%   This makes job names and output files informative instead of random.
% ----------------------------------------------------------

%% Sanitize main ID early

uniqID = sanitize_for_slurm(uniqID);

%% Grab defaults and parse overrides

gp = general_params_win;

% Extra sbatch-only options not defined in general_params_win
opts.sbatch_array  = '';
opts.sbatch_output = '';

% Split varargin into gp-compatible options and extra options
extraNames = {'sbatch_array', 'sbatch_output'};

varargin_gp = {};
varargin_extra = {};

i = 1;
while i <= numel(varargin)
    
    key = varargin{i};
    
    if i == numel(varargin)
        error('Optional input "%s" has no value.', string(key));
    end
    
    if any(strcmpi(char(key), extraNames))
        varargin_extra = [varargin_extra, varargin(i:i+1)]; %#ok<AGROW>
    else
        varargin_gp = [varargin_gp, varargin(i:i+1)]; %#ok<AGROW>
    end
    
    i = i + 2;
end

% Parse general_params_win-supported options.
% ParseOptionalInputs in your codebase expects the overrides as one cell array.
gp = ParseOptionalInputs(gp, varargin_gp);

% Parse extra sbatch-only options manually.
for ii = 1:2:numel(varargin_extra)
    opts.(char(varargin_extra{ii})) = varargin_extra{ii+1};
end

%% Decide Slurm job name

% If user explicitly supplied gp.sbatch_name, respect it.
% Otherwise, use uniqID so job/output names are informative.
if isfield(gp, 'sbatch_name') && ~isempty(gp.sbatch_name)
    jobName = sanitize_for_slurm(gp.sbatch_name);
else
    jobName = uniqID;
end

%% Where to save dynamic script

script_name = sprintf('run_%s.sh', uniqID);

file_path = [gp.local_bucket gp.dynamic_script_path];

if ~exist(file_path, 'dir')
    mkdir(file_path);
end

script_full_path = fullfile(file_path, script_name);

fid = fopen(script_full_path, 'wt');

if fid < 0
    error('Could not open script for writing: %s', script_full_path);
end

try
    %% ---------- Slurm header ----------
    
    fprintf(fid, '#!/bin/bash\n');
    
    fprintf(fid, '#SBATCH --job-name=%s\n', jobName);
    
    % Optional SLURM array line
    if ~isempty(opts.sbatch_array)
        fprintf(fid, '#SBATCH --array=%s\n', opts.sbatch_array);
    end
    
    % Output path: default old format unless custom output is provided.
    % Now the default uses informative jobName instead of EntertainingSpockNames.
    if ~isempty(opts.sbatch_output)
        fprintf(fid, '#SBATCH --output=%s\n', opts.sbatch_output);
    else
        fprintf(fid, '#SBATCH --output=out/%s_%%j.out\n', jobName);
    end
    
    fprintf(fid, '#SBATCH --nodes=1\n');
    fprintf(fid, '#SBATCH --ntasks=1\n');
    fprintf(fid, '#SBATCH --cpus-per-task=1\n');
    fprintf(fid, '#SBATCH --mem=%dG\n', gp.sbatch_memory);
    fprintf(fid, '#SBATCH --partition=all\n');
    
    % gp.sbatch_time is interpreted as minutes
    fprintf(fid, '#SBATCH --time=%02d:%02d:00\n', ...
        floor(gp.sbatch_time / 60), mod(gp.sbatch_time, 60));
    
    fprintf(fid, '#SBATCH --mail-type=END\n');
    fprintf(fid, '#SBATCH --mail-user=jp3025@princeton.edu\n\n');
    
    %% ---------- Runtime environment ----------
    
    fprintf(fid, 'module purge\n');
    fprintf(fid, 'module load matlab/R2021b\n');
    fprintf(fid, 'sleep 10\n\n');
    
    %% ---------- Working directory ----------
    
    fprintf(fid, 'cd "%s"\n\n', gp.sbatch_path);
    
    %% ---------- Build MATLAB call ----------
    
    argStr = cell(1, numel(input_val));
    
    for k = 1:numel(input_val)
        argStr{k} = sprintf(input_type{k}, input_val{k});
    end
    
    paramList = strjoin(string(argStr), ', ');
    
    matlabCmd = sprintf([ ...
        'matlab -nodisplay -nodesktop -nosplash -r "' ...
        'try; %s(%s); ' ...
        'catch ME; disp(getReport(ME,''extended'')); exit(1); end; ' ...
        'quit;"\n'], func_name, paramList);
    
    fprintf(fid, '%s', matlabCmd);
    
    fclose(fid);
    
    % Preserve previous line-ending behavior
    unix2dos(script_full_path, 1);
    
catch ME
    fclose(fid);
    error('Failed generating Scotty bash script:\n%s', ME.message);
end

end

%% Helper: sanitize Slurm/script names
function s = sanitize_for_slurm(s)

s = char(s);

% Slurm-safe: letters, numbers, underscore, dash only
s = regexprep(s, '[^\w\-]', '_');

% Avoid very long job/script names.
maxLen = 120;

if numel(s) > maxLen
    s = s(1:maxLen);
end

end