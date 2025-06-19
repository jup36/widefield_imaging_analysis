function script_name = WriteBashScriptMacScotty( ...
    uniqID, func_name, input_val, input_type, varargin)
% WriteBashScriptMacScotty
% ----------------------------------------------------------
% Create a Slurm batch script (bash) for the **Scotty** cluster
% that runs a MATLAB function in fully headless mode.
%
% Syntax
%   script_name = WriteBashScriptMacScotty(uniqID,func_name,...
%                   input_val,input_type,Name,Value,...)
%
% Required
%   uniqID      – string, becomes part of the .slurm filename
%   func_name   – string, MATLAB function to call on Scotty
%   input_val   – cell array of values to feed into func_name
%   input_type  – cell array of sprintf formats (e.g. '"%s"', '%d')
%
% Name-Value overrides (defaults in general_params_mac):
%   'sbatch_time'     – minutes    (default 15)
%   'sbatch_memory'   – gigabytes  (default 8)
%   'sbatch_matlabversion' – e.g. 'R2021b' (default)
%   'sbatch_name'     – job name (otherwise an “entertaining” one)
%   'sbatch_path'     – working dir on Scotty login node
%
% Returns
%   script_name – filename (not full path) of the generated script
%
% The script is written to:
%   fullfile( gp.local_bucket, gp.dynamic_script_path, script_name )
%
% ----------------------------------------------------------
% 2025-06-17  J.Park  •  Adapted for Scotty
% ----------------------------------------------------------

% Grab defaults and parse overrides
gp = general_params_mac;                    % your existing helper
gp = ParseOptionalInputs(gp,varargin);

% Where to save
script_name = sprintf('run_%s.sh', uniqID);
file_path = [gp.local_bucket gp.dynamic_script_path];
if ~exist(file_path, 'dir'), mkdir(file_path); end
fid = fopen(fullfile(file_path, script_name), 'wt');

try
    %% ---------- Slurm header ----------
    fprintf(fid,'#!/bin/bash\n');
    if isempty(gp.sbatch_name)
        jobName = EntertainingSpockNames;   % #ok<NASGU> (still fun!)
    else
        jobName = gp.sbatch_name;
    end
    fprintf(fid,'#SBATCH --job-name=%s\n', jobName);
    fprintf(fid,'#SBATCH --output=out/%s_%%j.out\n', jobName);
    fprintf(fid,'#SBATCH --nodes=1\n');
    fprintf(fid,'#SBATCH --ntasks=1\n');
    fprintf(fid,'#SBATCH --cpus-per-task=1\n');
    fprintf(fid,'#SBATCH --mem=%dG\n', gp.sbatch_memory);
    fprintf(fid,'#SBATCH --partition=all\n');
    fprintf(fid,'#SBATCH --time=%02d:%02d:00\n', ...
        floor(gp.sbatch_time/60), mod(gp.sbatch_time,60));
    fprintf(fid,'#SBATCH --mail-type=END\n');
    fprintf(fid,'#SBATCH --mail-user=<jp3025@princeton.edu>\n\n');

    %% ---------- Runtime environment ----------
    fprintf(fid,'module purge\n');
    fprintf(fid,'module load matlab/%s\n', 'R2021b');
    fprintf(fid,'sleep 10\n\n');

    %% ---------- Working directory ----------
    fprintf(fid,'cd "%s"\n\n', gp.sbatch_path);

    %% ---------- Build MATLAB call (robust version) -------------------------
    argStr = cell(1, numel(input_val));           % pre-allocate
    for k = 1:numel(input_val)
        argStr{k} = sprintf(input_type{k}, input_val{k});
    end
    paramList = strjoin(string(argStr), ', ');    % join as comma-separated list

    matlabCmd = sprintf([ ...
        'matlab -nodisplay -nodesktop -nosplash -r "' ...
        'try; %s(%s); ' ...
        'catch ME; disp(getReport(ME,''extended'')); exit(1); end; ' ...
        'quit;"\n'], func_name, paramList);

    fprintf(fid, '%s', matlabCmd);                % write exactly one MATLAB line

    fclose(fid);

    % Optional: convert line endings for macOS <-> Linux harmony
    unix2dos(fullfile(file_path, script_name), 1);

catch ME
    fclose(fid);
    error('Failed generating Scotty bash script:\n%s', ME.message);
end
end
