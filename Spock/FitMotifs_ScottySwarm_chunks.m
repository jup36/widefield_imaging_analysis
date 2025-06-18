function [swarm_id, save_fn] = FitMotifs_ScottySwarm_chunks( ...
                                    fn, dependency_id, s_conn, ...
                                    parameter_class, nChunks)
% Fit motifs on Scotty by fanning out one sbatch per data chunk
% ------------------------------------------------------------
% fn              – path to master *_processed.mat file
% dependency_id   – numeric string; combine job runs after this ID completes
% s_conn          – struct from ssh2_command_scotty('connect',...)
% parameter_class – e.g. 'general_params_dual'
% nChunks         – number of alternating train / test chunks
%
% Returns
%   swarm_id  – 1×nChunks cell of job-ID strings
%   save_fn   – paths where each *_fit_chunkX.mat will be written
% -------------------------------------------------------------

gp = loadobj(feval(parameter_class));

swarm_id = cell(1, nChunks);
save_fn  = cell(1, nChunks);

for i = 1:nChunks
    [fn_dir, fn_temp] = fileparts(fn);
    save_fn{i} = fullfile(fn_dir, sprintf('%s_fit_chunk%d.mat', fn_temp, i));

    %––– 1.  Write a Scotty-style Slurm script –––––––––––––––––––––––
    script_name = WriteBashScriptMacScotty(sprintf('motifchunk%d', i), ...
        'FitMotifs_Spock', ...
        {ConvertMacToBucketPath(fn), ...
         ConvertMacToBucketPath(save_fn{i}), ...
         i, parameter_class}, ...
        {"'%s'","'%s'","%d","'%s'"}, ...
        'sbatch_time', 600, 'sbatch_memory', 12, ...
        'sbatch_path', ...
          "/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/Spock/");

    %––– 2.  Submit with dependency on previous pipeline step ––––––––
    follow = ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts ; ' ...
        sprintf('sbatch --dependency=afterok:%s %s', dependency_id, script_name)];
    resp = ssh2_command_scotty(s_conn, follow);

    %––– 3.  Extract clean numeric job-ID ––––––––––––––––––––––––––––
    line = resp.command_result(contains(resp.command_result, "Submitted batch job"));
    if isempty(line)
        error('FitMotifs_ScottySwarm:NoJobID', ...
              'Failed to capture job-ID for chunk %d', i);
    end
    thisID = regexprep(line, '[^\d]', '');   % digits only
    swarm_id{i} = thisID;                    % store

    fprintf('Chunk %d submitted as job %s\n', i, thisID);
end

% (Optional) make a colon-joined list if caller needs a single string
% allSwarm = strjoin(swarm_id, ':');
end
