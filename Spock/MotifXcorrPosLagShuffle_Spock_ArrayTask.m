function MotifXcorrPosLagShuffle_Spock_ArrayTask(manifest_fn)
% MotifXcorrPosLagShuffle_Spock_ArrayTask
%
% Wrapper for SLURM array motif positive-lag xcorr + within-trial
% time-shuffle. Each array task reads SLURM_ARRAY_TASK_ID, runs one
% session, and SAVES the result to that session's own Matfiles folder
% (the underlying function itself just returns a struct; saving is this
% wrapper's job, matching FitMotifs_Spock_ArrayTask's pattern).

fprintf('\n============================================================\n');
fprintf('MotifXcorrPosLagShuffle_Spock_ArrayTask\n');
fprintf('Manifest: %s\n', manifest_fn);
fprintf('============================================================\n');

% Ensure shared utility functions (extract_date_animalID_header,
% GrabFiles_sort_trials, stack_trials_H, find_keyword_containing_folder,
% find_keyword_containing_files, trialTypeInfoAuditoryGngTbytDat,
% compatiblepath, etc.) are on path. Each array task launches a fresh,
% non-interactive `matlab -r` session on whatever compute node it lands
% on -- it does NOT inherit path changes made in an interactive login-node
% session (e.g. via addpath/savepath there), which is why
% 'extract_date_animalID_header' was undefined even after it appeared to
% be "added" manually. Adding it here guarantees every task sets its own
% path, every time, regardless of node or interactive-session state.
utilitiesDir = '/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/Utilities';
if exist(utilitiesDir, 'dir') == 7
    addpath(genpath(utilitiesDir));
    fprintf('Added utilities path: %s\n', utilitiesDir);
else
    warning('Expected utilities folder not found: %s', utilitiesDir);
end

if exist(manifest_fn, 'file') ~= 2
    error('Array manifest not found: %s', manifest_fn);
end

S = load(manifest_fn);

requiredFields = {'filePath_bucket', 'save_fn_bucket', 'fileKeyword', 'funcNV', 'nSessions'};
for i = 1:numel(requiredFields)
    if ~isfield(S, requiredFields{i})
        error('Manifest missing required field: %s', requiredFields{i});
    end
end

taskIDstr = getenv('SLURM_ARRAY_TASK_ID');
if isempty(taskIDstr)
    error('SLURM_ARRAY_TASK_ID is empty. This function must be run as a SLURM array task.');
end

sIdx = str2double(taskIDstr);
if isnan(sIdx) || sIdx < 1 || sIdx > S.nSessions
    error('Invalid SLURM_ARRAY_TASK_ID=%s for nSessions=%d.', taskIDstr, S.nSessions);
end

filePath   = S.filePath_bucket{sIdx};
fileKeyword = S.fileKeyword;
save_fn    = S.save_fn_bucket{sIdx};
funcNV     = S.funcNV;

if isfield(S, 'headerC') && numel(S.headerC) >= sIdx
    header = S.headerC{sIdx};
else
    header = sprintf('task%d', sIdx);
end

fprintf('\nArray task ID: %d/%d\n', sIdx, S.nSessions);
fprintf('Session header: %s\n', header);
fprintf('filePath:  %s\n', filePath);
fprintf('save_fn:   %s\n', save_fn);

%% -------------------- Open a parpool sized to the SLURM allocation --------------------
nWorkers = str2double(getenv('SLURM_CPUS_PER_TASK'));
if isnan(nWorkers) || nWorkers < 1
    nWorkers = feature('numcores');
    fprintf('SLURM_CPUS_PER_TASK not set; defaulting to feature(''numcores'')=%d\n', nWorkers);
end

pool = gcp('nocreate');
if isempty(pool)
    fprintf('Starting parpool with %d workers...\n', nWorkers);
    parpool('local', nWorkers);
elseif pool.NumWorkers ~= nWorkers
    fprintf('Existing parpool has %d workers (requested %d); leaving as-is.\n', ...
        pool.NumWorkers, nWorkers);
end

%% -------------------- Run --------------------
tStart = tic;

try
    result = motifH_perTrial_crosscorr_posLag_withinTrialTimeShuffle_func( ...
        filePath, fileKeyword, funcNV{:});
catch ME
    fprintf(2, '\nERROR running session %s (task %d/%d):\n%s\n', ...
        header, sIdx, S.nSessions, ME.message);
    rethrow(ME);
end

fprintf('\nCompleted computation for %s in %.1f sec. Saving...\n', header, toc(tStart));

%% -------------------- Save (this is the "save instead of return" step) --------------------
save_dir = fileparts(save_fn);
if exist(save_dir, 'dir') ~= 7
    mkdir(save_dir);
end

save(save_fn, 'result', 'header', 'filePath', 'fileKeyword', 'funcNV', '-v7.3');

fprintf('Saved result to:\n%s\n', save_fn);
fprintf('\nCompleted MotifXcorrPosLagShuffle_Spock_ArrayTask for session %d/%d (%s)\n', ...
    sIdx, S.nSessions, header);

end
