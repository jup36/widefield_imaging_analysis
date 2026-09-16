function MotifXcorrResidualPosLagShuffle_Spock_ArrayTask(manifest_fn)
% MotifXcorrResidualPosLagShuffle_Spock_ArrayTask
%
% Shared array-task wrapper for BOTH doPSTHSubtraction variants of the
% xcorr pipeline (PSTH-subtracted/"residual", AND raw-trace/no-subtraction
% -- see motifH_perTrial_xcorr_residual_posLag_withinTrialTimeShuff_func.m
% and Data_GNG_Scotty_MotifXcorrRsdPosLagShuff_ArraySubmit.m for the full
% rationale). Each array task reads SLURM_ARRAY_TASK_ID, runs one session
% through motifH_perTrial_xcorr_residual_posLag_withinTrialTimeShuff_func,
% and SAVES the result to that session's own Matfiles folder under a
% filename whose base keyword ("xcorrResidualPosLagShuffle" or
% "xcorrRawTracePosLagShuffle") is built by the submission script and
% handed to this wrapper pre-constructed via save_fn_bucket -- this
% wrapper itself does not need to know or care which variant it's running,
% since it just saves to whatever path it was given.
%
% NOTE: this ONE wrapper is reused for both variants (the submission
% script always points Scotty at this same function name regardless of
% doPSTHSubtraction) -- it is NOT residual-specific despite the function's
% name (kept for backward compatibility with already-deployed Scotty job
% scripts referencing it). The log messages below reflect the ACTUAL mode
% used for this task (extracted from funcNV), rather than hardcoding
% "RESIDUAL", so a raw-trace run's .out log can never be misread as a
% PSTH-subtracted run's log or vice versa.

fprintf('\n============================================================\n');
fprintf('MotifXcorrResidualPosLagShuffle_Spock_ArrayTask\n');
fprintf('Manifest: %s\n', manifest_fn);
fprintf('============================================================\n');

% Ensure shared utility functions (extract_date_animalID_header,
% GrabFiles_sort_trials, stack_trials_H, find_keyword_containing_folder,
% find_keyword_containing_files, trialTypeInfoAuditoryGngTbytDat,
% natsort/natsortfiles, compatiblepath, etc.) are on path. Each array task
% launches a fresh, non-interactive `matlab -r` session on whatever
% compute node it lands on -- it does NOT inherit path changes made in an
% interactive login-node session, which is why these were undefined even
% after appearing to be "added" manually.
%
% NOTE: this must cover the whole Widefield_Imaging_Analysis root, not
% just the Utilities subfolder (matches the raw-trace ArrayTask's fix:
% natsort.m/natsortfiles.m and stack_trials_H.m live in sibling folders
% that genpath(Utilities) alone could never reach).
widefieldRoot = '/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis';
if exist(widefieldRoot, 'dir') == 7
    addpath(genpath(widefieldRoot));
    fprintf('Added path (recursive): %s\n', widefieldRoot);
else
    warning('Expected Widefield_Imaging_Analysis root not found: %s', widefieldRoot);
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

% -------- extract doPSTHSubtraction from funcNV purely for accurate
% logging below (the actual function call just passes funcNV{:} through
% generically regardless -- this lookup does not change any behavior,
% it only prevents the log messages from hardcoding a mode label that
% could be wrong for the other variant) --------
doPSTHSubIdx = find(strcmpi(funcNV, 'doPSTHSubtraction'));
if ~isempty(doPSTHSubIdx)
    doPSTHSubtraction = funcNV{doPSTHSubIdx + 1};
    modeLabel = ternary_local(doPSTHSubtraction, 'RESIDUAL (PSTH-subtracted)', 'RAW-TRACE (no PSTH subtraction)');
else
    % Older manifest predating this parameter -- default matches the core
    % function's own default (true), so this label is accurate for any
    % manifest built before doPSTHSubtraction existed.
    modeLabel = 'RESIDUAL (PSTH-subtracted) [doPSTHSubtraction absent from funcNV -- assuming default true]';
end

fprintf('\nArray task ID: %d/%d\n', sIdx, S.nSessions);
fprintf('Session header: %s\n', header);
fprintf('Mode: %s\n', modeLabel);
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
    % NOTE function name: motifH_perTrial_xcorr_residual_posLag_withinTrialTimeShuff_func
    % (missing trailing "le" in "TimeShuffle" is intentional -- MATLAB
    % function-name length limit forced the shortened form; the .m
    % filename and the internal function declaration must match this
    % exactly or the call below will fail to resolve).
    result = motifH_perTrial_xcorr_residual_posLag_withinTrialTimeShuff_func( ...
        filePath, fileKeyword, funcNV{:});
catch ME
    fprintf(2, '\nERROR running session %s (task %d/%d, mode: %s):\n%s\n', ...
        header, sIdx, S.nSessions, modeLabel, ME.message);
    rethrow(ME);
end

fprintf('\nCompleted %s computation for %s in %.1f sec. Saving...\n', modeLabel, header, toc(tStart));

%% -------------------- Save (this is the "save instead of return" step) --------------------
save_dir = fileparts(save_fn);
if exist(save_dir, 'dir') ~= 7
    mkdir(save_dir);
end

save(save_fn, 'result', 'header', 'filePath', 'fileKeyword', 'funcNV', '-v7.3');

fprintf('Saved %s result to:\n%s\n', modeLabel, save_fn);
fprintf('\nCompleted MotifXcorrResidualPosLagShuffle_Spock_ArrayTask for session %d/%d (%s, mode: %s)\n', ...
    sIdx, S.nSessions, header, modeLabel);

end

%% ========================================================================
function out = ternary_local(cond, valTrue, valFalse)
if cond
    out = valTrue;
else
    out = valFalse;
end
end
