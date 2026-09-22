function GlobalDAMotifPixelXcorr_Spock_ArrayTask(manifest_fn)
%GLOBALDAMOTIFPIXELXCORR_SPOCK_ARRAYTASK
%   One SLURM array task = one session. Reads the manifest written by
%   Data_GNG_Scotty_GlobalDAMotifPixelXcorr_ArraySubmit, picks its own row
%   by SLURM_ARRAY_TASK_ID, and runs globalDA_motifPixel_xcorr_func.
%
%   The compute function does its own saving (via 'saveDir'/'saveKeyword'),
%   matching MotifXcorrResidualPosLagShuffle_Spock_ArrayTask. The manifest
%   carries the fully constructed save directory per session, so this task
%   builds no paths of its own.
%
%   MANIFEST FIELDS EXPECTED
%     fileDA_bucket{i}     <header>_green_dff_smCollect.mat
%     fileTbytDA_bucket{i} <header>_green_tbytDat_dff.mat
%     fileH_bucket{i}      refit file with hC + tbytDat
%     fileW_bucket         COMMON basis file (single char, not per session)
%     saveDir_bucket{i}    where the result goes
%     funcNV               name-value cell forwarded verbatim
%     nSessions, mIdC, headerC
%
%   Usage on Scotty:
%     GlobalDAMotifPixelXcorr_Spock_ArrayTask('/jukebox/.../manifest.mat')

tTask = tic;

%% -------------------- which array index am I? --------------------
taskIdStr = getenv('SLURM_ARRAY_TASK_ID');
if isempty(taskIdStr)
    error('SLURM_ARRAY_TASK_ID is not set -- this must run as an array task.');
end
taskId = str2double(taskIdStr);
assert(isfinite(taskId) && taskId >= 1, 'Bad SLURM_ARRAY_TASK_ID: "%s".', taskIdStr);

fprintf('\n============================================================\n');
fprintf('GlobalDAMotifPixelXcorr array task %d\n', taskId);
fprintf('Manifest: %s\n', manifest_fn);
fprintf('Host: %s | started %s\n', getenv('HOSTNAME'), datestr(now));
fprintf('============================================================\n');

%% -------------------- ensure the compute function is on the path --------------------
% The generated sbatch script adds Spock/ but not the TDR subtree where
% globalDA_motifPixel_xcorr_func lives, so the task must add it itself.
% Done here rather than in WriteBashScriptWinScotty because that helper is
% shared by every other submitter in this pipeline.
%
% which() is printed either way: if an older copy of the compute function
% is lurking elsewhere on the path, the .out log names the one that
% actually ran -- worth having, since a stale copy would silently produce
% mirrored correlograms or reject the 'saveDir' argument.
if isempty(which('globalDA_motifPixel_xcorr_func'))
    wfRoot = '/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis';
    addpath(genpath(fullfile(wfRoot, 'TDR')));
    fprintf('Added to path: %s\n', fullfile(wfRoot, 'TDR'));
end
assert(~isempty(which('globalDA_motifPixel_xcorr_func')), ...
    ['globalDA_motifPixel_xcorr_func not found on the MATLAB path even after ' ...
     'adding the TDR subtree -- check the file exists at ' ...
     '.../TDR/motifTransition/motifDAxcorr/ and is readable.']);
fprintf('Compute function: %s\n', which('globalDA_motifPixel_xcorr_func'));

%% -------------------- load manifest --------------------
M = load(manifest_fn);
req = {'fileDA_bucket','fileTbytDA_bucket','fileH_bucket','fileW_bucket', ...
       'saveDir_bucket','funcNV','nSessions','mIdC','headerC'};
for r = req
    assert(isfield(M, r{1}), 'Manifest is missing "%s".', r{1});
end
assert(taskId <= M.nSessions, ...
    'Task %d exceeds nSessions=%d in the manifest.', taskId, M.nSessions);

fileDA     = M.fileDA_bucket{taskId};
fileTbytDA = M.fileTbytDA_bucket{taskId};
fileH      = M.fileH_bucket{taskId};
fileW      = M.fileW_bucket;              % common across sessions
saveDir    = M.saveDir_bucket{taskId};

fprintf('Session %d/%d: %s  %s\n', taskId, M.nSessions, M.mIdC{taskId}, M.headerC{taskId});
fprintf('  DA     : %s\n', fileDA);
fprintf('  tbytDA : %s\n', fileTbytDA);
fprintf('  H      : %s\n', fileH);
fprintf('  W      : %s\n', fileW);
fprintf('  saveDir: %s\n', saveDir);

% Fail loudly here rather than inside the compute function: a missing input
% on the cluster is almost always a path-conversion problem, and naming the
% offending file makes that obvious from the .out log.
for f = {fileDA, fileTbytDA, fileH, fileW}
    assert(exist(f{1}, 'file') == 2, 'Input not found on the cluster:\n  %s', f{1});
end

%% -------------------- parallel pool --------------------
% The compute function's 'useParfor' only helps if a pool exists. Size it to
% the cores SLURM actually granted, not the node's total.
nCpu = str2double(getenv('SLURM_CPUS_PER_TASK'));
if ~isfinite(nCpu) || nCpu < 1, nCpu = 1; end
if nCpu > 1
    pc = parcluster('local');
    tmpDir = fullfile(tempdir, sprintf('mljob_%s_%d', getenv('SLURM_ARRAY_JOB_ID'), taskId));
    mkdir(tmpDir);
    pc.JobStorageLocation = tmpDir;   % per-task storage: concurrent array tasks
                                      % sharing the default location corrupt
                                      % each other's job files
    parpool(pc, nCpu);
    fprintf('Parallel pool: %d workers.\n', nCpu);
else
    fprintf('Running serially (SLURM_CPUS_PER_TASK=%g).\n', nCpu);
end

%% -------------------- run --------------------
S = globalDA_motifPixel_xcorr_func(fileDA, fileTbytDA, fileH, fileW, ...
    'saveDir', saveDir, ...
    M.funcNV{:});   %#ok<NASGU>

fprintf('\nTask %d finished in %.1f min.\n', taskId, toc(tTask)/60);
if isfield(S, 'meta') && isfield(S.meta, 'savedTo')
    fprintf('Result: %s\n', S.meta.savedTo);
else
    warning('Task %d: compute function returned without a savedTo path.', taskId);
end

delete(gcp('nocreate'));
end