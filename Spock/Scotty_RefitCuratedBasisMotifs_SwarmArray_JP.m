function Scotty_RefitCuratedBasisMotifs_SwarmArray_JP(filePathImg, fileKeyword, basis_dir, parameter_class, save_dir, varargin)
% Scotty_RefitCuratedBasisMotifs_SwarmArray_JP
%
% ARRAY VERSION.
%
% Submits one SLURM array job per session for curated spatiotemporal
% basis motif refitting. Each array task processes one chunk.
%
% Required helper on Scotty MATLAB path:
%   RefitCuratedBasisMotifs_ArrayTask_JP.m
%
% Required refit function on Scotty MATLAB path:
%   RefitCuratedBasisMotifs_JP.m
%
% Required fpCNMF function on Scotty MATLAB path:
%   fpCNMF_refit.m
%
% Required WriteBashScriptWinScotty support:
%   'sbatch_array'
%   'sbatch_output'
%
% Optional:
%   'fileProcessedFolder' : manually specified processed folder
%   'sbatch_time'         : minutes
%   'sbatch_memory'       : GB
%   'array_throttle'      : max concurrent array tasks
%   'scotty_host'         : default 'scotty'
%   'dateStr'             : output date string
%   'H_init_method'       : default 'projection'
%   'H_init_prctile_cap'  : default 99.5

%% Parse optional inputs
p = inputParser;
addParameter(p, 'fileProcessedFolder', [], @(x) isempty(x) || ischar(x) || isstring(x));
addParameter(p, 'sbatch_time', 59, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'sbatch_memory', 10, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'array_throttle', 2, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'scotty_host', 'scotty', @(x) ischar(x) || isstring(x));
addParameter(p, 'dateStr', datestr(now, 'mmddyy'), @(x) ischar(x) || isstring(x));
addParameter(p, 'H_init_method', 'projection', @(x) ischar(x) || isstring(x));
addParameter(p, 'H_init_prctile_cap', 99.5, @(x) isnumeric(x) && isscalar(x));
parse(p, varargin{:});

fileProcessedFolder = p.Results.fileProcessedFolder;
sbatch_time = p.Results.sbatch_time;
sbatch_memory = p.Results.sbatch_memory;
array_throttle = p.Results.array_throttle;
scotty_host = char(p.Results.scotty_host);
dateStr = char(p.Results.dateStr);
H_init_method = char(p.Results.H_init_method);
H_init_prctile_cap = p.Results.H_init_prctile_cap;

%% Connect to Scotty
keyFile = fullfile(getenv('USERPROFILE'), '.ssh', 'id_ed25519_scotty_matlab');
s_conn = ssh2_command_scotty('connect', scotty_host, keyFile);

if iscell(filePathImg)
    filePathImg = filePathImg{1};
end

filePathImg = compatiblepath(filePathImg);
basis_dir = compatiblepath(basis_dir);
save_dir = compatiblepath(save_dir);

if exist(save_dir, 'dir') ~= 7
    mkdir(save_dir);
end

%% Get params
gp = loadobj(feval(parameter_class));

%% Locate processed file
[~, fileheader] = fileparts(filePathImg);

if ~isempty(fileProcessedFolder)

    file_processed_folder = compatiblepath(char(fileProcessedFolder));

else

    file_processed_folder = find_keyword_containing_folder( ...
        compatiblepath([gp.local_bucket gp.processing_intermediates]), ...
        fileheader, ...
        'recursive', false);

    if numel(file_processed_folder) > 1
        valPrepropI = cell2mat(cellfun(@(a) ~contains(a, 'lag'), ...
            file_processed_folder, 'UniformOutput', false));

        file_processed_folder = file_processed_folder(valPrepropI);

        warning("More than one processed folder was found. Using non-lag folder.")
    end

    if isempty(file_processed_folder)
        error("No processed folder detected for %s", fileheader)
    end

    if iscell(file_processed_folder)
        file_processed_folder = file_processed_folder{1};
    end
end

file_processed = fullfile(file_processed_folder, [fileheader, '_processed' fileKeyword]);

if ~exist(file_processed, "file")
    warning("Preprocessed data not detected: %s", file_processed)
    return
end

fprintf('\nFound processed file:\n%s\n', file_processed);

%% Determine number of chunks
temp = load(file_processed, 'data_test');
nChunks = size(temp.data_test, 3);

fprintf('Detected %d chunks for curated fixed-W refitting.\n', nChunks);

%% Build expected output filenames
[~, processedBase, processedExt] = fileparts(file_processed);

save_fn_train = cell(nChunks, 1);
save_fn_test = cell(nChunks, 1);
save_fn_train_bucket = cell(nChunks, 1);
save_fn_test_bucket = cell(nChunks, 1);

for cur_chunk = 1:nChunks

    save_fn_train{cur_chunk} = fullfile(save_dir, ...
        sprintf('%s_refitChunkCurated_%s_%dtrain.mat', processedBase, dateStr, cur_chunk));

    save_fn_test{cur_chunk} = fullfile(save_dir, ...
        sprintf('%s_refitChunkCurated_%s_%dtest.mat', processedBase, dateStr, cur_chunk));

    save_fn_train_bucket{cur_chunk} = ConvertWinToBucketPath(save_fn_train{cur_chunk});
    save_fn_test_bucket{cur_chunk} = ConvertWinToBucketPath(save_fn_test{cur_chunk});
end

%% Build array manifest
manifestName = sprintf('%s_curatedRefit_arrayManifest_%s%s', ...
    fileheader, dateStr, fileKeyword);

manifest_fn = fullfile(save_dir, manifestName);

file_processed_bucket = ConvertWinToBucketPath(file_processed);
basis_dir_bucket = ConvertWinToBucketPath(basis_dir);
save_dir_bucket = ConvertWinToBucketPath(save_dir);

save(manifest_fn, ...
    'file_processed', ...
    'file_processed_bucket', ...
    'basis_dir', ...
    'basis_dir_bucket', ...
    'save_dir', ...
    'save_dir_bucket', ...
    'save_fn_train', ...
    'save_fn_test', ...
    'save_fn_train_bucket', ...
    'save_fn_test_bucket', ...
    'parameter_class', ...
    'fileKeyword', ...
    'fileheader', ...
    'dateStr', ...
    'H_init_method', ...
    'H_init_prctile_cap', ...
    'nChunks', ...
    '-v7.3');

fprintf('\nSaved curated refit array manifest:\n%s\n', manifest_fn);

manifest_fn_bucket = ConvertWinToBucketPath(manifest_fn);

%% Create one array sbatch script
jobStem = sanitize_for_slurm_local(sprintf('curatedRefit_%s_ARRAY', fileheader));

arraySpec = sprintf('1-%d%%%d', nChunks, array_throttle);
outputSpec = sprintf('out/%s_%%A_%%a.out', jobStem);

script_name = WriteBashScriptWinScotty(jobStem, ...
    'RefitCuratedBasisMotifs_ArrayTask_JP', ...
    {manifest_fn_bucket}, ...
    {"'%s'"}, ...
    'sbatch_time', sbatch_time, ...
    'sbatch_memory', sbatch_memory, ...
    'sbatch_name', jobStem, ...
    'sbatch_array', arraySpec, ...
    'sbatch_output', outputSpec, ...
    'sbatch_path', ...
    "/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/Spock/");

fprintf('\nGenerated array sbatch script:\n%s\n', script_name);
fprintf('Array spec: %s\n', arraySpec);
fprintf('Output spec: %s\n', outputSpec);

%% Submit array job once
remoteCmd = ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts ; ' ...
    sprintf('mkdir -p out ; sbatch %s', script_name)];

resp = ssh2_command_scotty(s_conn, remoteCmd);

if isfield(resp, 'command_result') && ~isempty(resp.command_result)
    disp(resp.command_result')
end

jidLine = resp.command_result(contains(resp.command_result, "Submitted batch job"));

if isempty(jidLine)
    error('No array job ID found in sbatch response.');
end

array_job_id = extractAfter(jidLine, "Submitted batch job ");
array_job_id = regexprep(array_job_id, '[^0-9]', '');

fprintf('\nSubmitted CURATED REFIT ARRAY job %s for %d chunks.\n', array_job_id, nChunks);
fprintf('Array throttle: %d concurrent tasks max.\n', array_throttle);
fprintf('Array spec: %s\n', arraySpec);

%% Save submission record
swarm_id = cell(nChunks, 1);

for i = 1:nChunks
    swarm_id{i} = sprintf('%s_%d', array_job_id, i);
end

submissionRecordName = sprintf('%s_curatedRefit_submissionRecord_%s%s', ...
    fileheader, dateStr, fileKeyword);

save(fullfile(save_dir, submissionRecordName), ...
    'array_job_id', ...
    'array_throttle', ...
    'arraySpec', ...
    'outputSpec', ...
    'script_name', ...
    'manifest_fn', ...
    'manifest_fn_bucket', ...
    'file_processed', ...
    'file_processed_bucket', ...
    'basis_dir', ...
    'basis_dir_bucket', ...
    'save_dir', ...
    'save_dir_bucket', ...
    'save_fn_train', ...
    'save_fn_test', ...
    'save_fn_train_bucket', ...
    'save_fn_test_bucket', ...
    'swarm_id', ...
    'parameter_class', ...
    'fileKeyword', ...
    'fileheader', ...
    'dateStr', ...
    'H_init_method', ...
    'H_init_prctile_cap', ...
    'nChunks', ...
    '-v7.3');

fprintf('\nSaved curated refit submission record:\n%s\n', ...
    fullfile(save_dir, submissionRecordName));

clearvars s_conn

end

%% Local helper
function out = sanitize_for_slurm_local(in)
% Make a safe SLURM job name stem.

out = char(in);
out = regexprep(out, '[^\w\-]', '_');

% Keep job names reasonably short
maxLen = 80;
if numel(out) > maxLen
    out = out(1:maxLen);
end

end