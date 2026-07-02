function Data_Dual_GNG_ScottyFunc_MotifGreenDA_ArraySubmitMotifs(filePathImg, fileKeyword, varargin)
% Data_DualPipeline_GNG_ScottyFunc_MotifVerRedCA_SubmitMotifsOnly
%
% Submit motif discovery jobs for one L/K combo using an existing shared
% preprocessed file.
%
% ARRAY VERSION:Gre
%   - Submits one SLURM array job per session/L/K combo.
%   - Each array task processes one chunk.
%   - This dramatically reduces scheduler load compared with submitting
%     one independent sbatch job per chunk.
%
% This function does NOT call ProcessAndSplitDataAuditoryGng.
% It expects preprocessing to already be complete.
%
% Required helper on Scotty MATLAB path:
%   FitMotifs_Spock_ArrayTask.m

%% Parse inputs
p = inputParser;

addRequired(p, 'filePathImg', @(x) ischar(x) || isstring(x));
addRequired(p, 'fileKeyword', @(x) ischar(x) || isstring(x));

addParameter(p, 'L', 10, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'K', 15, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'preproc_L', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'preproc_K', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'parameter_class', 'general_params_dual_L10', @(x) ischar(x) || isstring(x));

addParameter(p, 'sbatch_time', 600, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'sbatch_memory', 16, @(x) isnumeric(x) && isscalar(x) && x > 0);

% New array-specific parameter
addParameter(p, 'array_throttle', 5, @(x) isnumeric(x) && isscalar(x) && x > 0);

parse(p, filePathImg, fileKeyword, varargin{:});

filePathImg = char(p.Results.filePathImg);
fileKeyword = char(p.Results.fileKeyword);

L = p.Results.L;
K = p.Results.K;
preproc_L = p.Results.preproc_L;
preproc_K = p.Results.preproc_K;

parameter_class = char(p.Results.parameter_class);

sbatch_time = p.Results.sbatch_time;
sbatch_memory = p.Results.sbatch_memory;
array_throttle = p.Results.array_throttle;

fprintf('\n============================================================\n');
fprintf('Motif-only ARRAY submission\n');
fprintf('filePathImg: %s\n', filePathImg);
fprintf('L = %d, K = %d\n', L, K);
fprintf('Shared preprocessing label: L%d K%d\n', preproc_L, preproc_K);
fprintf('Array throttle: max %d concurrent tasks per array\n', array_throttle);
fprintf('============================================================\n');

%% Connect to Scotty
keyFile = fullfile(getenv('USERPROFILE'), '.ssh', 'id_ed25519_scotty_matlab');
s_conn = ssh2_command_scotty('connect', 'scotty', keyFile);

%% Load base params and create L/K-specific parameter class
gp = loadobj(feval(parameter_class));

gp.L = L;
gp.K = K;

baseParamClass = parameter_class;
newParamClass = sprintf('general_params_dual_L%d_K%d', L, K);

write_lk_parameter_class(baseParamClass, newParamClass, L, K);

parameter_class = newParamClass;

fprintf('\nUsing motif parameter class: %s\n', parameter_class);

%% Define paths
[~, fileheader] = fileparts(filePathImg);
headerS = extract_date_animalID_header(filePathImg);

file_processed_parent_dir = fullfile(gp.local_bucket, gp.processing_intermediates, ...
    'LKcombo', headerS);

file_processed = fullfile(file_processed_parent_dir, ...
    sprintf('%s_processed%s', fileheader, fileKeyword));

if exist(file_processed, 'file') ~= 2
    error(['Shared preprocessed file not found:\n%s\n\n' ...
           'Run preprocess mode first and confirm the preprocessing job completed.'], ...
           file_processed);
end

fprintf('\nFound shared preprocessed file:\n%s\n', file_processed);

%% Load DFF list to determine nChunks
filePathImg_dffList = fullfile(file_processed_parent_dir, ...
    sprintf('%s_list%s', fileheader, fileKeyword));

if exist(filePathImg_dffList, 'file') ~= 2
    error(['Shared DFF list file not found:\n%s\n\n' ...
           'Run preprocess mode first.'], filePathImg_dffList);
end

tmp = load(filePathImg_dffList, 'fnC');

if ~isfield(tmp, 'fnC')
    error('fnC not found in shared DFF list file:\n%s', filePathImg_dffList);
end

fnC = tmp.fnC;
nChunks = size(fnC, 1);

fprintf('Detected %d chunks from shared DFF list.\n', nChunks);

%% L/K-specific motif output directory
motif_output_dir = fullfile(file_processed_parent_dir, ...
    sprintf('%s_motif_lag%d_k%d', fileheader, L, K));

if exist(motif_output_dir, 'dir') ~= 7
    mkdir(motif_output_dir);
end

fprintf('\nMotif outputs will be saved under:\n%s\n', motif_output_dir);

%% Build chunk output filenames
save_fn = cell(nChunks, 1);
swarm_motifs = cell(nChunks, 1);

[~, processedBase, processedExt] = fileparts(file_processed);

for i = 1:nChunks
    save_fn{i} = fullfile(motif_output_dir, ...
        sprintf('%s_fit_L%d_K%d_chunk%d%s', processedBase, L, K, i, processedExt));

    swarm_motifs{i} = save_fn{i};
end

%% Build array manifest
% Store bucket paths because the manifest will be read on Scotty.
file_processed_bucket = ConvertWinToBucketPath(file_processed);
save_fn_bucket = cellfun(@ConvertWinToBucketPath, save_fn, 'UniformOutput', false);

manifestName = sprintf('%s_arrayManifest_L%d_K%d%s', fileheader, L, K, fileKeyword);
manifest_fn = fullfile(motif_output_dir, manifestName);

save(manifest_fn, ...
    'file_processed', ...
    'file_processed_bucket', ...
    'save_fn', ...
    'save_fn_bucket', ...
    'parameter_class', ...
    'newParamClass', ...
    'baseParamClass', ...
    'L', ...
    'K', ...
    'preproc_L', ...
    'preproc_K', ...
    'file_processed_parent_dir', ...
    'motif_output_dir', ...
    'filePathImg_dffList', ...
    'nChunks', ...
    '-v7.3');

fprintf('\nSaved array manifest:\n%s\n', manifest_fn);

manifest_fn_bucket = ConvertWinToBucketPath(manifest_fn);

%% Create one array sbatch script
jobStem = sanitize_for_slurm(sprintf('motif_%s_L%d_K%d_ARRAY', fileheader, L, K));

script_name = WriteBashScriptWinScotty(jobStem, ...
    'FitMotifs_Spock_ArrayTask', ...
    {manifest_fn_bucket}, ...
    {"'%s'"}, ...
    'sbatch_time', sbatch_time, ...
    'sbatch_memory', sbatch_memory, ...
    'sbatch_name', jobStem, ...
    'sbatch_path', ...
      "/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/Spock/");

%% Patch the generated sbatch script to add --array
% We insert:
%   #SBATCH --array=1-nChunks%array_throttle
%
% This keeps at most array_throttle chunks running simultaneously.
arrayDirective = sprintf('#SBATCH --array=1-%d%%%d', nChunks, array_throttle);
outputDirective = sprintf('#SBATCH --output=out/%s_%%A_%%a.out', jobStem);

dynamicScriptsDir = '/jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts';

patchCmd = sprintf([ ...
    'cd %s ; ' ...
    'mkdir -p out ; ' ...
    'if ! grep -q "^#SBATCH --array" "%s"; then ' ...
    'awk ''BEGIN{done=0} /^#SBATCH/ && done==0 {print; print "%s"; done=1; next} {print}'' "%s" > "%s.tmp" && mv "%s.tmp" "%s"; ' ...
    'fi ; ' ...
    'sed -i ''s|^#SBATCH --output=.*|%s|'' "%s"'], ...
    dynamicScriptsDir, ...
    script_name, ...
    arrayDirective, ...
    script_name, script_name, script_name, script_name, ...
    outputDirective, ...
    script_name);

respPatch = ssh2_command_scotty(s_conn, patchCmd);
if isfield(respPatch, 'command_result') && ~isempty(respPatch.command_result)
    disp(respPatch.command_result);
end

%% Submit the array job once
remoteCmd = ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts ; ' ...
    sprintf('sbatch %s', script_name)];

resp = ssh2_command_scotty(s_conn, remoteCmd);

disp(resp.command_result{end});

jidLine = resp.command_result(contains(resp.command_result, "Submitted batch job"));

if isempty(jidLine)
    error('No array job ID found in sbatch response.');
end

array_job_id = extractAfter(jidLine, "Submitted batch job ");
array_job_id = regexprep(array_job_id, '[^0-9]', '');

fprintf('\nSubmitted ARRAY job %s for %d chunks.\n', array_job_id, nChunks);
fprintf('Array throttle: %d concurrent tasks max.\n', array_throttle);
fprintf('Array directive: %s\n', arrayDirective);

%% Create swarm_id-style labels for compatibility
% These are not separate job IDs anymore. They are array task IDs.
swarm_id = cell(nChunks, 1);
for i = 1:nChunks
    swarm_id{i} = sprintf('%s_%d', array_job_id, i);
end

%% Save motif submission record
[~, filePathImg_name] = fileparts(filePathImg);

saveName = sprintf('%s_motifList_L%d_K%d%s', ...
    filePathImg_name, L, K, fileKeyword);

save(fullfile(motif_output_dir, saveName), ...
    'swarm_motifs', ...
    'swarm_id', ...
    'array_job_id', ...
    'array_throttle', ...
    'arrayDirective', ...
    'script_name', ...
    'manifest_fn', ...
    'manifest_fn_bucket', ...
    'save_fn', ...
    'save_fn_bucket', ...
    'parameter_class', ...
    'newParamClass', ...
    'baseParamClass', ...
    'L', ...
    'K', ...
    'preproc_L', ...
    'preproc_K', ...
    'file_processed', ...
    'file_processed_bucket', ...
    'file_processed_parent_dir', ...
    'motif_output_dir', ...
    'filePathImg_dffList', ...
    'nChunks', ...
    '-v7.3');

fprintf('\nSaved motif ARRAY submission record:\n%s\n', fullfile(motif_output_dir, saveName));

clearvars s_conn

end