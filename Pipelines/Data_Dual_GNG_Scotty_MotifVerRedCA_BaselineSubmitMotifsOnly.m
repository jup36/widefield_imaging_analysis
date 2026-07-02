function Data_Dual_GNG_Scotty_MotifVerRedCA_BaselineSubmitMotifsOnly(filePathImg, fileKeyword, varargin)
% Data_Dual_GNG_Scotty_MotifVerRedCA_BaselineSubmitMotifsOnly
%
% Motif-only submission for baseline sessions.
%
% This assumes preprocessing has already been completed and that the shared
% processed file exists under:
%
%   gp.local_bucket / gp.processing_intermediates / LKcombo / headerS
%
% Example expected processed file:
%
%   Z:\Rodent Data\Wide Field Microscopy\ExampleData\Preprocessed\LKcombo\m1873_050525\
%       m1873_050525_baseline_dayX-X_img_processed_red_dff_combined.mat
%
% This function:
%
%   1) Loads base parameter class
%   2) Generates L/K-specific parameter class
%   3) Finds the already preprocessed file
%   4) Submits motif fitting using FitMotifs_ScottySwarm_chunks
%   5) Saves a motif submission record under the LKcombo/headerS folder
%
% Example:
%
%   Data_Dual_GNG_Scotty_MotifVerRedCA_BaselineSubmitMotifsOnly( ...
%       filePath_img, ...
%       '_red_dff_combined.mat', ...
%       'L', 5, ...
%       'K', 10, ...
%       'preproc_L', 1, ...
%       'preproc_K', 1, ...
%       'sbatch_time', 600, ...
%       'sbatch_memory', 16);

%% Parse inputs
p = inputParser;

addRequired(p, 'filePathImg', @(x) ischar(x) || isstring(x));
addRequired(p, 'fileKeyword', @(x) ischar(x) || isstring(x));

addParameter(p, 'L', 10, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'K', 10, @(x) isnumeric(x) && isscalar(x) && x > 0);

% These define which shared preprocessing file to use.
% Usually L1 K1 for the shared preprocessing step.
addParameter(p, 'preproc_L', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'preproc_K', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);

addParameter(p, 'base_parameter_class', 'general_params_dual_L1_K1', @(x) ischar(x) || isstring(x));

addParameter(p, 'sbatch_time', 600, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'sbatch_memory', 16, @(x) isnumeric(x) && isscalar(x) && x > 0);

% For baseline this should usually be auto-detected as 1.
% You can override if needed.
addParameter(p, 'n_fit_chunks', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));

parse(p, filePathImg, fileKeyword, varargin{:});

filePathImg = char(p.Results.filePathImg);
fileKeyword = char(p.Results.fileKeyword);

lagVal = p.Results.L;
kVal = p.Results.K;

preproc_L = p.Results.preproc_L;
preproc_K = p.Results.preproc_K;

base_parameter_class = char(p.Results.base_parameter_class);

sbatch_time = p.Results.sbatch_time;
sbatch_memory = p.Results.sbatch_memory;

n_fit_chunks_user = p.Results.n_fit_chunks;

fprintf('\n============================================================\n');
fprintf('Baseline motif-only submission\n');
fprintf('filePathImg: %s\n', filePathImg);
fprintf('fileKeyword: %s\n', fileKeyword);
fprintf('Motif L/K: L%d K%d\n', lagVal, kVal);
fprintf('Shared preprocessing label: L%d K%d\n', preproc_L, preproc_K);
fprintf('============================================================\n');

%% Connect to Scotty
keyFile = fullfile(getenv('USERPROFILE'), '.ssh', 'id_ed25519_scotty_matlab');
s_conn = ssh2_command_scotty('connect', 'scotty', keyFile);

%% Load base params
gp = loadobj(feval(base_parameter_class));

%% Generate motif parameter class
motifParamClass = sprintf('general_params_dual_L%d_K%d', lagVal, kVal);

fprintf('\nGenerating motif parameter class: %s\n', motifParamClass);

write_lk_parameter_class(base_parameter_class, motifParamClass, lagVal, kVal);

parameter_class = motifParamClass;

fprintf('Using motif parameter class: %s\n', parameter_class);

%% Resolve shared preprocessing folder
[~, fileheader] = fileparts(filePathImg);
headerS = extract_date_animalID_header(filePathImg);

file_processed_parent_dir = fullfile(gp.local_bucket, gp.processing_intermediates, ...
    'LKcombo', headerS);

fprintf('\nExpected shared preprocessing parent directory:\n%s\n', file_processed_parent_dir);

if exist(file_processed_parent_dir, 'dir') ~= 7
    error(['Shared preprocessing parent directory does not exist.\n' ...
           'Expected:\n%s\n\n' ...
           'This usually means preprocessing did not complete or was saved somewhere else.'], ...
           file_processed_parent_dir);
end

%% Find shared processed file
file_processed = find_shared_processed_file( ...
    file_processed_parent_dir, ...
    fileheader, ...
    fileKeyword, ...
    preproc_L, ...
    preproc_K);

fprintf('\nUsing shared processed file:\n%s\n', file_processed);

%% Sanity check processed file contents
processedInfo = whos('-file', file_processed);

varNames = {processedInfo.name};

requiredVars = {'data_train', 'data_test', 'nanpxs', 'gp'};
missingVars = requiredVars(~ismember(requiredVars, varNames));

if ~isempty(missingVars)
    error('Processed file is missing required variable(s): %s\nFile:\n%s', ...
        strjoin(missingVars, ', '), file_processed);
end

trainInfo = whos('-file', file_processed, 'data_train');
testInfo  = whos('-file', file_processed, 'data_test');

if isempty(trainInfo) || isempty(testInfo)
    error('Could not inspect data_train/data_test in processed file:\n%s', file_processed);
end

trainSize = trainInfo.size;
testSize = testInfo.size;

fprintf('\nProcessed data dimensions:\n');
fprintf('  data_train: %s\n', mat2str(trainSize));
fprintf('  data_test : %s\n', mat2str(testSize));

if numel(trainSize) >= 3
    nTrainChunks = trainSize(3);
else
    nTrainChunks = 1;
end

if numel(testSize) >= 3
    nTestChunks = testSize(3);
else
    nTestChunks = 1;
end

if nTrainChunks ~= nTestChunks
    error('data_train and data_test have different chunk counts: %d vs %d', ...
        nTrainChunks, nTestChunks);
end

if isempty(n_fit_chunks_user)
    n_fit_chunks = nTrainChunks;
else
    n_fit_chunks = n_fit_chunks_user;
end

fprintf('\nBaseline motif fitting will submit %d chunk(s).\n', n_fit_chunks);

if n_fit_chunks > nTrainChunks
    error('Requested n_fit_chunks=%d, but processed file only has %d train/test chunk(s).', ...
        n_fit_chunks, nTrainChunks);
end

%% Define L/K-specific motif output directory and save header

% Folder header should come from the image/session header.
% Example:
%   m1045_122424_task_day4-8_img -> m1045_122424_base_day4-8_img
folder_header = localMakeBaselineSaveHeader(fileheader);

motif_output_dir = fullfile(file_processed_parent_dir, ...
    sprintf('%s_motif_lag%d_k%d', folder_header, lagVal, kVal));

if exist(motif_output_dir, 'dir') ~= 7
    mkdir(motif_output_dir);
end

% File save header should come from the processed file base name.
% Example:
%   m1045_122424_task_day4-8_img_processed_red_dff_combined
%   -> m1045_122424_base_day4-8_img_processed_red_dff_combined
[~, processedBase] = fileparts(file_processed);
save_header = localMakeBaselineSaveHeader(processedBase);

fprintf('\nMotif outputs will be saved under:\n%s\n', motif_output_dir);
fprintf('Motif output file header:\n%s\n', save_header);

%% Submit motif fitting
%
% Preprocessing is already complete, so submit motif jobs immediately
% with no dependency.

fprintf('\nSubmitting baseline motif fitting with no dependency because preprocessing is already complete...\n');

dependency_job_id = '';

[swarm_id, save_fn] = FitMotifs_ScottySwarm_chunks( ...
    file_processed, ...
    dependency_job_id, ...
    s_conn, ...
    parameter_class, ...
    n_fit_chunks, ...
    'save_dir', motif_output_dir, ...
    'save_header', save_header, ...
    'L', lagVal, ...
    'K', kVal, ...
    'sbatch_time', sbatch_time, ...
    'sbatch_memory', sbatch_memory);

fprintf('\nSubmitted baseline motif jobs.\n');

%% Save motif submission record
%
% Save the submission record inside the same L/K-specific motif output folder.
% Keep the original baseline fileheader in the record filename for traceability.

saveName = fullfile(motif_output_dir, ...
    sprintf('%s_baselineMotifSubmit_L%d_K%d%s', fileheader, lagVal, kVal, fileKeyword));

save(saveName, ...
    'swarm_id', ...
    'save_fn', ...
    'file_processed', ...
    'file_processed_parent_dir', ...
    'motif_output_dir', ...
    'filePathImg', ...
    'fileKeyword', ...
    'fileheader', ...
    'folder_header', ...
    'save_header', ...
    'headerS', ...
    'parameter_class', ...
    'base_parameter_class', ...
    'lagVal', ...
    'kVal', ...
    'preproc_L', ...
    'preproc_K', ...
    'n_fit_chunks', ...
    'dependency_job_id', ...
    '-v7.3');

fprintf('\nSaved baseline motif submission record:\n%s\n', saveName);

end

%% Helper: find shared processed file
function file_processed = find_shared_processed_file(file_processed_parent_dir, fileheader, fileKeyword, preproc_L, preproc_K)

% Primary current convention
candidate1 = fullfile(file_processed_parent_dir, ...
    sprintf('%s_processed%s', fileheader, fileKeyword));

% Older/commented convention
candidate2 = fullfile(file_processed_parent_dir, ...
    sprintf('%s_lag%d_k%d_processed%s', fileheader, preproc_L, preproc_K, fileKeyword));

if exist(candidate1, 'file') == 2
    file_processed = candidate1;
    return
end

if exist(candidate2, 'file') == 2
    file_processed = candidate2;
    return
end

% Fallback search
searchPattern = fullfile(file_processed_parent_dir, ...
    sprintf('%s*processed%s', fileheader, fileKeyword));

d = dir(searchPattern);

% Exclude submission records and lists
if ~isempty(d)
    names = {d.name};

    keepI = ~contains(lower(names), 'submit') & ...
            ~contains(lower(names), 'list') & ...
            ~contains(lower(names), 'motif');

    d = d(keepI);
end

if isempty(d)
    error(['Could not find shared processed file.\n\n' ...
           'Tried:\n%s\n%s\n\n' ...
           'Also searched:\n%s'], ...
           candidate1, candidate2, searchPattern);
end

if numel(d) > 1
    fprintf('\nMultiple candidate processed files found:\n');
    for i = 1:numel(d)
        fprintf('  %02d: %s\n', i, fullfile(d(i).folder, d(i).name));
    end
    error('Ambiguous processed file match. Please inspect LKcombo folder.');
end

file_processed = fullfile(d(1).folder, d(1).name);

end


%% Helper: make baseline save header
function out = localMakeBaselineSaveHeader(in)
% Convert task/baseline naming to a compact baseline-safe output label.
%
% Examples:
%   m1045_122424_task_day4-8_img      -> m1045_122424_base_day4-8_img
%   m1045_122424_baseline_day4-8_img  -> m1045_122424_base_day4-8_img
%   m1045_122424_baseline_img         -> m1045_122424_base_img

out = char(in);

out = regexprep(out, '_task_', '_base_');
out = regexprep(out, '_baseline_', '_base_');

% Fallbacks for cases like "..._task" or "..._baseline" without following underscore
out = regexprep(out, '_task$', '_base');
out = regexprep(out, '_baseline$', '_base');

end