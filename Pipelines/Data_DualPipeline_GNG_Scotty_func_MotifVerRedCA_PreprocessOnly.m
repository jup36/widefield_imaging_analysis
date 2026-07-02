function Data_DualPipeline_GNG_Scotty_func_MotifVerRedCA_PreprocessOnly(filePathImg, fileKeyword, varargin)
% Submit one preprocessing job per session.
%
% This creates a shared processed file under:
%
%   gp.processing_intermediates/LKcombo/headerS/
%
% using the naming convention:
%
%   <fileheader>_processed<fileKeyword>
%
% Example:
%   m1044_121824_task_day4-4_img_processed_red_dff_combined.mat
%
% This version handles:
%
%   TASK:
%       multiple split red DFF files
%       -> paired into train/test fnC
%       -> ProcessAndSplitDataAuditoryGng
%
%   BASELINE:
%       one continuous split red DFF file
%       -> passed directly
%       -> ProcessAndSplitData

%% Parse inputs
p = inputParser;

addRequired(p, 'filePathImg', @(x) ischar(x) || isstring(x));
addRequired(p, 'fileKeyword', @(x) ischar(x) || isstring(x));

addParameter(p, 'preproc_L', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'preproc_K', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'parameter_class', 'general_params_dual_L10', @(x) ischar(x) || isstring(x));

addParameter(p, 'data_mode', 'auto', @(x) any(strcmpi(char(x), {'auto', 'task', 'baseline'})));

addParameter(p, 'sbatch_time', 300, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'sbatch_memory', 16, @(x) isnumeric(x) && isscalar(x) && x > 0);

parse(p, filePathImg, fileKeyword, varargin{:});

filePathImg = char(p.Results.filePathImg);
fileKeyword = char(p.Results.fileKeyword);

preproc_L = p.Results.preproc_L;
preproc_K = p.Results.preproc_K;
parameter_class = char(p.Results.parameter_class);

data_mode = lower(char(p.Results.data_mode));

sbatch_time = p.Results.sbatch_time;
sbatch_memory = p.Results.sbatch_memory;

fprintf('\n============================================================\n');
fprintf('Preprocess-only submission\n');
fprintf('filePathImg: %s\n', filePathImg);
fprintf('fileKeyword: %s\n', fileKeyword);
fprintf('Preprocessing label: L%d K%d\n', preproc_L, preproc_K);
fprintf('Requested data mode: %s\n', data_mode);
fprintf('============================================================\n');

%% Resolve data mode
if strcmpi(data_mode, 'auto')
    if contains(lower(filePathImg), 'baseline')
        data_mode_resolved = 'baseline';
    else
        data_mode_resolved = 'task';
    end
else
    data_mode_resolved = data_mode;
end

fprintf('\nResolved data mode: %s\n', data_mode_resolved);

%% Connect to Scotty
keyFile = fullfile(getenv('USERPROFILE'), '.ssh', 'id_ed25519_scotty_matlab');
s_conn = ssh2_command_scotty('connect', 'scotty', keyFile);

%% Load base params and create preprocessing parameter class
gp = loadobj(feval(parameter_class));

gp.L = preproc_L;
gp.K = preproc_K;

baseParamClass = parameter_class;
preprocParamClass = sprintf('general_params_dual_L%d_K%d', preproc_L, preproc_K);

write_lk_parameter_class(baseParamClass, preprocParamClass, preproc_L, preproc_K);

parameter_class = preprocParamClass;

fprintf('\nUsing preprocessing parameter class: %s\n', parameter_class);

%% Define session headers
[~, fileheader] = fileparts(filePathImg);
headerS = extract_date_animalID_header(filePathImg);

%% Shared preprocessing parent directory
% IMPORTANT:
% This preserves your original LK-combo save structure.
file_processed_parent_dir = fullfile(gp.local_bucket, gp.processing_intermediates, ...
    'LKcombo', headerS);

if exist(file_processed_parent_dir, 'dir') ~= 7
    mkdir(file_processed_parent_dir);
end

fprintf('\nShared preprocessing parent directory:\n%s\n', file_processed_parent_dir);

%% Define shared processed file
file_processed = fullfile(file_processed_parent_dir, ...
    sprintf('%s_processed%s', fileheader, fileKeyword));

fprintf('\nShared processed file will be:\n%s\n', file_processed);

%% Get DFF file list safely
% Folder-first, file-second search.
% This avoids picking up files such as:
%   *_motifList_red_dff_combined.mat
file_list_dff = find_dff_files_from_split_red_folders(filePathImg, fileKeyword);

if isempty(file_list_dff)
    error('No valid split-folder DFF files found using keyword %s in:\n%s', ...
        fileKeyword, filePathImg);
end

fprintf('\nValid DFF files selected:\n');
for ff = 1:numel(file_list_dff)
    fprintf('  %02d: %s\n', ff, file_list_dff{ff});
end

%% Build preprocessing input according to task vs baseline
filePathImg_dffList = fullfile(file_processed_parent_dir, ...
    sprintf('%s_list%s', fileheader, fileKeyword));

switch lower(data_mode_resolved)

    case 'baseline'

        % ------------------------------------------------------------
        % BASELINE:
        % One continuous DFF file.
        % Use ProcessAndSplitData, which internally chunks train/test.
        % ------------------------------------------------------------

        if numel(file_list_dff) > 1
            warning(['Baseline mode detected, but more than one split-folder DFF file was found.\n' ...
                     'Using the first file only:\n%s'], file_list_dff{1});
        end

        filePathImg_dffSingle = file_list_dff{1};

        save(filePathImg_dffList, ...
            'filePathImg_dffSingle', ...
            'file_list_dff', ...
            'preproc_L', ...
            'preproc_K', ...
            'parameter_class', ...
            'filePathImg', ...
            'data_mode_resolved', ...
            '-v7.3');

        fprintf('\nSaved baseline DFF record:\n%s\n', filePathImg_dffList);

        preprocessFunctionName = 'ProcessAndSplitData';

        preprocessArgs = { ...
            ConvertWinToBucketPath(filePathImg_dffSingle), ...
            ConvertWinToBucketPath(file_processed), ...
            parameter_class};

        fprintf('\nBaseline preprocessing input:\n%s\n', filePathImg_dffSingle);
        fprintf('Using preprocessing function: %s\n', preprocessFunctionName);

    case 'task'

        % ------------------------------------------------------------
        % TASK:
        % Multiple split DFF files.
        % Pair files as train/test.
        % Use ProcessAndSplitDataAuditoryGng.
        % ------------------------------------------------------------

        nPairs = floor(numel(file_list_dff) / 2);

        if nPairs < 1
            error('Not enough DFF files found to make train/test pairs. Found %d file(s).\nFolder:\n%s', ...
                numel(file_list_dff), filePathImg);
        end

        if mod(numel(file_list_dff), 2) ~= 0
            warning(['Odd number of DFF files found (%d). Using maximum possible pairs (%d) ' ...
                     'and ignoring the last unpaired file:\n%s'], ...
                     numel(file_list_dff), nPairs, file_list_dff{end});
        end

        fprintf('\nFound %d DFF files. Using %d train/test pairs.\n', ...
            numel(file_list_dff), nPairs);

        fnC = cell(nPairs, 2);

        for f = 1:nPairs
            fnC{f,1} = ConvertWinToBucketPath(file_list_dff{f*2-1});
            fnC{f,2} = ConvertWinToBucketPath(file_list_dff{f*2});
        end

        save(filePathImg_dffList, ...
            'fnC', ...
            'file_list_dff', ...
            'preproc_L', ...
            'preproc_K', ...
            'parameter_class', ...
            'filePathImg', ...
            'data_mode_resolved', ...
            '-v7.3');

        fprintf('\nSaved shared DFF list:\n%s\n', filePathImg_dffList);

        preprocessFunctionName = 'ProcessAndSplitDataAuditoryGng';

        preprocessArgs = { ...
            ConvertWinToBucketPath(filePathImg_dffList), ...
            ConvertWinToBucketPath(file_processed), ...
            parameter_class};

        fprintf('Using preprocessing function: %s\n', preprocessFunctionName);

    otherwise

        error('Unknown data mode: %s', data_mode_resolved);

end

%% Submit preprocessing job
jobStem = sanitize_for_slurm(sprintf('preproc_%s', fileheader));

script_name = WriteBashScriptWinScotty(jobStem, ...
    preprocessFunctionName, ...
    preprocessArgs, ...
    {"'%s'", "'%s'", "'%s'"}, ...
    'sbatch_time', sbatch_time, ...
    'sbatch_memory', sbatch_memory, ...
    'sbatch_name', jobStem, ...
    'sbatch_path', "/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/Preprocessing/");

remoteCmd = ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts ; ' ...
    sprintf('sbatch %s', script_name)];

fprintf('\nSubmitting preprocessing job:\n%s\n', remoteCmd);

resp = ssh2_command_scotty(s_conn, remoteCmd);
disp(resp.command_result{end});

jidLine = resp.command_result(contains(resp.command_result, "Submitted batch job"));

if isempty(jidLine)
    error('No job ID found in sbatch response!');
end

preproc_job_id = extractAfter(jidLine, "Submitted batch job ");
preproc_job_id = regexprep(preproc_job_id, '[^0-9]', '');

fprintf('\nSubmitted preprocessing job ID: %s\n', preproc_job_id);

%% Save preprocessing submission record
saveName = fullfile(file_processed_parent_dir, ...
    sprintf('%s_preprocessSubmit_L%d_K%d%s', fileheader, preproc_L, preproc_K, fileKeyword));

saveVars = { ...
    'preproc_job_id', ...
    'file_processed', ...
    'file_processed_parent_dir', ...
    'filePathImg_dffList', ...
    'file_list_dff', ...
    'parameter_class', ...
    'preproc_L', ...
    'preproc_K', ...
    'filePathImg', ...
    'data_mode_resolved', ...
    'preprocessFunctionName'};

if strcmpi(data_mode_resolved, 'task')
    save(saveName, saveVars{:}, 'fnC', '-v7.3');
else
    save(saveName, saveVars{:}, 'filePathImg_dffSingle', '-v7.3');
end

fprintf('\nSaved preprocessing submission record:\n%s\n', saveName);

clearvars s_conn

fprintf('\nPreprocess-only submission complete.\n');
fprintf('Expected processed file:\n%s\n', file_processed);
fprintf('============================================================\n\n');

end

%% Helper: sanitize Slurm names
function s = sanitize_for_slurm(s)
% Make a string safe for Slurm job names and script names.

s = char(s);
s = regexprep(s, '[^\w\-]', '_');

maxLen = 80;
if numel(s) > maxLen
    s = s(1:maxLen);
end

end

%% Helper: folder-first DFF search
function fileC = find_dff_files_from_split_red_folders(filePathImg, fileKeyword)
% Find DFF files in two steps:
%
%   1) Find split red folders directly under filePathImg.
%      Expected examples:
%
%         m1045_122424_baseline_1_red
%         m1045_122424_task_day4-8_1_red
%         m1045_122424_task_day4-8_2_red
%
%   2) Inside each split red folder, find files matching:
%
%         *fileKeyword
%
% This prevents accidental inclusion of files such as:
%
%         *_motifList_red_dff_combined.mat
%
% saved directly under filePathImg.

fileC = {};

if exist(filePathImg, 'dir') ~= 7
    error('filePathImg does not exist:\n%s', filePathImg);
end

%% Step 1: find direct child folders
d = dir(filePathImg);
d = d([d.isdir]);

folderNames = {d.name};
folderNames = folderNames(~ismember(folderNames, {'.', '..'}));

folderPaths = cellfun(@(x) fullfile(filePathImg, x), folderNames, 'UniformOutput', false);

% Strong expected pattern:
%     *_<number>_red
splitRedI = cellfun(@(x) ~isempty(regexp(x, '_\d+_red$', 'once')), folderNames);

splitRedFolders = folderPaths(splitRedI);

% Fallback: direct folders containing "_red"
if isempty(splitRedFolders)
    fprintf('\nNo folders matched *_<number>_red under:\n%s\n', filePathImg);
    fprintf('Falling back to direct child folders containing "_red".\n');

    splitRedI = cellfun(@(x) contains(lower(x), '_red'), folderNames);
    splitRedFolders = folderPaths(splitRedI);
end

if isempty(splitRedFolders)
    warning('No split red folders found under:\n%s', filePathImg);
    return
end

try
    splitRedFolders = sort_nat(splitRedFolders);
catch
    splitRedFolders = sort(splitRedFolders);
end

fprintf('\nCandidate split red folders:\n');
for i = 1:numel(splitRedFolders)
    fprintf('  %02d: %s\n', i, splitRedFolders{i});
end

%% Step 2: find DFF inside each split red folder only
for i = 1:numel(splitRedFolders)

    thisFolder = splitRedFolders{i};

    dffD = dir(fullfile(thisFolder, ['*' fileKeyword]));

    if isempty(dffD)
        fprintf('  No DFF file found in split folder:\n  %s\n', thisFolder);
        continue
    end

    dffNames = {dffD.name};

    % Extra safety exclusions
    keepI = ~contains(lower(dffNames), 'motiflist') & ...
            ~contains(lower(dffNames), 'processed') & ...
            ~contains(lower(dffNames), 'refitchunks') & ...
            ~contains(lower(dffNames), 'preprocesssubmit');

    dffD = dffD(keepI);

    if isempty(dffD)
        fprintf('  Only excluded files found in split folder:\n  %s\n', thisFolder);
        continue
    end

    dffFull = arrayfun(@(x) fullfile(x.folder, x.name), dffD, 'UniformOutput', false);

    try
        dffFull = sort_nat(dffFull);
    catch
        dffFull = sort(dffFull);
    end

    if numel(dffFull) > 1
        warning('Multiple DFF files found in split folder. Using first after sorting:\n%s', thisFolder);
    end

    fileC{end+1, 1} = dffFull{1};

end

%% Final cleanup
fileC = fileC(~cellfun(@isempty, fileC));

fileC = fileC(cellfun(@(f) exist(f, 'file') == 2, fileC));

[~, ia] = unique(fileC, 'stable');
fileC = fileC(sort(ia));

try
    fileC = sort_nat(fileC);
catch
    fileC = sort(fileC);
end

fprintf('\nDFF files selected from split red folders:\n');
for i = 1:numel(fileC)
    fprintf('  %02d: %s\n', i, fileC{i});
end

end

% function Data_DualPipeline_GNG_Scotty_func_MotifVerRedCA_PreprocessOnly(filePathImg, fileKeyword, varargin)
% % Submit one preprocessing job per session.
% %
% % This creates a shared processed file under:
% %
% %   gp.processing_intermediates/LKcombo/headerS/
% %
% % using the naming convention:
% %
% %   <fileheader>_lag<preproc_L>_k<preproc_K>_processed<fileKeyword>
% %
% % Example:
% %   m1044_121824_task_day4-4_img_lag1_k1_processed_red_dff_combined.mat
% 
% %% Parse inputs
% p = inputParser;
% 
% addRequired(p, 'filePathImg', @(x) ischar(x) || isstring(x));
% addRequired(p, 'fileKeyword', @(x) ischar(x) || isstring(x));
% 
% addParameter(p, 'preproc_L', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
% addParameter(p, 'preproc_K', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
% addParameter(p, 'parameter_class', 'general_params_dual_L10', @(x) ischar(x) || isstring(x));
% addParameter(p, 'sbatch_time', 300, @(x) isnumeric(x) && isscalar(x) && x > 0);
% addParameter(p, 'sbatch_memory', 16, @(x) isnumeric(x) && isscalar(x) && x > 0);
% 
% parse(p, filePathImg, fileKeyword, varargin{:});
% 
% filePathImg = char(p.Results.filePathImg);
% fileKeyword = char(p.Results.fileKeyword);
% 
% preproc_L = p.Results.preproc_L;
% preproc_K = p.Results.preproc_K;
% parameter_class = char(p.Results.parameter_class);
% 
% sbatch_time = p.Results.sbatch_time;
% sbatch_memory = p.Results.sbatch_memory;
% 
% fprintf('\n============================================================\n');
% fprintf('Preprocess-only submission\n');
% fprintf('filePathImg: %s\n', filePathImg);
% fprintf('fileKeyword: %s\n', fileKeyword);
% fprintf('Preprocessing label: L%d K%d\n', preproc_L, preproc_K);
% fprintf('============================================================\n');
% 
% %% Connect to Scotty
% keyFile = fullfile(getenv('USERPROFILE'), '.ssh', 'id_ed25519_scotty_matlab');
% s_conn = ssh2_command_scotty('connect', 'scotty', keyFile);
% 
% %% Load base params and create preprocessing parameter class
% gp = loadobj(feval(parameter_class));
% 
% gp.L = preproc_L;
% gp.K = preproc_K;
% 
% baseParamClass = parameter_class;
% preprocParamClass = sprintf('general_params_dual_L%d_K%d', preproc_L, preproc_K);
% 
% write_lk_parameter_class(baseParamClass, preprocParamClass, preproc_L, preproc_K);
% 
% parameter_class = preprocParamClass;
% 
% fprintf('\nUsing preprocessing parameter class: %s\n', parameter_class);
% 
% %% Get DFF file list
% file_list_dff = GrabFiles_subfolders(fileKeyword, filePathImg);
% 
% if isempty(file_list_dff)
%     error('No DFF files found using keyword %s in:\n%s', fileKeyword, filePathImg);
% end
% 
% nPairs = floor(numel(file_list_dff) / 2);
% 
% if nPairs < 1
%     error('Not enough DFF files found to make train/test pairs. Found %d file(s).\nFolder:\n%s', ...
%         numel(file_list_dff), filePathImg);
% end
% 
% if mod(numel(file_list_dff), 2) ~= 0
%     warning(['Odd number of DFF files found (%d). Using maximum possible pairs (%d) ' ...
%              'and ignoring the last unpaired file:\n%s'], ...
%              numel(file_list_dff), nPairs, file_list_dff{end});
% end
% 
% fprintf('\nFound %d DFF files. Using %d train/test pairs.\n', ...
%     numel(file_list_dff), nPairs);
% 
% fnC = cell(nPairs, 2);
% 
% for f = 1:nPairs
%     fnC{f,1} = ConvertWinToBucketPath(file_list_dff{f*2-1});
%     fnC{f,2} = ConvertWinToBucketPath(file_list_dff{f*2});
% end
% 
% %% Define session headers
% [~, fileheader] = fileparts(filePathImg);
% headerS = extract_date_animalID_header(filePathImg);
% 
% %% Shared preprocessing parent directory
% file_processed_parent_dir = fullfile(gp.local_bucket, gp.processing_intermediates, ...
%     'LKcombo', headerS);
% 
% if exist(file_processed_parent_dir, 'dir') ~= 7
%     mkdir(file_processed_parent_dir);
% end
% 
% %% Save shared DFF list in parent directory
% filePathImg_dffList = fullfile(file_processed_parent_dir, ...
%     sprintf('%s_list%s', fileheader, fileKeyword));
% 
% save(filePathImg_dffList, ...
%     'fnC', ...
%     'file_list_dff', ...
%     'preproc_L', ...
%     'preproc_K', ...
%     'parameter_class', ...
%     'filePathImg', ...
%     '-v7.3');
% 
% fprintf('\nSaved shared DFF list:\n%s\n', filePathImg_dffList);
% 
% %% Define shared processed file
% file_processed = fullfile(file_processed_parent_dir, ...
%     sprintf('%s_processed%s', fileheader, fileKeyword));
% 
% fprintf('\nShared processed file will be:\n%s\n', file_processed);
% 
% %% Submit preprocessing job
% jobStem = sanitize_for_slurm(sprintf('preproc_%s', fileheader));
% 
% script_name = WriteBashScriptWinScotty(jobStem, ...
%     'ProcessAndSplitDataAuditoryGng', ...
%     {ConvertWinToBucketPath(filePathImg_dffList), ...
%      ConvertWinToBucketPath(file_processed), ...
%      parameter_class}, ...
%     {"'%s'", "'%s'", "'%s'"}, ...
%     'sbatch_time', sbatch_time, ...
%     'sbatch_memory', sbatch_memory, ...
%     'sbatch_name', jobStem, ...
%     'sbatch_path', "/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/Preprocessing/");
% 
% remoteCmd = ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts ; ' ...
%     sprintf('sbatch %s', script_name)];
% 
% resp = ssh2_command_scotty(s_conn, remoteCmd);
% disp(resp.command_result{end});
% 
% jidLine = resp.command_result(contains(resp.command_result, "Submitted batch job"));
% 
% if isempty(jidLine)
%     error('No job ID found in sbatch response!');
% end
% 
% preproc_job_id = extractAfter(jidLine, "Submitted batch job ");
% preproc_job_id = regexprep(preproc_job_id, '[^0-9]', '');
% 
% fprintf('\nSubmitted preprocessing job ID: %s\n', preproc_job_id);
% 
% %% Save preprocessing submission record
% saveName = fullfile(file_processed_parent_dir, ...
%     sprintf('%s_preprocessSubmit_L%d_K%d%s', fileheader, preproc_L, preproc_K, fileKeyword));
% 
% save(saveName, ...
%     'preproc_job_id', ...
%     'file_processed', ...
%     'file_processed_parent_dir', ...
%     'filePathImg_dffList', ...
%     'fnC', ...
%     'file_list_dff', ...
%     'parameter_class', ...
%     'preproc_L', ...
%     'preproc_K', ...
%     'filePathImg', ...
%     '-v7.3');
% 
% fprintf('\nSaved preprocessing submission record:\n%s\n', saveName);
% 
% clearvars s_conn
% 
% end
% 
% 
% function s = sanitize_for_slurm(s)
% % Make a string safe for Slurm job names and script names.
% 
% s = char(s);
% s = regexprep(s, '[^\w\-]', '_');
% 
% % Slurm job names can get annoying if very long
% maxLen = 80;
% if numel(s) > maxLen
%     s = s(1:maxLen);
% end
% 
% end