function Data_DualPipeline_GNG_Scotty_func_MotifVerRedCA_LKcombo(filePathImg, fileKeyword, varargin)
% Data_DualPipeline_GNG_Scotty_func_MotifVerRedCA_LKcombo
%
% This script runs fpCNMF motif discovery across user-specified L/K values.
%
% It performs:
%   1. Input parsing for L and K
%   2. Creation of an L/K-specific parameter class
%   3. Preprocessing / train-test splitting using ProcessAndSplitDataAuditoryGng
%   4. Motif fitting using FitMotifs_ScottySwarm_chunks
%
% filePathImg:
%   Folder containing dff_combined stacks.
%
% Example:
%   filePathImg = '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1045_jRGECO_GRABda/m1045_122424/task/m1045_122424_task_day4-8_img';
%
% fileKeyword:
%   Keyword to search for dff matfiles.
%
% Example:
%   fileKeyword = '_red_dff_combined.mat';
%
% Name-value inputs:
%   'L'               - motif lag length
%   'K'               - motif number / upper bound
%   'parameter_class' - base parameter class
%
% Example call:
%   Data_DualPipeline_GNG_Scotty_func_MotifVerRedCA_LKcombo( ...
%       filePathImg, '_red_dff_combined.mat', ...
%       'L', 5, 'K', 10);

%% Parse inputs
p = inputParser;

addRequired(p, 'filePathImg', @(x) ischar(x) || isstring(x));
addRequired(p, 'fileKeyword', @(x) ischar(x) || isstring(x));

addParameter(p, 'L', 10, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'K', 15, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'parameter_class', 'general_params_dual_L10', @(x) ischar(x) || isstring(x));

parse(p, filePathImg, fileKeyword, varargin{:});

filePathImg = char(p.Results.filePathImg);
fileKeyword = char(p.Results.fileKeyword);

L = p.Results.L;
K = p.Results.K;
parameter_class = char(p.Results.parameter_class);

fprintf('\n============================================================\n');
fprintf('Running red-CA motif discovery\n');
fprintf('L = %d, K = %d\n', L, K);
fprintf('Base parameter class: %s\n', parameter_class);
fprintf('filePathImg: %s\n', filePathImg);
fprintf('fileKeyword: %s\n', fileKeyword);
fprintf('============================================================\n');

%% Connect to Scotty
% First time per session connect Scotty using key-based login
keyFile = fullfile(getenv('USERPROFILE'), '.ssh', 'id_ed25519_scotty_matlab');
s_conn = ssh2_command_scotty('connect', 'scotty', keyFile);

%% Load base params and override L/K locally
gp = loadobj(feval(parameter_class));

gp.L = L;
gp.K = K;

fprintf('\nUpdated local gp values:\n');
fprintf('  gp.L = %d\n', gp.L);
fprintf('  gp.K = %d\n', gp.K);

%% Write modified parameter class file for this L/K combo
baseParamClass = parameter_class;
newParamClass = sprintf('general_params_dual_L%d_K%d', L, K);

writtenParamFiles = write_lk_parameter_class(baseParamClass, newParamClass, L, K);

parameter_class = newParamClass;

fprintf('\nCreated/updated parameter class: %s\n', parameter_class);
fprintf('Parameter class files written:\n');
for ww = 1:numel(writtenParamFiles)
    fprintf('  %s\n', writtenParamFiles{ww});
end

%% Get DFF file list
file_list_dff = GrabFiles_subfolders(fileKeyword, filePathImg);
% Use GrabFiles_sort_trials instead if chronological ordering becomes an issue.

if isempty(file_list_dff)
    error('No DFF files found using keyword %s in:\n%s', fileKeyword, filePathImg);
end

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

fnC = cell(nPairs, 2); % train, test

for f = 1:nPairs
    fnC{f,1} = ConvertWinToBucketPath(file_list_dff{f*2-1}); % train dff stack path
    fnC{f,2} = ConvertWinToBucketPath(file_list_dff{f*2});   % test dff stack path
end

%% Define file/session headers and save DFF list
[~, fileheader] = fileparts(filePathImg);
headerS = extract_date_animalID_header(filePathImg);

% Make this L/K-specific to avoid collisions across parallel submissions
filePathImg_dffList = fullfile(filePathImg, ...
    sprintf('%s_list_L%d_K%d%s', fileheader, gp.L, gp.K, fileKeyword));

save(filePathImg_dffList, ...
    'fnC', ...
    'file_list_dff', ...
    'L', ...
    'K', ...
    'parameter_class');

fprintf('\nSaved L/K-specific DFF list:\n%s\n', filePathImg_dffList);

%% Deconvolution / normalization / train-test split
folderName = sprintf('%s_motif_lag%d_k%d', fileheader, gp.L, gp.K);

file_processed_dir = fullfile(gp.local_bucket, gp.processing_intermediates, ...
    'LKcombo', headerS, folderName);

if exist(file_processed_dir, 'dir') ~= 7
    mkdir(file_processed_dir);
end

fileName = sprintf('%s_lag%d_k%d_processed%s', ...
    fileheader, gp.L, gp.K, fileKeyword);

file_processed = fullfile(file_processed_dir, fileName);

fprintf('\nProcessed data output will be:\n%s\n', file_processed);

%% Write preprocessing bash script
script_name = WriteBashScriptWinScotty(sprintf('%d', 1), ...
    'ProcessAndSplitDataAuditoryGng', ...
    {ConvertWinToBucketPath(filePathImg_dffList), ...
     ConvertWinToBucketPath(file_processed), ...
     parameter_class}, ...
    {"'%s'", "'%s'", "'%s'"}, ...
    'sbatch_time', 600, ...
    'sbatch_memory', 16, ...
    'sbatch_path', "/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/Preprocessing/");

%% Submit preprocessing job
remoteCmd = ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts ; ' ...
    sprintf('sbatch %s', script_name)];

resp = ssh2_command_scotty(s_conn, remoteCmd);

disp(resp.command_result{end});

%% Get preprocessing job ID
jidLine = resp.command_result(contains(resp.command_result, "Submitted batch job"));

if isempty(jidLine)
    error('No job ID found in sbatch response!');
end

temp_job_id = extractAfter(jidLine, "Submitted batch job ");
temp_job_id = regexprep(temp_job_id, '[^0-9]', '');

fprintf('\nSubmitted preprocessing job ID: %s\n', temp_job_id);

%% Fit motifs and cross-validate
nChunks = size(fnC, 1);

fprintf('\nSubmitting motif fitting swarm with %d chunks.\n', nChunks);
fprintf('Using parameter class: %s\n', parameter_class);

[swarm_id, swarm_motifs] = FitMotifs_ScottySwarm_chunks( ...
    file_processed, ...
    temp_job_id, ...
    s_conn, ...
    parameter_class, ...
    nChunks);

%% Save swarm IDs and motif job list
[~, filePathImg_name] = fileparts(filePathImg);

% L/K-specific save name to prevent overwrite across combinations
saveName = sprintf('%s_motifList_L%d_K%d%s', ...
    filePathImg_name, gp.L, gp.K, fileKeyword);

save(fullfile(filePathImg, saveName), ...
    'swarm_motifs', ...
    'swarm_id', ...
    'parameter_class', ...
    'newParamClass', ...
    'baseParamClass', ...
    'writtenParamFiles', ...
    'L', ...
    'K', ...
    'file_processed', ...
    'file_processed_dir', ...
    'filePathImg_dffList', ...
    'file_list_dff', ...
    '-v7.3');

fprintf('\nSaved motif swarm list:\n%s\n', fullfile(filePathImg, saveName));

%% Clean up SSH connection variable
clearvars s_conn

fprintf('\nDone submitting L = %d, K = %d for:\n%s\n', L, K, filePathImg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function writtenFiles = write_lk_parameter_class(baseParamClass, newParamClass, L, K)
        % Write a new parameter class file based on an existing parameter class.
        %
        % Writes the generated class file to:
        %   1) every folder where baseParamClass exists on the MATLAB path
        %   2) the Slurm-accessible ParameterClasses folder
        %
        % This handles the case where the local copy is shadowed by the Z: copy.

        %% Locate all visible copies of base parameter class
        baseFiles = which(baseParamClass, '-all');

        if isempty(baseFiles)
            error('Could not find base parameter class file on MATLAB path: %s', baseParamClass);
        end

        if ischar(baseFiles) || isstring(baseFiles)
            baseFiles = cellstr(baseFiles);
        end

        % Defensive cleanup in case display annotations are included
        for ii = 1:numel(baseFiles)
            thisFile = baseFiles{ii};
            pctIdx = strfind(thisFile, ' % ');
            if ~isempty(pctIdx)
                thisFile = strtrim(thisFile(1:pctIdx(1)-1));
            end
            baseFiles{ii} = thisFile;
        end

        fprintf('\nAll detected copies of %s:\n', baseParamClass);
        for ii = 1:numel(baseFiles)
            fprintf('  %s\n', baseFiles{ii});
        end

        % Use first copy as template source.
        % This matches MATLAB's active class resolution.
        templateFile = baseFiles{1};

        %% Get all folders containing baseParamClass
        baseParamDirs = cell(size(baseFiles));

        for ii = 1:numel(baseFiles)
            [baseParamDirs{ii}, ~, ~] = fileparts(baseFiles{ii});
        end

        %% Fixed Slurm-accessible ParameterClasses directory
        slurmParamDir = compatiblepath( ...
            'Z:\Rodent Data\Wide Field Microscopy\Widefield_Imaging_Analysis\ParameterClasses');

        %% Destination folders: all base-class folders + Slurm folder
        destDirs = [baseParamDirs(:); {slurmParamDir}];

        % Remove duplicate folders while preserving order
        destDirsNorm = cellfun(@(x) lower(char(x)), destDirs, 'UniformOutput', false);
        [~, uniqueIdx] = unique(destDirsNorm, 'stable');
        destDirs = destDirs(uniqueIdx);

        fprintf('\nGenerated parameter class will be written to:\n');
        for dd = 1:numel(destDirs)
            fprintf('  %s\n', destDirs{dd});
        end

        %% Check all destination folders before writing anything
        for dd = 1:numel(destDirs)

            thisDir = destDirs{dd};

            if exist(thisDir, 'dir') ~= 7
                error(['Required parameter-class output folder is not accessible:\n%s\n\n' ...
                       'Cannot write generated parameter class.'], thisDir);
            end

            testFile = fullfile(thisDir, sprintf('__write_test_%s.tmp', ...
                datestr(now, 'yyyymmdd_HHMMSSFFF')));

            fid = fopen(testFile, 'w');

            if fid == -1
                error(['Parameter-class output folder exists but is not writable:\n%s\n\n' ...
                       'Check drive mounting, permissions, or network access.'], thisDir);
            end

            fprintf(fid, 'write test\n');
            fclose(fid);

            if exist(testFile, 'file') ~= 2
                error('Write test failed unexpectedly in folder:\n%s', thisDir);
            end

            delete(testFile);

        end

        %% Read template class file
        txt = fileread(templateFile);

        %% Replace classdef name
        txt = regexprep( ...
            txt, ...
            ['classdef\s+', baseParamClass], ...
            ['classdef ', newParamClass], ...
            'once');

        %% Replace K assignment
        txt = regexprep( ...
            txt, ...
            '(\n\s*K\s*=\s*)[0-9]+(\s*[;\n])', ...
            sprintf('$1%d$2', K), ...
            'once');

        %% Replace L assignment
        txt = regexprep( ...
            txt, ...
            '(\n\s*L\s*=\s*)[0-9]+(\s*[;\n])', ...
            sprintf('$1%d$2', L), ...
            'once');

        %% Safety checks
        if ~contains(txt, ['classdef ', newParamClass])
            error('Failed to update classdef name from %s to %s.', ...
                baseParamClass, newParamClass);
        end

        if ~contains(txt, sprintf('K = %d', K)) && ~contains(txt, sprintf('K=%d', K))
            warning('Could not verify K replacement for generated class: %s', newParamClass);
        end

        if ~contains(txt, sprintf('L = %d', L)) && ~contains(txt, sprintf('L=%d', L))
            warning('Could not verify L replacement for generated class: %s', newParamClass);
        end

        %% Write generated class to all destination folders
        writtenFiles = cell(numel(destDirs), 1);

        for dd = 1:numel(destDirs)

            thisFile = fullfile(destDirs{dd}, [newParamClass, '.m']);

            fid = fopen(thisFile, 'w');

            if fid == -1
                error('Could not open generated parameter class for writing:\n%s', thisFile);
            end

            fwrite(fid, txt);
            fclose(fid);

            if exist(thisFile, 'file') ~= 2
                error('File write appeared to succeed, but file was not found afterward:\n%s', thisFile);
            end

            writtenFiles{dd} = thisFile;

        end

        %% Refresh MATLAB path/class cache
        for dd = 1:numel(destDirs)
            addpath(destDirs{dd});
        end

        rehash;

        %% Report
        fprintf('\nParameter class written successfully to:\n');
        for dd = 1:numel(writtenFiles)
            fprintf('  %s\n', writtenFiles{dd});
        end

    end

end