function Scotty_RefitCuratedBasisMotifs_Swarm_JP(filePathImg, fileKeyword, basis_dir, parameter_class, save_dir, varargin)
% Camden MacDowell
% customized by Junchol Park
%
% Summary: This function commands to run 'RefitBasisMotifs.m' on Spock.
%
% INPUTS
%   filePathImg      : image file path, usually cell or char/string
%   fileKeyword      : e.g. '_red_dff_combined.mat'
%   basis_dir        : directory/file where clusterW_output is saved
%   parameter_class  : e.g. 'general_params_dual'
%   save_dir         : directory to save refit results
%
% OPTIONAL NAME-VALUE INPUTS
%   'fileProcessed'  : full path to the preprocessed file containing data_test.
%                      If provided, this bypasses automatic lookup.
%
% E.g.,
%   Scotty_RefitCuratedBasisMotifs_Swarm_JP( ...
%       filePathImg, fileKeyword, basis_dir, parameter_class, save_dir, ...
%       'fileProcessed', file_processed);

% ------------------------------------------------------------
% Parse optional inputs
% ------------------------------------------------------------
p = inputParser;
p.addParameter('fileProcessedFolder', [], @(x) isempty(x) || ischar(x) || isstring(x));
p.parse(varargin{:});

fileProcessedFolder = p.Results.fileProcessedFolder;

% ------------------------------------------------------------
% First time per session connect scotty
% ------------------------------------------------------------
keyFile = fullfile(getenv('USERPROFILE'), '.ssh', 'id_ed25519_scotty_matlab');
s_conn = ssh2_command_scotty('connect', 'scotty', keyFile);

if iscell(filePathImg)
    filePathImg = filePathImg{1};
end

filePathImg = compatiblepath(filePathImg);
basis_dir = compatiblepath(basis_dir);
save_dir = compatiblepath(save_dir);

% ------------------------------------------------------------
% Get params
% ------------------------------------------------------------
gp = loadobj(feval(parameter_class));

% ------------------------------------------------------------
% Get file_processed
% ------------------------------------------------------------
[~, fileheader] = fileparts(filePathImg);
if ~isempty(fileProcessedFolder)

    % Use user-provided preprocessed file directly
    file_processed_folder = compatiblepath(char(fileProcessedFolder));

else
    % Automatically locate preprocessed file

    % get 'file_processed' to access the train and test datasets used for
    % motif fitting
    file_processed_folder = find_keyword_containing_folder( ...
        compatiblepath([gp.local_bucket gp.processing_intermediates]), ...
        fileheader, ...
        'recursive', false);

    if numel(file_processed_folder) > 1
        warning("More than one processed folder was found!")
    end

    if iscell(file_processed_folder)
        file_processed_folder = file_processed_folder{1};
    end
end

file_processed = fullfile(file_processed_folder, ...
[fileheader, '_processed' fileKeyword]); % fileKeyword = '_red_dff_combined.mat'


% ------------------------------------------------------------
% Validate file_processed
% ------------------------------------------------------------
if ~exist(file_processed, "file")
    error("Preprocessed data not detected: %s", file_processed)
end

%% Main
temp = load(file_processed, 'data_test');

for cur_chunk = 1:size(temp.data_test, 3)

    script_name = WriteBashScriptWinScotty( ...
        sprintf('refitchunk%d', cur_chunk), ...
        'RefitCuratedBasisMotifs_JP', ...
        {ConvertWinToBucketPath(file_processed), ...
         ConvertWinToBucketPath(basis_dir), ...
         cur_chunk, ...
         parameter_class, ...
         ConvertWinToBucketPath(save_dir)}, ...
        {"'%s'", "'%s'", '%d', "'%s'", "'%s'"}, ...
        'sbatch_time', 59, ...
        'sbatch_memory', 10, ...
        'sbatch_path', "/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/Spock/");

    % Submit a job
    remoteCmd = ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts ; ' ...
        sprintf('sbatch %s', script_name)];

    resp = ssh2_command_scotty(s_conn, remoteCmd);
    disp(resp.command_result{end});

end

end