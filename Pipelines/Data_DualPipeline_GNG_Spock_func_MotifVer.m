function Data_DualPipeline_GNG_Spock_func_MotifVer(filePathImg, fileKeyword)
%This script runs fpCNMF (FitMotifs_SpockSwarm.m, FitMotifs_Spock.m) and data 
% preprocessing required for fpCNMF (ProcessAndSplitDataAuditoryGng.m)
% This script has been modified to accommodate the current imaging scheme
% for auditory Go/No-go task (Feb, 2025), where data are acquired in
% chunks (e.g., 30 chuncks). As fpCNMF uses cross-validation, 30 acquired
% chunks can be split into 15 train and test datasets in an altenating
% fashion. 

% filePathImg: a folder that contains 'dff_combined' stacks.
%  E.g. filePathImg = '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1045_jRGECO_GRABda/m1045_122424/task/m1045_122424_task_day4-8_img';  
% fileKeyword: a keyword to search for the dff matfiles. 
%  E.g. fileKeyword = '_red_dff_combined.mat'; 

if exist('s_conn', 'var')~=1
    % Open ssh connection
    username = input(' Spock Username: ', 's'); % Use PU (NOT PNI) id and password as of 11/1/24
    password = passcode();
    s_conn = ssh2_config('spock.princeton.edu',username,password);
end

%load general params (this is for anything after preprocessing)
parameter_class = 'general_params_dual';
gp = loadobj(feval(parameter_class)); %this is not needed here, but demonstrates how I load this class in other functions by just passing the string.

% Get a cell array (fnC) that contains filePaths to 'dff_combined' matfiles
file_list_dff = GrabFiles_subfolders(fileKeyword, filePathImg); % use GrabFiles_sort_trials to sort both files and folders

fnC = cell(floor(length(file_list_dff)/2), 2); % train (1st col) and test (2nd col) set directories

for f = 1:size(fnC, 1)
    fnC{f,1} = ConvertMacToBucketPath(file_list_dff{f*2-1}); % train dff stack path
    fnC{f,2} = ConvertMacToBucketPath(file_list_dff{f*2}); % test dff stack path
end

[~, fileheader] = fileparts(filePathImg{1}); 
filePathImg_dffList = fullfile(filePathImg{1}, [fileheader, '_list', fileKeyword]); 
save(filePathImg_dffList, 'fnC')

%% Deconvolution/normalization and Motif Fitting. Results in cross validated motifs
% optional restart: selected combined dff files
% file_list_preprocessed = GrabFiles('501\w*dff_combined.mat');
file_processed = [gp.local_bucket_mac gp.processing_intermediates_mac fileheader '_processed' fileKeyword]; % path in temporary folder with the train/test split data

%deconvolve and split the data 
script_name = WriteBashScriptMac(sprintf('%d', 1),'ProcessAndSplitDataAuditoryGng', ...
    {ConvertMacToBucketPath(filePathImg_dffList),...
    ConvertMacToBucketPath(file_processed), ...
    'general_params_dual'},...
    {"'%s'","'%s'","'%s'"},...
    'sbatch_time',5,'sbatch_memory',16,...
    'sbatch_path',"/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/Preprocessing/");
%Run job
response = ssh2_command(s_conn,...
    ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory
    sprintf('sbatch %s',script_name)]);
%get job id
temp_job_id = erase(response.command_result{1},'Submitted batch job '); % passed as dependency_id to ensure 'ProcessAndSplitData' is done before 'FitMotifs_Spock'

%Fit motifs and cross-validate (parallellize by chunk)
nChunks = size(fnC, 1); 
[swarm_id, swarm_motifs] = FitMotifs_SpockSwarm_chunks(file_processed, temp_job_id, s_conn, 'general_params_dual', nChunks); %  

% save swarm_id and swarm_motifs
[~, filePathImg_name] = fileparts(filePathImg{1}); 
saveName = [filePathImg_name, '_', 'motifList', fileKeyword]; 
save(fullfile(filePathImg{1}, saveName), 'swarm_motifs', 'swarm_id', '-v7.3'); 
clearvars s_conn

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function p = parse_input_ctxWideImageDualPreprocess(filePath, vargs)
        % parse input, and extract name-value pairs
        default_redoManual = false;         % logic to redo manual curation
        default_dffMethod = 'movingavg';    % default method for dff is moving average
        default_winF = 30;       % default window width (30 s) for F
        default_channelToProcess = 'both';

        p = inputParser; % create parser object
        addRequired(p,'filePath')
        addParameter(p,'redoManual', default_redoManual)
        addParameter(p, 'dffMethod', default_dffMethod)
        addParameter(p, 'winF', default_winF)
        addParameter(p, 'channelToProcess', default_channelToProcess)

        parse(p, filePath, vargs{:})
    end

    function visualizeDff(dff)
        % Assuming dff is a 64x64x860 matrix
        figureHandle = figure; % Create a figure window

        for idx = 1:size(dff, 3) % Loop through all frames
            if ~isvalid(figureHandle) % Check if the figure is closed
                disp('Visualization aborted by user.');
                break;
            end

            imagesc(dff(:, :, idx)); % Visualize the current frame
            clim([-2 2]); % Set the color limits
            colormap('parula'); % Set a colormap (optional)
            colorbar; % Add a colorbar for reference
            title(sprintf('Frame %d', idx)); % Display the frame number
            pause(0.15); % Pause for 1 second
        end
    end

end



