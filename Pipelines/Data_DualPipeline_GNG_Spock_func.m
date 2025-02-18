function Data_DualPipeline_GNG_Spock_func(filePathImg, fileKeyword, varargin)
%This script runs preprocessing of an ome.tiff image stack.
% filePathImg must be a folder that contains ome.tiff stack(s).

p = parse_input_ctxWideImageDualPreprocess(filePathImg, varargin);
%p = parse_input_ctxWideImageDualPreprocess(filePathImg, {'redoManual', false, 'dffMethod', 'movingavg', 'winF', 30});

if ~ismember(p.Results.dffMethod, {'mean', 'median', 'mode', 'movingavg', 'baseline'})
    error('Unknown method was selected for dff!')
end

if exist('s_conn', 'var')~=1
    % Open ssh connection
    username = input(' Spock Username: ', 's'); % Use PU (NOT PNI) id and password as of 11/1/24
    password = passcode();
    s_conn = ssh2_config('spock.princeton.edu',username,password);
end

%Add paths
%addpath(genpath('/Volumes/buschman/Rodent Data/Wide Field Microscopy/fpCNMF'));
%addpath(genpath('/Volumes/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis'));

%configure preprocessing options
opts = ConfigurePreProcessingDual('crop_w', 400, 'crop_h', 400, 'x_bregma_margin', 210, 'y_bregma_margin', 160, ...
    'fps', 10, 'vasc_std', 1.5, 'save_uncorrected', 0, 'method', p.Results.dffMethod, 'method_window', p.Results.winF);
opts.x_bregma_margin = 280-(540-opts.crop_w)/2; % Adjusted to match the brain-outline mask to new image size (PNI284)
opts.y_bregma_margin = 230-(540-opts.crop_h)/2; % Adjusted to match the brain-outline mask to new image size

%load general params (this is for anything after preprocessing)
parameter_class = 'general_params_mac';
gp_mac = loadobj(feval(parameter_class)); %this is not needed here, but demonstrates how I load this class in other functions by just passing the string.

%% Manual Portion
[file_list_first_stack, folder_list_raw] = GrabFiles_subfolders(fileKeyword, filePathImg); % use GrabFiles_sort_trials to sort both files and folders

path_prepro_log = GrabFiles_sort_trials('prepro_log', 0, folder_list_raw(1));

if isempty(path_prepro_log) || p.Results.redoManual
    %Grab reference images for each. Preload so no delay between loop.
    ref_img = GetReferenceImage(file_list_first_stack{1}, opts.fixed_image); % use the first frame of the first trial (576x576)

    %manual allignment
    prepro_log = ManualAlignmentAdjust(ref_img, opts);
    % 1. crop (move the window and 2-click)
    % 2. midline (draw a line with two dots and 2-click)
    % 3. bregma (drop a dot) IMPORTANT! To better align the mask with the data place the bregma a bit more anteriorly than the joint point of the sutures

    %mask vasculature and manual cleanup (optional)
    prepro_log.vasc_std = 5; % note that the default vasc_std is 2
    prepro_log = MaskVasculature_JP(prepro_log.cropped_aligned_img, prepro_log); % This edited function uses 'showImgAndMask' with brighter visualization of images
    close; % TODO: remake the mask itself, implement mask repositioning relative to the bregma coordinates!

    %no transformation
    prepro_log.tform = [];
    prepro_log.output_size = [];

    % Skipping 'RegisterReferenceImages.m', as we're not combining data across
    % sessions, in which case registering to a common reference frame would be
    % necessary.

    %save off the options to each folder
    save([folder_list_raw{1} filesep 'prepro_log'],'prepro_log')
else
    load(fullfile(path_prepro_log{1}), 'prepro_log')
end


%% Run PreProcess on Spock
file_list_preprocessed = cell(1,numel(folder_list_raw));

for cur_fold = 1:numel(folder_list_raw)
    [file_list_raw, ~] = GrabFiles('.tif', 0, folder_list_raw(cur_fold)); % note that there's only one file per folder in this experiment
    [opts_list, ~] = GrabFiles('prepro_log.mat', 0, folder_list_raw(1)); % use opts from the 1st trial

    %% 'Spock_Preprocessing_Pipeline.m' Create spock bash script for each file and run it
    job_id = cell(1,numel(file_list_raw));
    for cur_file = 1:numel(file_list_raw)
        input_val = {ConvertMacToBucketPath(file_list_raw{cur_file}), ConvertMacToBucketPath(opts_list{1})};
        script_name = WriteBashScriptMac(sprintf('%d_%d', cur_fold, 1), ...
            'Spock_Preprocessing_Pipeline',input_val,{"'%s'","'%s'"},...
            'sbatch_time',15,'sbatch_memory',8);  % bash script to run for preprocessing

        % Run job
        response = ssh2_command(s_conn,...
            ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory
            sprintf('sbatch %s',script_name)]);

        % get job id
        job_id{cur_file} = erase(response.command_result{1},'Submitted batch job ');
        if cur_file ~= numel(file_list_raw)
            job_id{cur_file} = [job_id{cur_file} ','];
        end
    end

    %% 'Spock_CombineStacksHemoCorrect.m' Once each folder is done, combine all the stacks and do hemocorrection
    [~,header] = fileparts(ConvertMacToBucketPath(folder_list_raw{cur_fold}));
    file_list_preprocessed{cur_fold} = [folder_list_raw{cur_fold} filesep header '_dff_combined.mat'];
    script_name = WriteBashScriptMac(sprintf('%d_%d_combine', cur_fold, cur_file), ...
        'Spock_CombineStacksHemoCorrect', ...
        {ConvertMacToBucketPath(folder_list_raw{cur_fold}), ConvertMacToBucketPath(file_list_preprocessed{cur_fold}), 'general_params_dual'}, ...
        {"'%s'","'%s'","'%s'"});

    % Run job with dependency
    response = ssh2_command(s_conn,...
        ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory
        sprintf('sbatch --dependency=afterok:%s %s',[job_id{:}],script_name)]); % Note the dependency: Spock_CombineStacksHemoCorrect will run only after the completion of Spock_Preprocessing_Pipeline job

end

%% 
% %% Deconvolution/normalization and Motif Fitting. Results in cross validated motifs
% % %optional restart: selected combined dff files
% % file_list_preprocessed = GrabFiles('501\w*dff_combined.mat');
% job_id = cell(1,numel(file_list_preprocessed));
% file_list_motifs = cell(1,numel(file_list_preprocessed));
% file_list_processed = cell(1,numel(file_list_preprocessed));
% for cur_file = 1:numel(file_list_preprocessed)
%     [~, fn_temp] = fileparts(file_list_preprocessed{cur_file});
%     file_list_processed{cur_file} = [gp_mac.local_bucket gp_mac.processing_intermediates fn_temp '_processed.mat']; %path in temporary folder with the train/test split data
% 
%     %deconvolve and split the data
%     script_name = WriteBashScriptMac(sprintf('%d',cur_file),'ProcessAndSplitData', ...
%         {ConvertMacToBucketPath(file_list_preprocessed{cur_file}),...
%         ConvertMacToBucketPath(file_list_processed{cur_file}), ...
%         'general_params_example'},...
%         {"'%s'","'%s'","'%s'"},...
%         'sbatch_time',5,'sbatch_memory',16,...
%         'sbatch_path',"/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/Preprocessing/");
%     %Run job
%     response = ssh2_command(s_conn,...
%         ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory
%         sprintf('sbatch %s',script_name)]);
%     %get job id
%     temp_job_id = erase(response.command_result{1},'Submitted batch job '); % passed as dependency_id to ensure 'ProcessAndSplitData' is done before 'FitMotifs_Spock'
% 
%     %Fit motifs and cross-validate (parallellize by chunk)
%     [swarm_id, swarm_motifs] = FitMotifs_SpockSwarm(file_list_processed{cur_file},temp_job_id,s_conn,'general_params_dual');
%     job_id{cur_file} = [swarm_id{:}];
%     file_list_motifs{cur_file} = swarm_motifs;
% end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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



