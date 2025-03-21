function Data_Pipeline_GNG_Spock_func(filePath, varargin)

p = parse_input_ctxWideImagePreprocess(filePath, varargin);
%p = parse_input_ctxWideImagePreprocess(filePath, {'redoManual', false, 'dffMethod', 'movingavg', 'winF', 30}); 

if ~ismember(p.Results.dffMethod, {'mean', 'median', 'mode', 'movingavg', 'baseline'})
    error('Unknown method was selected for dff!')
end

filePathImg = GrabFiles_sort_trials('_img', 0, {filePath});
if isempty(filePathImg)
    filePathImg  = uigetdir(filePath, 'Select a folder containing the images');
end

if exist('s_conn', 'var')~=1
    % Open ssh connection
    username = input(' Spock Username: ', 's');
    password = passcode();
    s_conn = ssh2_config('spock.princeton.edu',username,password);
end

%Add paths
%addpath(genpath('/Volumes/buschman/Rodent Data/Wide Field Microscopy/fpCNMF'));
%addpath(genpath('/Volumes/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis'));

%configure preprocessing options
opts = ConfigurePreProcessing('crop_w', 540, 'vasc_std', 2, 'save_uncorrected', 0, 'method', p.Results.dffMethod, 'method_window', p.Results.winF);

%load general params (this is for anything after preprocessing)
parameter_class = 'general_params_example';
gp = loadobj(feval(parameter_class)); %this is not needed here, but demonstrates how I load this class in other functions by just passing the string. 

%% Manual Portion
%select folders to process and grab the first file from each rec.
%EXAMPLE DATA: Select 'subfolders' and then select '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA001/DA001_072623/DA001_072623_img'
[file_list_first_stack, folder_list_raw] = GrabFiles_sort_trials('_img_', 1, filePathImg(1)); % use GrabFiles_sort_trials to sort both files and folders 

path_prepro_log = GrabFiles_sort_trials('prepro_log', 0, folder_list_raw(1));

if isempty(path_prepro_log) || p.Results.redoManual
    %Grab reference images for each. Preload so no delay between loop.
    ref_img = GetReferenceImage(file_list_first_stack{1},opts.fixed_image); % use the first frame of the first trial

    %manual allignment
    prepro_log = ManualAlignmentAdjust(ref_img,opts);
    % 1. crop (move the window and 2-click)
    % 2. midline (draw a line with two dots and 2-click)
    % 3. bregma (drop a dot)

    %mask vasculature and manual cleanup (optional)
    prepro_log.vasc_std = 1; % note that the default vasc_std is 2
    prepro_log = MaskVasculature_JP(prepro_log.cropped_alligned_img,prepro_log); % This edited function uses 'showImgAndMask' with brighter visualization of images
    close;

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
    [file_list_raw,~] = GrabFiles('.tif',0,folder_list_raw(cur_fold)); % note that there's only one file per folder in this experiment 
    [opts_list,~] = GrabFiles('prepro_log.mat',0,folder_list_raw(1)); % use opts from the 1st trial

    %% 'Spock_Preprocessing_Pipeline.m' Create spock bash script for each file and run it
    job_id = cell(1,numel(file_list_raw));
    for cur_file = 1:numel(file_list_raw)
        input_val = {ConvertMacToBucketPath(file_list_raw{cur_file}), ConvertMacToBucketPath(opts_list{1})};
        script_name = WriteBashScriptMac(sprintf('%d_%d', cur_fold, 1),'Spock_Preprocessing_Pipeline',input_val,{"'%s'","'%s'"},...
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
    
    %% 'Spock_CombineStacksBVcorrect.m' Once each folder is done, combine all the stacks and do hemocorrection
    [~,header] = fileparts(ConvertMacToBucketPath(folder_list_raw{cur_fold}));
    file_list_preprocessed{cur_fold} = [folder_list_raw{cur_fold} filesep header '_dff_combined.mat'];
    script_name = WriteBashScriptMac(sprintf('%d_%d_combine', cur_fold, cur_file), ...
        'Spock_CombineStacksBVcorrect',{ConvertMacToBucketPath(folder_list_raw{cur_fold}), ConvertMacToBucketPath(file_list_preprocessed{cur_fold}), 'general_params_example'},{"'%s'","'%s'","'%s'"});

    % Run job with dependency
    response = ssh2_command(s_conn,...
        ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory
        sprintf('sbatch --dependency=afterok:%s %s',[job_id{:}],script_name)]);    

end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = parse_input_ctxWideImagePreprocess(filePath, vargs)
    % parse input, and extract name-value pairs
    default_redoManual = false;         % logic to redo manual curation
    default_dffMethod = 'movingavg';    % default method for dff is moving average
    default_winF = 30;       % default window width (30 s) for F

    p = inputParser; % create parser object
    addRequired(p,'filePath')
    addParameter(p,'redoManual', default_redoManual)
    addParameter(p, 'dffMethod', default_dffMethod)
    addParameter(p, 'winF', default_winF)

    parse(p, filePath, vargs{:})
end



end



