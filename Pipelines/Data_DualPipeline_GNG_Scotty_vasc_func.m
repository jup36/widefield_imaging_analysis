function Data_DualPipeline_GNG_Scotty_vasc_func(filePathImg, fileKeyword, varargin)
%This script runs preprocessing of an ome.tiff image stack.
% filePathImg must be a folder that contains ome.tiff stack(s).

p = parse_input_ctxWideImageDualPreprocess(filePathImg, varargin);
% Example:
% p = parse_input_ctxWideImageDualPreprocess(filePathImg, ...
%     {'redoManual', false, 'dffMethod', 'movingavg', 'winF', 30, ...
%      'byPassManualAlign', true});

if ~ismember(p.Results.dffMethod, {'mean', 'median', 'mode', 'movingavg', 'baseline'})
    error('Unknown method was selected for dff!')
end

% First time per session connect scotty
s_conn = ssh2_command_scotty('connect','scotty');   % no key path needed

%configure preprocessing options
opts = ConfigurePreProcessingDual('crop_w', 400, 'crop_h', 400, 'x_bregma_margin', 210, 'y_bregma_margin', 160, ...
    'fps', 10, 'vasc_std', 1.5, 'save_uncorrected', 0, 'method', p.Results.dffMethod, 'method_window', p.Results.winF);
opts.x_bregma_margin = 280-(540-opts.crop_w)/2; % Adjusted to match the brain-outline mask to new image size (PNI284)
opts.y_bregma_margin = 230-(540-opts.crop_h)/2; % Adjusted to match the brain-outline mask to new image size

%% Manual Portion
[file_list_first_stack, folder_list_raw] = GrabFiles_subfolders(fileKeyword, filePathImg); % use GrabFiles_sort_trials to sort both files and folders
path_prepro_log = GrabFiles_sort_trials('prepro_log', 0, folder_list_raw(1));

if isempty(path_prepro_log) || p.Results.redoManual
    % Grab reference image from first stack
    ref_img = GetReferenceImage(file_list_first_stack{1}, opts.fixed_image); % use the first frame of the first trial (576x576)

    % manual alignment
    prepro_log = ManualAlignmentAdjust(ref_img, opts);
    % 1. crop (move the window and 2-click)
    % 2. midline (draw a line with two dots and 2-click)
    % 3. bregma (drop a dot)
    % IMPORTANT! To better align the mask with the data place the bregma a
    % bit more anteriorly than the joint point of the sutures

    % mask vasculature and manual cleanup
    prepro_log.vasc_std = 10; % note that the default vasc_std is 2
    prepro_log = MaskVasculature_JP(prepro_log.cropped_aligned_img, prepro_log);
    close;

    % no transformation
    prepro_log.tform = [];
    prepro_log.output_size = [];

    % save options to first folder
    save([folder_list_raw{1} filesep 'prepro_log'],'prepro_log')

elseif p.Results.byPassManualAlign && ~isempty(path_prepro_log)
    % Load existing prepro_log and bypass ManualAlignmentAdjust
    load(fullfile(path_prepro_log{1}), 'prepro_log')

    % Re-open vasculature definition step only
    prepro_log.vasc_std = 10; % keep same behavior as manual branch
    prepro_log = MaskVasculature_JP(prepro_log.cropped_aligned_img, prepro_log);
    close;

    % preserve transform fields
    if ~isfield(prepro_log, 'tform')
        prepro_log.tform = [];
    end
    if ~isfield(prepro_log, 'output_size')
        prepro_log.output_size = [];
    end

    % save updated prepro_log back
    save(fullfile(folder_list_raw{1}, 'prepro_log'), 'prepro_log')

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
        input_val = {ConvertWinToBucketPath(file_list_raw{cur_file}), ConvertWinToBucketPath(opts_list{1})};
        script_name = WriteBashScriptWinScotty(sprintf('%d_%d', cur_fold, 1), ...
            'Spock_Preprocessing_Pipeline',input_val,{"'%s'","'%s'"});  % bash script to run for preprocessing

        % Submit a job
        remoteCmd = ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts ; ' ...
            sprintf('sbatch %s', script_name)];

        resp = ssh2_command_scotty(s_conn, remoteCmd);
        disp(resp.command_result{end});

        % get job id        
        jidLine = resp.command_result(contains(resp.command_result, "Submitted batch job"));

        if isempty(jidLine)
            error('No job ID found in sbatch response!');
        end

        jobID = extractAfter(jidLine, "Submitted batch job ");
        jobID = regexprep(jobID, '[^0-9]', '');
        job_id{cur_file} = jobID;
    end

    %% 'Spock_CombineStacksHemoCorrect1.m'
    [~,header] = fileparts(ConvertWinToBucketPath(folder_list_raw{cur_fold}));
    file_list_preprocessed{cur_fold} = [folder_list_raw{cur_fold} filesep header '_dff_combined.mat'];
    script_name = WriteBashScriptWinScotty(sprintf('%d_%d_combine', cur_fold, cur_file), ...
        'Spock_CombineStacksHemoCorrect2', ...
        {ConvertWinToBucketPath(folder_list_raw{cur_fold}), ConvertWinToBucketPath(file_list_preprocessed{cur_fold}), 'general_params_dual'}, ...
        {"'%s'","'%s'","'%s'"});

    follow = ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts ; ' ...
        sprintf('sbatch --dependency=afterok:%s %s', job_id{cur_file}, script_name)];
    resp = ssh2_command_scotty(s_conn, follow);
    disp(resp.command_result{end});
end
clearvars s_conn

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function p = parse_input_ctxWideImageDualPreprocess(filePath, vargs)
        % parse input, and extract name-value pairs
        default_redoManual = false;         % logic to redo manual curation
        default_dffMethod = 'movingavg';    % default method for dff is moving average
        default_winF = 30;                  % default window width (30 s) for F
        default_channelToProcess = 'both';
        default_byPassManualAlign = false;  % if true, skip ManualAlignmentAdjust and redo vasculature only

        p = inputParser;
        addRequired(p,'filePath')
        addParameter(p,'redoManual', default_redoManual)
        addParameter(p, 'dffMethod', default_dffMethod)
        addParameter(p, 'winF', default_winF)
        addParameter(p, 'channelToProcess', default_channelToProcess)
        addParameter(p, 'byPassManualAlign', default_byPassManualAlign)

        parse(p, filePath, vargs{:})
    end

    function visualizeDff(dff)
        figureHandle = figure;

        for idx = 1:size(dff, 3)
            if ~isvalid(figureHandle)
                disp('Visualization aborted by user.');
                break;
            end

            imagesc(dff(:, :, idx));
            clim([-2 2]);
            colormap('parula');
            colorbar;
            title(sprintf('Frame %d', idx));
            pause(0.15);
        end
    end

end