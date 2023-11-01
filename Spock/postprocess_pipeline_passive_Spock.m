function postprocess_pipeline_passive_Spock(filePath)
%filePath = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA003/DA003_083023/DA003_083023_img';
% Start from an img folder where ome tiff stacks are saved per folder. 
% Preprocessing requirements: 
%   1) parse_tbyt_passive_behavior.m (output: '*tbytDat.mat', the trial-by-trial data structure.)
%   2) Data_pipeline_passive_Spock.m 
%       2-1) Spock_Preprocessing_Pipeline.m (output: '*ome_stack.mat', crop, align, spatial binning etc.)
%       2-2) Spock_CombineStacksBVcorrect.m (output: '*dff_combined.mat', combine stacks (optional),
%       hemodynamic correction)
% 

%% Open ssh connection
username = input(' Spock Username: ', 's');
password = passcode();
s_conn = ssh2_config('spock.princeton.edu',username,password);

%% load tbytDat
parentDir = fileparts(filePath);
[~, header] = fileparts(parentDir);
[fileBeh, folderBeh] = GrabFiles_sort_trials('tbytDat', 1, {parentDir});
load(fullfile(fileBeh{1}), 'tbytDat'); 

%% file directory for trials (use this folder to save dff, frames.png, and video.mp4)
filePathTrials = fullfile(parentDir, strcat(header, '_trials'));
if exist(filePathTrials, 'dir') == 0
    mkdir(filePathTrials);
end

%% load preprocessed dffs
[file_list_dff, folder_list_dff] = GrabFiles_sort_trials('dff_combined.mat', 1, {filePath}); % use GrabFiles_sort_trials to sort both files and folders
dffIC = cellfun(@(a) regexp(a, '_(\d{1,3})_dff', 'tokens', 'once'), file_list_dff, 'un', 0); 
dffI = cell2mat(cellfun(@(a) str2double(a{1}), dffIC, 'un', 0)); % dff file index 

%% process tbytDat further and save it
for tt = 1:length(tbytDat)
    if ~isempty(tbytDat(tt).cmosExp)
        % get path for the corresponding dff
        dffPath = file_list_dff{tbytDat(tt).cmosExpTrainI==dffI};     
        
        % organize frame info so that relevant frames can be indexed within the dff stack 
        tbytDat(tt).frameT = tbytDat(tt).cmosExp(1:2:end); % frame time (needs to alternate due to interleaved violet frames for hemodynamic correction)
        temp1stFrameI = floor(tbytDat(tt).cmosExpPulsesOfTrain{1}./2); % 1st frame to take in dff (indices must be divided by 2 since dff already taken excluding violet frames) 
        tbytDat(tt).frameI = temp1stFrameI:temp1stFrameI+length(tbytDat(tt).frameT)-1; 

        % map temporal events to cmosExp pulses
        tbytDat(tt).frameLickI = check_timestamp_overlap(tbytDat(tt).frameT, tbytDat(tt).licks);
        tbytDat(tt).frameWaterI = check_timestamp_overlap(tbytDat(tt).frameT, tbytDat(tt).water);

        % timestamp each frame relative to the stim onset
        if tbytDat(tt).evtType < 3  % 1 or 2: visual (common, uncommon)
            tbytDat(tt).frameStimOnI = check_timestamp_overlap(tbytDat(tt).frameT, tbytDat(tt).evtOn);
            tbytDat(tt).frameStimOffI = check_timestamp_overlap(tbytDat(tt).frameT, tbytDat(tt).evtOff);
            tbytDat(tt).frameStimI = zeros(length(tbytDat(tt).frameT), 1);
            tbytDat(tt).frameStimI(find(tbytDat(tt).frameStimOnI, 1):find(tbytDat(tt).frameStimOffI, 1), 1) = 1;
        end
        tbytDat(tt).frameTrel = tbytDat(tt).frameT-tbytDat(tt).evtOn;
        
        % store relevant paths 
        tbytDat(tt).dffPath = file_list_dff{tbytDat(tt).cmosExpTrainI};
        tbytDat(tt).dffPathSpock = ConvertMacToBucketPath(dffPath);
        tbytDat(tt).dffSavePath = fullfile(filePathTrials, sprintf('block_%d_trial_%d', tbytDat(tt).cmosExpTrainI, tt));
        tbytDat(tt).dffSavePathSpock = ConvertMacToBucketPath(tbytDat(tt).dffSavePath); 
    end
end
clearvars tt
save(fullfile(folderBeh{1}, [header, '_tbytDat_dff.mat']), 'tbytDat')

% Create spock bash script for each file and run it
input_val = {ConvertMacToBucketPath(fileBeh{1})};
script_id = string(datetime('now', 'Format', 'MM-dd-yy-HH-mm'));
script_name = WriteBashScriptMac(script_id, 'dffPostprocess_Spock', input_val, {"'%s'"}, ...
    'sbatch_time',15,'sbatch_memory',8);  % bash script to run for preprocessing

% Run job
response = ssh2_command(s_conn,...
    ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory
    sprintf('sbatch %s',script_name)]);



job_id = cell(1, length(tbytDat));
for tt = 1:length(tbytDat)
    if ~isempty(tbytDat(tt).cmosExp)

        %Create spock bash script for each file and run it
        input_val = {ConvertMacToBucketPath(fileBeh{1})};
        script_id = string(datetime('now', 'Format', 'MM-dd-yy-HH-mm'));
        script_name = WriteBashScriptMac(script_id, 'dffPostprocess_Spock', input_val, {"'%s'"}, ...
            'sbatch_time',15,'sbatch_memory',8);  % bash script to run for preprocessin

        %Run job
        response = ssh2_command(s_conn,...
            ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory
            sprintf('sbatch %s',script_name)]);

        %get job id
        job_id{tt} = erase(response.command_result{1},'Submitted batch job ');



    end
end
clearvars tt







%Run job
response = ssh2_command(s_conn,...
    ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory
    sprintf('sbatch %s',script_name)]);

%get job id
job_id = erase(response.command_result{1},'Submitted batch job ');

% take corresponding frames with for each trial with 2D gaussian filtering
for tt = 1:length(tbytDat)
    if ~isempty(tbytDat(tt).cmosExp)
        dff = dffC{1, tbytDat(tt).cmosExpTrainI};

        tbytDat(tt).dff = dff(:,:,tbytDat(tt).frameI); % aligned dff
        tbytDat(tt).dffsm = applyImgaussfilt(tbytDat(tt).dff);


    end
    fprintf('processed dff of trial#%d\n', tt)
end




for cur_fold = 1:numel(folder_list_raw)
    [file_list_raw,~] = GrabFiles('.tif',0,folder_list_raw(cur_fold)); % note that there's only one file per folder in this experiment 
    [opts_list,~] = GrabFiles('prepro_log.m',0,folder_list_raw(1)); % use opts from the 1st trial

    %Create spock bash script for each file and run it
    job_id = cell(1,numel(file_list_raw));
    % for cur_file = 1:numel(file_list_raw)
    input_val = {ConvertMacToBucketPath(file_list_raw{1}), ConvertMacToBucketPath(opts_list{1})};
    script_name = WriteBashScriptMac(sprintf('%d_%d', cur_fold, 1),'Spock_Preprocessing_Pipeline',input_val,{"'%s'","'%s'"},...
        'sbatch_time',15,'sbatch_memory',8);  % bash script to run for preprocessing

    %Run job
    response = ssh2_command(s_conn,...
        ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory
        sprintf('sbatch %s',script_name)]);

    %get job id
    job_id{cur_fold} = erase(response.command_result{1},'Submitted batch job ');

    %Once each folder is done, combine all the stacks and do hemocorrection
    [~,header] = fileparts(ConvertMacToBucketPath(folder_list_raw{cur_fold}));
    file_list_preprocessed{cur_fold} = [folder_list_raw{cur_fold} filesep header '_dff_combined.mat'];
    script_name = WriteBashScriptMac(sprintf('%d_%d_combine', cur_fold, 1), ...
        'Spock_CombineStacksBVcorrect',{ConvertMacToBucketPath(folder_list_raw{cur_fold}), ConvertMacToBucketPath(file_list_preprocessed{cur_fold}), 'general_params_example'},{"'%s'","'%s'","'%s'"});

    % Run job with dependency
    response = ssh2_command(s_conn,...
        ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory
        sprintf('sbatch --dependency=afterok:%s %s',[job_id{:}],script_name)]);    

end









function record_dff_frames_Spock(tbytDat, filePathPng)
        % e.g., dff = tbytDat(21).dffsm
        %% record frames first
        close all;
        for t = 1:length(tbytDat) % trial
            if ~isempty(tbytDat(t).dffsm)
                pngSubDir = fullfile(filePathPng, sprintf('block_%d_trial_%d', tbytDat(t).cmosExpTrainI, t));

                if exist(pngSubDir, 'dir') == 7

                else
                    mkdir(pngSubDir);
                    for i = 1:size(tbytDat(t).dffsm,3)
                        pngName = sprintf('frame_%d.png', i);
                        figHandle = imageFrameWithNaNs(tbytDat(t).dffsm(:, :, i), [-2 2]); hold on;

                        % Turn off the axes
                        axis off;
                        % Set the figure's background to white
                        set(gcf, 'Color', 'w');

                        if ~isempty(tbytDat(t).frameStimI) && tbytDat(t).frameStimI(i) % for visual trials
                            if tbytDat(t).evtType == 1 % common visual stim (draw 45 degree lines at the upper left corner)
                                insertgrating45(figHandle, tbytDat(t).dffsm(:, :, i))
                            elseif tbytDat(t).evtType == 2 % common visual stim (draw 135 degree lines at the upper left corner)
                                insertgrating135(figHandle, tbytDat(t).dffsm(:, :, i))
                            end
                            hold off;
                        end

                        % Capture the current figure with dots
                        frameLabled = getframe(gca);
                        frameLabled = frameLabled.cdata;

                        % save the figure with dots

                        imwrite(frameLabled, fullfile(pngSubDir, pngName));

                        fprintf('Frame #%d of trial #%d is labeled and saved.\n', i, t);
                        close all
                    end
                end
            end
        end
end


end