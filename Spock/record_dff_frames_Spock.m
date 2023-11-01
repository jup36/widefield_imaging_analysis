
% Open ssh connection
username = input(' Spock Username: ', 's');
password = passcode();
s_conn = ssh2_config('spock.princeton.edu',username,password);








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