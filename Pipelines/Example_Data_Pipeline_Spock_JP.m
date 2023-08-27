
% Open ssh connection
username = input(' Spock Username: ', 's');
password = passcode();
s_conn = ssh2_config('spock.princeton.edu',username,password);

%Add paths
addpath(genpath('/Volumes/buschman/Rodent Data/Wide Field Microscopy/fpCNMF'));
%addpath(genpath('/Volumes/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis'));

%configure preprocessing options
opts = ConfigurePreProcessing('crop_w',540,'vasc_std',2,'save_uncorrected',0);

%load general params (this is for anything after preprocessing)
parameter_class = 'general_params_example';
gp = loadobj(feval(parameter_class)); %this is not needed here, but demonstrates how I load this class in other functions by just passing the string. 

bdat = load(fullfile('/Volumes/buschman/Users/Caroline/NADA_dynamics/data/DA001_072623_behv_data.mat'), 'data'); 
bdat = bdat.('data'); 

%% Manual Portion
%select folders to process and grab the first file from each rec.
%EXAMPLE DATA: Select 'subfolders' and then select '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA001/DA001_072623/DA001_072623_img'
[file_list_first_stack,folder_list_raw] = GrabFiles_sort_trials('Pos0.ome.tif',1, ... % use GrabFiles_sort_trials to sort both files and folders 
    {'/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA001/DA001_072623/DA001_072623_img'},true);

%Grab reference images for each. Preload so no delay between loop.
ref_img = GetReferenceImage(file_list_first_stack{1},opts.fixed_image); % use the first frame of the first trial 

%manual allignment 
prepro_log = ManualAlignment(ref_img,opts);
% 1. crop (move the window and 2-click)
% 2. midline (draw a line with two dots and 2-click)
% 3. bregma (drop a dot)

%mask vasculature and manual cleanup (optional)
prepro_log = MaskVasculature(...
    prepro_log.cropped_alligned_img,prepro_log);
close; 
%no transformation
prepro_log.tform = []; 
prepro_log.output_size = [];

% Skipping 'RegisterReferenceImages.m', as we're not combining data across
% sessions, in which case registering to a common reference frame would be
% necessary. 

%save off the options to each folder
save([folder_list_raw{1} filesep 'prepro_log'],'prepro_log') 

%% Run PreProcess on Spock
file_list_preprocessed = cell(1,numel(folder_list_raw));

for cur_fold = 1:numel(folder_list_raw)
    [file_list_raw,~] = GrabFiles('.tif',0,folder_list_raw(cur_fold));
    [opts_list,~] = GrabFiles('prepro_log.m',0,folder_list_raw(1)); % use opts from the 1st trial 
    
    %Create spock bash script for each file and run it
    job_id = cell(1,numel(file_list_raw));
    for cur_file = 1:numel(file_list_raw)
        input_val = {ConvertMacToBucketPath(file_list_raw{cur_file}), ConvertMacToBucketPath(opts_list{1})};
        script_name = WriteBashScriptMac(sprintf('%d_%d',cur_fold,cur_file),'Spock_Preprocessing_Pipeline',input_val,{"'%s'","'%s'"},...
            'sbatch_time',15,'sbatch_memory',8);  %
        
        %Run job
        response = ssh2_command(s_conn,...
            ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory
            sprintf('sbatch %s',script_name)]); 
        
        %get job id
        job_id{cur_file} = erase(response.command_result{1},'Submitted batch job ');
        if cur_file ~=numel(file_list_raw)
            job_id{cur_file} = [job_id{cur_file} ','];
        end
    end    
end

%% Preprocessing. Results in a single hemo corrected, masked recording for each day in the 'preprocessed' folder
%gather recordings
cur_fold = 1; % trial #3, 5, 13 has normal signal. 
file_list_preprocessed = cell(1,numel(folder_list_raw));
[file_list_raw,~] = GrabFiles('.tif',0,folder_list_raw(cur_fold));
[opts_list,~] = GrabFiles('prepro_log.m',0,folder_list_raw(1)); % always point to the first folder 
opts = load(opts_list{1});
opts = opts.prepro_log;
[path, fn] = fileparts(file_list_raw{1});

%Process
stack = PreProcess(file_list_raw{1},opts);

%make Dff %NOTE: METHOD MAY MAKE TRIAL-based recs weird
%stack_b = stack(:,:,1:2:end);
%dff = makeDFF(stack_b, opts); 

%basic dff
stack_b = stack(:,:,1:2:end);  
avg_proj = nanmean(stack_b(:,:,1:15),3); % baseline average (0.5 sec)
dff = (double(stack_b)-avg_proj)./avg_proj; 

%examine the data
temp = reshape(dff,size(dff,1)*size(dff,2),size(dff,3));
imagesc(temp)

%histrogram
figure; hold on; 
histogram(dff(:))

%play our video; 
close all;
for i = 1:size(dff,3)
   imagesc(dff(:,:,i),[-0.03, 0.03]);
   title(sprintf('frame %d of %d',i,size(dff,3)));
   pause(0.1)
end

%if you want to see what the preprocessed data looks like then run
% InspectPreprocessedData(PreprocessedDataFilepath,'preprocessed')

