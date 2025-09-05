function Scotty_RefitBasisMotifs_Swarm_JP(filePathImg, fileKeyword, basis_dir, parameter_class, save_dir)
%Camden MacDowell
% customized by Junchol Park

% Summary: This function commands to run 'RefitBasisMotifs.m' on Spock.
%   Input: 
%       - The train and test datasets ('data_train', 'data_test') generated for motif fitting 
%         is necessary; see below how the relevant file path is reconstructed using 'filePathImg', 'fileKeyword' and 'parameter_class' inputs.  
%       - basis_dir: directory where 'clusterW_output' is saved. 
%       - save_dir: directory to save refit results.  
%       
% E.g., 
% filePathImg = 
%       {'/Volumes/buschman/Rodent
%       Data/dualImaging_parkj/m1045_jRGECO_GRABda/m1045_122424/task/m1045_122424_task_day4-8_img'}; % usually a cell 
% fileKeyword = '_red_dff_combined.mat'; 
% basis_dir = '/Volumes/buschman/Rodent Data/dualImaging_parkj/collectData/clusterW_output_CAmotifs_072125.mat';
% parameter_class = 'general_params_dual';
% save_dir = '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1045_jRGECO_GRABda/m1045_122424/task/Matfiles';

% First time per session connect scotty
s_conn = ssh2_command_scotty('connect','scotty');   % no key path needed

if iscell(filePathImg)
    filePathImg = filePathImg{1}; 
end

filePathImg = compatiblepath(filePathImg); 
basis_dir = compatiblepath(basis_dir); 

% get params
gp = loadobj(feval(parameter_class));

% get 'file_processed' to access the train and test datasets used for motif fitting
[~, fileheader] = fileparts(filePathImg); 

file_processed_folder = find_keyword_containing_folder([gp.local_bucket gp.processing_intermediates], fileheader, 'recursive', false); 
if numel(file_processed_folder)>1
    warning("More than one processed folder was found!")
end

if iscell(file_processed_folder)
    file_processed_folder = file_processed_folder{1}; 
end

file_processed = fullfile(file_processed_folder, [fileheader, '_processed' fileKeyword]); % fileKeyword = '_red_dff_combined.mat'

if ~exist(file_processed, "file")
    error("Preprocessed data not detected!")
end

%% Main
temp = load(file_processed, 'data_test');
for cur_chunk = 1:size(temp.data_test,3)
    script_name = WriteBashScriptWinScotty(sprintf('refitchunk%d',cur_chunk),'RefitBasisMotifs_JP',{ConvertWinToBucketPath(file_processed),...
        ConvertWinToBucketPath(basis_dir),cur_chunk,parameter_class,ConvertWinToBucketPath(save_dir)},...
        {"'%s'","'%s'",'%d',"'%s'","'%s'"},...
        'sbatch_time',59,'sbatch_memory',10,...
        'sbatch_path',"/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/Spock/");

    % Submit a job
    remoteCmd = ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts ; ' ...
    sprintf('sbatch %s', script_name)];

    resp = ssh2_command_scotty(s_conn, remoteCmd);
    disp(resp.command_result{end});
end