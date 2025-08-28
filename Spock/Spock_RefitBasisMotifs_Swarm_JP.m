function Spock_RefitBasisMotifs_Swarm_JP(filePathImg, fileKeyword, basis_dir, parameter_class, save_dir)
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
% basis_dir = '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1045_jRGECO_GRABda/clusterW_output_CAmotifs.mat';
% parameter_class = 'general_params_dual';
% save_dir = '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1045_jRGECO_GRABda/m1045_122424/task/Matfiles';

if exist('s_conn', 'var')~=1
    % Open ssh connection
    username = input(' Spock Username: ', 's'); % Use PU (NOT PNI) id and password as of 11/1/24
    password = passcode();
    s_conn = ssh2_config('spock.princeton.edu',username,password);
end

if iscell(filePathImg)
    filePathImg = filePathImg{1};  
end

%Add paths
if ispc
    addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\fpCNMF'));
    addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\Widefield_Imaging_Analysis'));
else
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/fpCNMF'));
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis'));    
end

% get 'file_processed' to access the train and test datasets used for motif fitting
gp = loadobj(feval(parameter_class));
[~, fileheader] = fileparts(filePathImg);
file_processed = [gp.local_bucket_mac gp.processing_intermediates_mac fileheader '_processed' fileKeyword]; % path in temporary folder with the train/test split data

%generate swarm
temp = load(file_processed, 'data_test');
for cur_chunk = 1:size(temp.data_test,3)
    script_name = WriteBashScriptMac(sprintf('refitchunk%d',cur_chunk),'RefitBasisMotifs_JP',{ConvertMacToBucketPath(file_processed),...
        ConvertMacToBucketPath(basis_dir),cur_chunk,parameter_class,ConvertMacToBucketPath(save_dir)},...
        {"'%s'","'%s'",'%d',"'%s'","'%s'"},...
        'sbatch_time',59,'sbatch_memory',10,...
        'sbatch_path',"/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/Spock/");
    ssh2_command(s_conn,...
        ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory
        sprintf('sbatch %s',script_name)]);

    % if ~isempty(job_id)
    %     ssh2_command(s_conn,...
    %         ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory
    %         sprintf('sbatch --dependency=afterok:%s %s',job_id,script_name)]);
    % else
    %     ssh2_command(s_conn,...
    %         ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory
    %         sprintf('sbatch %s',script_name)]);
    % end
end

