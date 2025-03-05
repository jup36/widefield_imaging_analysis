function [swarm_id,save_fn] = FitMotifs_SpockSwarm_chunks(fn,dependency_id,s_conn,parameter_class,nChunks)
%Camden MacDowell
%once file fn is done (dependency), parallellizes the fitting of each data
%chunk in fn

gp = loadobj(feval(parameter_class)); 

%generate swarm
swarm_id = cell(1,nChunks);
save_fn = cell(1,nChunks);
for i = 1:nChunks
    [~, fn_temp] = fileparts(fn);  
    save_fn{i} = [gp.local_bucket_mac gp.processing_intermediates_mac fn_temp sprintf('_fit_chunk%d.mat',i)];
        
    script_name = WriteBashScriptMac(sprintf('motifchunk%d',i), ...
        'FitMotifs_Spock',{ConvertMacToBucketPath(fn), ...
        ConvertMacToBucketPath(save_fn{i})...
        , i, parameter_class},{"'%s'","'%s'","%d","'%s'"},...
        'sbatch_time',600,'sbatch_memory',12,...
        'sbatch_path',"/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/Spock/"); %cutoff for spock priority is 4hrs (240s), then 48 hours 
    
    response = ssh2_command(s_conn,...
        ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory
        sprintf('sbatch --dependency=afterok:%s %s',dependency_id,script_name)]);
    
    swarm_id{i} = erase(response.command_result{1},'Submitted batch job ');
    if i ~=nChunks
        swarm_id{i} = [swarm_id{i} ','];
    end
end