filePath_base = '/Volumes/buschman/Rodent Data/Wide Field Microscopy/ExampleData/Preprocessed'; 

filePath_W = {'m1045_121124_task_day4_img_processed_red_dff_combined_fit_chunk*.mat', ...
              'm1045_121324_task_day4-1_img_processed_red_dff_combined_fit_chunk*.mat', ...
              'm1045_121624_task_day4-2_img_processed_red_dff_combined_fit_chunk*.mat', ...
              'm1045_121824_task_day4-4_img_processed_red_dff_combined_fit_chunk*.mat', ...
              'm1045_122424_task_day4-8_img_processed_red_dff_combined_fit_chunk*.mat'}; 

% filePath_W = {'m1045_121124_task_day4_img_processed_green_dff_combined_fit_chunk*.mat', ...
%               'm1045_121324_task_day4-1_img_processed_green_dff_combined_fit_chunk*.mat', ...
%               'm1045_121624_task_day4-2_img_processed_green_dff_combined_fit_chunk*.mat', ...
%               'm1045_121824_task_day4-4_img_processed_green_dff_combined_fit_chunk*.mat', ...
%               'm1045_122424_task_day4-8_img_processed_green_dff_combined_fit_chunk*.mat'}; 

% filePath_W = {'m1045_121324_task_day4-1_img_processed_red_dff_combined_fit_chunk*.mat', ...
%               'm1045_121824_task_day4-4_img_processed_red_dff_combined_fit_chunk*.mat', ...
%               'm1045_122424_task_day4-8_img_processed_red_dff_combined_fit_chunk*.mat'}; 

parameter_class = 'general_params_dual';
gp = loadobj(feval(parameter_class));

%% Grab Ws
wC = {}; 
nanpxsC = {}; 

for f = 1:length(filePath_W)
    filePath_chunks = findFilePattern(filePath_base, filePath_W{f}); 
    for ff = 1:length(filePath_chunks)
        load(filePath_chunks{ff}, 'w', 'nanpxs')
        wC = [wC, {w}]; 
        nanpxsC = [nanpxsC, {nanpxs}]; 
        fprintf("Finished loading chunk #%d of file #%d!\n", ff, f); 
    end
end

% Compare each matrix to the first matrix
is_equal_nanpxsC = cellfun(@(x) isequal(x, nanpxsC{1}), nanpxsC);

% Check if all nan pixels are equal
assert(all(is_equal_nanpxsC));

[W_basis, kval, ovr_q, cluster_idx, idx_knn, tcorr_mat, handles, lag_mat, lags, nanpxs, shift] = ClusterW(wC, gp, nanpxsC{1}); 
% temporarily save ClusterW output
save(fullfile('/Volumes/buschman/Rodent Data/dualImaging_parkj/m1045_jRGECO_GRABda', 'clusterW_output_CAmotifs'), ...
    'W_basis', 'kval', 'ovr_q', 'cluster_idx', 'idx_knn', 'tcorr_mat', 'lag_mat', 'lags', 'nanpxs', 'shift');  


PlotMotifs
RefitBasisMotifs