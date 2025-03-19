function refitBasisMotifsJP_local(filePathImg, fileKeyword, basis_dir, parameter_class, save_dir)
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
filePathImg = {'/Volumes/buschman/Rodent Data/dualImaging_parkj/m1045_jRGECO_GRABda/m1045_122424/task/m1045_122424_task_day4-8_img'}; % usually a cell 
fileKeyword = '_red_dff_combined.mat'; 
basis_dir = '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1045_jRGECO_GRABda/clusterW_output_CAmotifs.mat';
parameter_class = 'general_params_dual';
save_dir = '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1045_jRGECO_GRABda/m1045_122424/task/Matfiles';

% get 'file_processed' to access the train and test datasets used for motif fitting
gp = loadobj(feval(parameter_class));
[~, fileheader] = fileparts(filePathImg{1});
file_processed = [gp.local_bucket_mac gp.processing_intermediates_mac fileheader '_processed' fileKeyword]; % path in temporary folder with the train/test split data

%load the basis motifs
tempW = load(basis_dir,'W_basis');
%temp = load(basis_dir,'W_basis','noise_clusters');
%ignore noise clusters when refitting (you could also refit these and then remove them and their representation in the data
W_basis = tempW.W_basis;
%W(:,temp.noise_clusters,:)=[];

%load the data_train and the data_test data
fprintf('\n\tLoadings Data\n'); 
temp = load(file_processed,'data_test','data_train','nanpxs');
chunkN = size(temp.data_test, 3); 

for chunk = 1:chunkN

    data_train = temp.data_train(:,:,chunk);
    data_test = temp.data_test(:,:,chunk);
    nanpxs = temp.nanpxs;
    [~,name] = fileparts(file_processed);
    name = [name, sprintf('_refitchunk_%d',chunk)];

    %recondition, smooth, and flatten
    if numel(gp.smt_kernel)>2 %autodetermine appropriate smoothing value
        fprintf('\n\t Autofitting Smoothing Value');
        gp.smt_kernel = AutoFitSmoothingLevel(cat(2,data_train,data_test),nanpxs,W,gp);
        data_train = reshape(SpatialGaussian(conditionDffMat(data_train',nanpxs),gp.smt_kernel),gp.pixel_dim(1)*gp.pixel_dim(2),size(data_train,2));
        data_test = reshape(SpatialGaussian(conditionDffMat(data_test',nanpxs),gp.smt_kernel),gp.pixel_dim(1)*gp.pixel_dim(2),size(data_test,2));
    elseif numel(gp.smt_kernel)==2
        fprintf('\n\t Set Smoothing Value');
        data_train = reshape(SpatialGaussian(conditionDffMat(data_train',nanpxs),gp.smt_kernel),gp.pixel_dim(1)*gp.pixel_dim(2),size(data_train,2));
        data_test = reshape(SpatialGaussian(conditionDffMat(data_test',nanpxs),gp.smt_kernel),gp.pixel_dim(1)*gp.pixel_dim(2),size(data_test,2));
    else
        fprintf('\n\t No Smoothing');
        data_train = reshape(conditionDffMat(data_train',nanpxs),gp.pixel_dim(1)*gp.pixel_dim(2),size(data_train,2));
        data_test = reshape(conditionDffMat(data_test',nanpxs),gp.pixel_dim(1)*gp.pixel_dim(2),size(data_test,2));
    end

    %Only work on shared pixels
    [~, bad_pxl] = SharedPixels(W_basis,cat(2,data_test,data_train));
    data_train(bad_pxl,:) = [];
    data_test(bad_pxl,:) = [];
    W = W_basis; 
    W(bad_pxl,:,:) = [];

    %fit motif weightings to data
    [w,H] = fpCNMF(data_train,'non_penalized_iter',...
        gp.non_penalized_iter,'penalized_iter',gp.penalized_iter_refit,...
        'speed','fast','verbose',0,'lambda',0,'w_update_iter',0,...
        'ortho_H',gp.ortho_H,'W',W,'sparse_H',0);
    stats_refit = CNMF_Stats(w,H,data_train,0);
    stats_refit.smoothingkernel = gp.smt_kernel;

    %get residuals
    residuals = data_train-tensor_convolve(w,H);

    %save off
    save([save_dir filesep name 'train.mat'],'w','H','stats_refit','bad_pxl','residuals');
    fprintf("Completed refitting train chunk #%d/%d\n", chunk, chunkN); 

    [w,H] = fpCNMF(data_test,'non_penalized_iter',...
        gp.non_penalized_iter,'penalized_iter',gp.penalized_iter_refit,...
        'speed','fast','verbose',0,'lambda',0,'w_update_iter',0,...
        'ortho_H',gp.ortho_H,'W',W,'sparse_H',0);
    stats_refit = CNMF_Stats(w,H,data_test,0);
    stats_refit.smoothingkernel = gp.smt_kernel;

    %get residuals
    residuals = data_test-tensor_convolve(w,H);

    save([save_dir filesep name 'test.mat'],'w','H','stats_refit','bad_pxl','residuals');
    fprintf("Completed refitting test chunk #%d/%d\n", chunk, chunkN); 
end