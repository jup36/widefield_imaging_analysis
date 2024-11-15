addpath("Z:\Users\Manuel\DEMO")
%upload dff and mask
t_TEMP = t_TEMP_resize{2}; %change dimension for each registration
for i = 1:size(t_TEMP,3)
    Mask(:,:,i) = Mask(:,:,1);
end

for i = 1:size(t_TEMP,3)
    t_TEMP(~Mask)= NaN; 
end %nan pixels out of mask

dff = t_TEMP;
gp = 'genparams'; %parameter class handle to pass through functions

[~, nanpxs, data_train, data_test] = ProcessAndSplitData(dff,[],gp); %patience.... takes a few min
%I added a bunch of normalization options. Haven't explored all of these in great detail, 
% but wanted to added them in case as different data structures may warrant different methods of normalization. 

%Visualize the data.  
data_train_stack = conditionDffMat(data_train',nanpxs);

%for reproducibility
rng('default');
opts = loadobj(feval(gp));
lambda = 0.0077;


%Fit Motifs To Training Data And Collect Statistics
W_temp = cell(1,opts.repeat_fits);
H_temp = cell(1,opts.repeat_fits);
stats_train_temp =cell(1,opts.repeat_fits);
for cur_fit = 1:opts.repeat_fits %fit multiple times due to random initialization. Each iteration takes ~5min
    if opts.verbose; fprintf('\nFitting Training Data Round %d of %d',cur_fit,opts.repeat_fits); end
    [W_temp{cur_fit},H_temp{cur_fit},stats_train_temp{cur_fit}] = fpCNMF(data_train,'L',opts.L,'K',opts.K,'non_penalized_iter',...
        opts.non_penalized_iter,'penalized_iter',opts.penalized_iter,...
        'speed','fast','verbose',opts.verbose,'lambda',lambda,...
        'ortho_H',opts.ortho_H,'w_update_iter',opts.w_update_iter,...
        'sparse_H',opts.sparse_H);  
    %Remove Empty Motifs 
    [W_temp{cur_fit},H_temp{cur_fit}] = RemoveEmptyMotifs(W_temp{cur_fit},H_temp{cur_fit});
end

%choose best fit
idx = InternallyValidateWs(data_train,W_temp,H_temp,'AIC',1);
W = W_temp{idx}; H = H_temp{idx}; stats_train = stats_train_temp{idx}; 

%reorder by decreasing explained variance
[~,idx] = sort(stats_train.loadings{:},'descend');
W = W(:,idx,:); H = H(idx,:);
stats_train.lambda = lambda;


%montage motifs and save
for i = 1:stats_train.n_motifs
    motifs(:,:,:,i) = conditionDffMat(squeeze(W(:,i,:))',nanpxs);
end

montage(motifs(:,:,:,7), 'Size', [1,10],'DisplayRange', [0 1]); colormap turbo %change 4th dimension for each motif



