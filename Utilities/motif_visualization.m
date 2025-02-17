filePathM = '/Volumes/buschman/Rodent Data/Wide Field Microscopy/ExampleData/Preprocessed/m1092_100924_baseline_1_red_dff_combined_processed_fit_chunk1.mat'; 
filePathPre = '/Volumes/buschman/Rodent Data/Wide Field Microscopy/ExampleData/Preprocessed/m1092_100924_baseline_1_red_dff_combined_processed.mat'; 
savePathM = '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1092_jRGECO/m1092_100924/baseline/Figure'; 
header = extract_date_animalID_header(filePathM); 
load(filePathM, 'w', 'h', 'nanpxs', 'stats_train'); 
preS = load(filePathPre, 'data_train', 'data_test');

%% montage motifs and save each motif
for i = 1:stats_train.n_motifs
    motifs(:,:,:,i) = conditionDffMat(squeeze(w(:,i,:))',nanpxs);
    montage(motifs(:,:,:,i), 'Size', [1,10],'DisplayRange', [0 1]); colormap turbo; colorbar %change 4th dimension for each motif
    %print(fullfile(savePathM, [header, sprintf('_montage_motif%d', i)]), '-bestfit', '-dpdf', '-vector')
    fprintf("printed montage of motif#%d\n", i); 
end

%% tensor convolve for full reconstruction 
wh_conv = tensor_convolve(w,h);
wh = conditionDffMat(wh_conv', nanpxs); 
montage(wh, 'Size', [24,50],'DisplayRange', [0 0.4]); colormap turbo; colorbar %change 4th dimension for each motif
print(fullfile(savePathM, [header, '_montage_full_reconstruction']), '-bestfit', '-dpdf', '-vector')

%% manual visualization of motif expression comparing between the raw and the reconstructed
motifN = 4; 
data_train_rs = conditionDffMat(preS.data_train', nanpxs); % reshape train data back to 64x64
ht = h(motifN,:)'; 
ht(:,2) = 1:length(h); 
ht_sort = sortrows(ht, -1); 

hRank = 1; 
tp = ht_sort(hRank, 2); 
figure; montage(data_train_rs(:, :, tp-6:tp+3), 'Size', [1,10], 'DisplayRange', [0 0.3]); colormap turbo
title(sprintf("motif%d_timepoint%d_hRank%d_raw", motifN, tp, hRank), 'interpreter', 'none'); 
print(fullfile(savePathM, strcat(header, '_', sprintf("motif%d_timepoint%d_hRank%d_raw", motifN, tp, hRank))), '-dpdf', '-vector', '-bestfit')

figure; montage(wh(:, :, tp-6:tp+3), 'Size', [1,10], 'DisplayRange', [0 0.3]); colormap turbo
title(sprintf("motif%d_timepoint%d_hRank%d_reconstruct", motifN, tp, hRank), 'interpreter', 'none'); 
print(fullfile(savePathM, strcat(header, '_', sprintf("motif%d_timepoint%d_hRank%d_reconstruct", motifN, tp, hRank))), '-dpdf', '-vector', '-bestfit')