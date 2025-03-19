filePathM = '/Volumes/buschman/Rodent Data/Wide Field Microscopy/ExampleData/Preprocessed/m1045_122424_task_day4-8_img_processed_red_dff_combined_fit_chunk5.mat'; 
filePathPre = '/Volumes/buschman/Rodent Data/Wide Field Microscopy/ExampleData/Preprocessed/m1045_122424_task_day4-8_img_processed_red_dff_combined.mat'; 
savePathM = '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1045_jRGECO_GRABda/m1045_122424/task/Figure'; 
header = extract_date_animalID_header(filePathM); 
load(filePathM, 'w', 'h', 'nanpxs', 'stats_train'); 
preS = load(filePathPre, 'data_train', 'data_test');

tokens = regexp(filePathM, 'chunk(\d{1,3})', 'tokens');
chunkN = str2double(tokens{1}{1});

data_train_chunk = squeeze(preS.data_train(:, :, chunkN)); 
data_test_chunk = squeeze(preS.data_test(:, :, chunkN)); 

%% montage motifs and save each motif
for i = 1:stats_train.n_motifs
    motifs(:,:,:,i) = conditionDffMat(squeeze(w(:,i,:))',nanpxs);
    figure; montage(motifs(:,:,:,i), 'Size', [1,10],'DisplayRange', [0 0.5]); colormap turbo; colorbar %change 4th dimension for each motif
    %print(fullfile(savePathM, [header, sprintf('_montage_motif%d', i)]), '-bestfit', '-dpdf', '-vector')
    fprintf("printed montage of motif#%d\n", i); 
end

%% Create the figure and visualize all motifs in one montage
% Preallocate a 4D array to hold all motifs
allMotifs = zeros(64, 64, 1, 10 * stats_train.n_motifs); % 1 channel for grayscale

% Fill the allMotifs array
for i = 1:stats_train.n_motifs
    motifs = conditionDffMat(squeeze(w(:,i,:))', nanpxs); % Get 64x64x10 motif
    allMotifs(:,:,:, (i-1)*10 + (1:10)) = motifs; % Stack frames in correct order
end

% Create the figure and visualize all motifs in one montage
figure;
montage(allMotifs, 'Size', [stats_train.n_motifs, 10], 'DisplayRange', [0 0.3]); 
colormap turbo;
colorbar;
title('All Motifs');
print(fullfile(savePathM, [header, '_montage_allMotifs_', sprintf('chunk#%d', chunkN)]), '-bestfit', '-dpdf', '-vector')

%% tensor convolve for full reconstruction 
wh_conv = tensor_convolve(w,h);
wh = conditionDffMat(wh_conv', nanpxs); 
figure; montage(wh, 'Size', [24,50],'DisplayRange', [0 0.2]); colormap turbo; colorbar %change 4th dimension for each motif
print(fullfile(savePathM, [header, '_montage_full_reconstruction']), '-bestfit', '-dpdf', '-vector')

%% manual visualization of motif expression comparing between the raw and the reconstructed
motifN = 5; 
data_train_rs = conditionDffMat(data_train_chunk', nanpxs); % reshape train data back to 64x64
ht = h(motifN,:)'; 
ht(:,2) = 1:length(h); 
ht_sort = sortrows(ht, -1); 

hRank = 1; 
tp = ht_sort(hRank, 2); 
figure; montage(data_train_rs(:, :, tp-6:tp+3), 'Size', [1,10], 'DisplayRange', [0 0.3]); colormap turbo
title(sprintf("motif%d_timepoint%d_hRank%d_raw", motifN, tp, hRank), 'interpreter', 'none'); 
%print(fullfile(savePathM, strcat(header, '_', sprintf("motif%d_timepoint%d_hRank%d_raw", motifN, tp, hRank))), '-dpdf', '-vector', '-bestfit')

figure; montage(wh(:, :, tp-6:tp+3), 'Size', [1,10], 'DisplayRange', [0 0.3]); colormap turbo
title(sprintf("motif%d_timepoint%d_hRank%d_reconstruct", motifN, tp, hRank), 'interpreter', 'none'); 
%print(fullfile(savePathM, strcat(header, '_', sprintf("motif%d_timepoint%d_hRank%d_reconstruct", motifN, tp, hRank))), '-dpdf', '-vector', '-bestfit')