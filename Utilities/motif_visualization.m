 filePathM = '/Volumes/buschman/Rodent Data/Wide Field Microscopy/ExampleData/Preprocessed/m1092_100924_baseline_1_red_dff_combined_processed_fit_chunk1.mat'; 
 savePathM = '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1092_jRGECO/m1092_100924/baseline/Figure'; 
 header = extract_date_animalID_header(filePathM); 
 load(filePathM, 'w', 'h', 'nanpxs', 'stats_train')

%% montage motifs and save each motif
for i = 1:stats_train.n_motifs
    motifs(:,:,:,i) = conditionDffMat(squeeze(w(:,i,:))',nanpxs);
    montage(motifs(:,:,:,i), 'Size', [1,10],'DisplayRange', [0 1]); colormap turbo; colorbar %change 4th dimension for each motif
    print(fullfile(savePathM, [header, sprintf('_montage_motif%d', i)]), '-bestfit', '-dpdf', '-vector')
    fprintf("printed montage of motif#%d\n", i); 
end

%% tensor convolve
wh_conv = tensor_convolve(w,h);
wh = conditionDffMat(wh_conv', nanpxs); 
montage(wh, 'Size', [24,50],'DisplayRange', [0 0.4]); colormap turbo; colorbar %change 4th dimension for each motif
print(fullfile(savePathM, [header, '_montage_full_reconstruction']), '-bestfit', '-dpdf', '-vector')
