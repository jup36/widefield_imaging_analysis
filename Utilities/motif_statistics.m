 filePathM = '/Volumes/buschman/Rodent Data/Wide Field Microscopy/ExampleData/Preprocessed/m1092_100924_baseline_1_red_dff_combined_processed_fit_chunk1.mat'; 
 filePathPre = '/Volumes/buschman/Rodent Data/Wide Field Microscopy/ExampleData/Preprocessed/m1092_100924_baseline_1_red_dff_combined_processed.mat'; 
 savePathM = '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1092_jRGECO/m1092_100924/baseline/Figure'; 
 header = extract_date_animalID_header(filePathM); 
 load(filePathM, 'w', 'h', 'nanpxs', 'stats_train', 'stats_test')
 load(filePathPre, 'data_test')
 
% visualize the results of the stats_train
fprintf('\nNumber of Motifs: %d\n', stats_train.n_motifs)
fprintf('\nPercent Variance Explained: %0.2g%%\n', stats_train.pev*100)

%plot the relative pev of each motif
loadings = [stats_train.loadings{:}];
figure; bar(loadings); xlabel('motifs'); ylabel('PEV'); title('Motif Loadings')

%and the spatial correlation per frame
figure; plot(stats_train.rho_frame,'linewidth',2,'color','k'); 
title('spatial correlation')
ylabel('Correlation'); 
set(gca,'xticklabel',round(get(gca,'xtick')/10)); 
xlabel('time (s)'); 
%set(gca,'TickDir', 'out')
axis tight

%for comparison, let's see how well PCA captures the variance in this dataset
[~, score, ~, ~, explained, ~] = pca(data_test);
figure; hold on; plot(cumsum(explained),'color','k','linewidth',2); 
plot([0 100],[stats_train.pev*100,stats_train.pev*100],'linestyle','--','color','r','linewidth',2);
plot([stats_train.n_motifs,stats_train.n_motifs],[0 100],'linestyle','--','color','r','linewidth',2);
set(gca,'xlim',[0 100]); title('Comparing Motif fit to PCA')
xlabel('PCs'); 
ylabel('PEV'); 

%visualize some of the PCs as a gutcheck. Sometimes deconovlution can cause individual pixels to dominate variance. 
% If so, you may want to turn on the gp.w_pc_denoise option during processing  
temp = conditionDffMat(score',nanpxs);
%temp = SpatialGaussian(temp,[1 1],'sigma'); %smooth for visualization
tiledlayout(3,3);
for i = 1:9; nexttile; imagesc(temp(:,:,i)); title(sprintf('PC %d',i));  axis square; axis off; end

%let's visualize the motifs. 
figure;
for i = 1:stats_train.n_motifs
    temp = conditionDffMat(squeeze(w(:,i,:))',nanpxs);
    %smooth for visualization
    %temp = SpatialGaussian(temp,[1 1],'sigma');
    for j = 1:size(w,3)
        cvals = prctile(temp(:),99.9);
        imagesc(temp(:,:,j),[0 cvals]); title(sprintf('Motif %d Frame %d',i,j)) 
        pause(0.1); 
    end   
    pause(1); %in between motifs
end

%Plot the motif temporal weightings. You'll see that some motifs occur multiple times, others just once. 
% In general, spontaneous activity is relatively sparse
figure; hold on; plot(h'); legend; ylabel('Motif Activity'); set(gca,'xticklabel',round(get(gca,'xtick')/10)); xlabel('time (s)');
legend({'motif1', 'motif2', 'motif3', 'motif4', 'motif5', 'motif6', 'motif7'}, 'Location', 'best');

