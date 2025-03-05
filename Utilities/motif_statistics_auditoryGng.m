filePathChunks = findFilePattern('/Volumes/buschman/Rodent Data/Wide Field Microscopy/ExampleData/Preprocessed/',...
    'm1045_122424_task_day4-8_img_processed_green_dff_combined_fit_chunk*.mat'); % cell array containing paths to all chunks
filePathM = '/Volumes/buschman/Rodent Data/Wide Field Microscopy/ExampleData/Preprocessed/m1045_122424_task_day4-8_img_processed_green_dff_combined_fit_chunk1.mat';
filePathPre = '/Volumes/buschman/Rodent Data/Wide Field Microscopy/ExampleData/Preprocessed/m1045_122424_task_day4-8_img_processed_green_dff_combined.mat';
savePathM = '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1045_jRGECO_GRABda/m1045_122424/task/Figure';
header = extract_date_animalID_header(filePathM);
load(filePathM, 'w', 'h', 'nanpxs', 'stats_train', 'stats_test')
load(filePathPre, 'data_train', 'data_test')

%% Visualize number of motifs and PEV for each train and test splits
for f = 1:length(filePathChunks) 
    load(filePathChunks{f}, 'w', 'h', 'nanpxs', 'stats_train', 'stats_test'); 
    
    n_motifC_train{1, f} = stats_train.n_motifs; 
    pevC_train{1, f} = stats_train.pev; 

    n_motifC_test{1, f} = stats_test.n_motifs; 
    pevC_test{1, f} = stats_test.pev; 
end

% plot number of motifs
figure;
b = bar([cell2mat(n_motifC_train); cell2mat(n_motifC_test)]', 'grouped');

% Set colors for each group
b(1).FaceColor = 'b'; % Blue for train
b(2).FaceColor = 'r'; % green for test

% Labels and title
xlabel('Train/Test Splits');
ylabel('Motif Count');
title('Train vs. Test Motif Counts');
legend({'Train', 'Test'});
grid on;

% plot PEV
figure;
b = bar([cell2mat(pevC_train); cell2mat(pevC_test)]', 'grouped');

% Set colors for each group
b(1).FaceColor = 'b'; % Blue for train
b(2).FaceColor = 'r'; % green for test

% Labels and title
xlabel('Train/Test Splits');
ylabel('PEV');
title('Train vs. Test Motif PEV');
legend({'Train', 'Test'});
grid on;

%% visualize the results of the stats_train
fprintf('\nNumber of Motifs: %d\n', stats_test.n_motifs)
fprintf('\nPercent Variance Explained: %0.2g%%\n', stats_test.pev*100)

%%plot the relative pev of each motif
loadings = [stats_test.loadings{:}];
figure; bar(loadings); xlabel('motifs'); ylabel('PEV'); title('Motif Loadings')
%set(gca, 'TickDir', 'out')
% print(fullfile(savePathM, 'pevMotifs_heldout'), '-dpdf', '-bestfit', '-vector')

%and the spatial correlation per frame
figure; plot(stats_test.rho_frame,'linewidth',2,'color','k');
title('spatial correlation')
ylabel('Correlation');
set(gca,'xticklabel',round(get(gca,'xtick')/10));
xlabel('time (s)');
ylim([0 1]);
%set(gca,'TickDir', 'out')
axis tight
% print(fullfile(savePathM, 'spatialCorrMotifs_heldout'), '-dpdf', '-bestfit', '-vector')


%% for comparison, let's see how well PCA captures the variance in this dataset]
% Perform PCA on training data
[pcW_train, score_train, ~, ~, explained_train] = pca(data_train_chunk);

% Project test data onto the learned principal components
score_test = data_test_chunk * pcW_train; % Projecting test data using train PCs

% Compute total variance of test data
total_var_test = sum(var(data_test_chunk, 0, 1)); % Variance across pixels (dim 1)

% Compute variance captugreen by each PC in test data
var_test_PCs = var(score_test, 0, 1); % Variance along PC dimensions

% Recalculate explained variance for test data
explained_test = 100 * var_test_PCs / total_var_test; % Percentage variance explained

% Plot explained variance for test data
figure;
plot(cumsum(explained_test), '-o', 'LineWidth', 2);
xlabel('Number of Principal Components');
ylabel('Cumulative Explained Variance (%)');
title('PCA Explained Variance (Test Set)');
grid on;

figure; hold on; plot(cumsum(explained_test),'color','k','linewidth',2);
plot([0 100],[stats_test.pev*100,stats_test.pev*100],'linestyle','--','color','r','linewidth',2);
plot([stats_test.n_motifs,stats_test.n_motifs],[0 100],'linestyle','--','color','r','linewidth',2);
set(gca,'xlim',[0 100]); title('Comparing Motif fit to PCA')
xlabel('PCs');
ylabel('PEV');

%visualize some of the PCs as a gutcheck. Sometimes deconovlution can cause individual pixels to dominate variance.
% If so, you may want to turn on the gp.w_pc_denoise option during processing
temp = conditionDffMat(score_train', nanpxs);
%temp = SpatialGaussian(temp,[1 1],'sigma'); %smooth for visualization
tiledlayout(3,3);
for i = 1:9; nexttile; imagesc(temp(:,:,i)); title(sprintf('PC %d',i));  axis square; axis off; end

%let's visualize the motifs.
figure;
for i = 1:stats_test.n_motifs
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
set(gca, 'TickDir', 'out')
% print(fullfile(savePathM, 'motif_temporal_weighting'), '-dpdf', '-bestfit', '-vector')

