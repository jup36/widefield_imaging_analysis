
filePath = '/Volumes/buschman/Rodent Data/dualImaging_parkj/collectData/clusterW_output_CAmotifs_tcorr_mat_FitPhenoK_072125.mat';
saveFigDir = '/Volumes/buschman/Rodent Data/dualImaging_parkj/collectFigure/visualCluster'; 
if ispc
    filePath = 'Z:\Rodent Data\dualImaging_parkj\collectData\clusterW_output_CAmotifs_tcorr_mat_FitPhenoK_072125.mat';
    saveFigDir = 'Z:\Rodent Data\dualImaging_parkj\collectFigure\visualCluster';
end

load(fullfile(filePath), 'ovr_q', 'opts', 'num_clust')
kvals = opts.k_range;

fp = fig_params;
figure; hold on;

yyaxis left
plot(kvals, num_clust, 'o-', 'color', 'k', 'linewidth', 1.25);
yticks(0:50:max(num_clust)+50);  % set ticks every 50

yyaxis right
plot(kvals, ovr_q, 'o-', 'color', 'b', 'linewidth', 1.25);
ylabel('Modularity (Q)');

xlabel('Number of Neighbors');
title('Autofitting PhenoCluster KNN', ...
    'FontSize', fp.font_size, 'Fontweight', fp.font_weight);

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'b';
fp.FormatAxes(gca);

fh = gcf;

timestampStr = datestr(now, 'mmddyy_HHMMSS');  % Date and time string
figSaveName = sprintf('phenoK_modularity_over_kNN_CaMotifs_%s', timestampStr);
print(fullfile(saveFigDir, figSaveName),'-dpdf','-painters','-bestfit')