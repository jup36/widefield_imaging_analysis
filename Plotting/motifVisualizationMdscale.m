% Inputs:
%   tcorr_mat    – NxN similarity matrix based on cross-correlation
%   cluster_idxC – 1xN cell array of cluster labels (length K)
%% whereabouts 
filePath = '/Volumes/buschman/Rodent Data/dualImaging_parkj/collectData'; 
if ispc
   filePath = 'Z:\Rodent Data\dualImaging_parkj\collectData'; 
end
load(fullfile(filePath, 'clusterW_output_DAmotifs.mat'), 'tcorr_mat', 'cluster_idxC')

%% Perform MDS (get cooridnates)
Ymds = motifMdsY(tcorr_mat); 

%% Visualize and Print
for j = 1:numel(cluster_idxC)
    motifVisualizationMdsFunc(Ymds, cluster_idxC{1, j}, 1)
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Y = motifMdsY(tcorr_mat)
% Performs multi-dimensional scaling on the temporal correlation (similarity) matrix 

% Convert to dissimilarity matrix
D = 1 - tcorr_mat;
D(isnan(D)) = 1;                     % Treat NaN as max dissimilarity
D = min(max(D, 0), 2);               % Clip to [0, 2]
D(1:size(D,1)+1:end) = 0;            % Force diagonal to zero
D = (D + D') / 2;                    % Ensure symmetry

% Classical MDS
Y = mdscale(D, 2, 'Criterion', 'metricstress');
end

function motifVisualizationMdsFunc(Y, cluster_idx, figSaveLogic)
% Visualize motifs using MDS (multidimensional scaling) based on temporal correlation
% Inputs:
%   tcorr_mat   – [N × N] similarity matrix (e.g., from temporal cross-correlation)
%   cluster_idx – [N × 1] cluster assignments

saveFigDir = '/Volumes/buschman/Rodent Data/dualImaging_parkj/collectFigure/visualCluster';
if ispc
    saveFigDir = 'Z:\Rodent Data\dualImaging_parkj\collectFigure\visualCluster';
end

% Cluster info
cluster_ids = unique(cluster_idx);
K = numel(cluster_ids);
nMotifs = size(cluster_idx, 1); 

% Compute centroids in MDS space
centroids = zeros(K, 2);
for c = 1:K
    centroids(c, :) = mean(Y(cluster_idx == cluster_ids(c), :), 1);
end

% Use PCA to determine hue ordering
coeff = pca(centroids);
centroid_proj = centroids * coeff(:, 1:2);  % project to 2D
centroid_proj = normalize(centroid_proj, 'range');
hues = centroid_proj(:, 1);  % use 1st PC for hue
sat = 0.6; val = 0.95;
hsv_colors = hsv2rgb([hues, repmat(sat, K, 1), repmat(val, K, 1)]);

% Map each cluster to a color
cluster_color_map = containers.Map(num2cell(cluster_ids), num2cell(hsv_colors, 2));
motif_colors = zeros(nMotifs, 3);
for j = 1:nMotifs
    motif_colors(j, :) = cluster_color_map(cluster_idx(j));
end

% Plot
figure;
scatter(Y(:,1), Y(:,2), 10, motif_colors, 'filled'); hold on;
title(sprintf('MDS: %d Clusters', K));
axis equal off;

% Annotate cluster centers
for c = 1:K
    text(centroids(c,1), centroids(c,2), sprintf('c%d', cluster_ids(c)), ...
        'Color', 'k', 'FontSize', 10, 'HorizontalAlignment', 'center');
end

%sgtitle('MDS of Motif Similarity Based on Cross-Correlation');

if figSaveLogic
    timestampStr = datestr(now, 'mmddyy_HHMMSS');  % Date and time string
    figSaveName = sprintf('clusterMDS_total%dmotifs_%s', K, timestampStr);
    print(fullfile(saveFigDir, figSaveName),'-dpdf','-painters','-bestfit')
end


end
