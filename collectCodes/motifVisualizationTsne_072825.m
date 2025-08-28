% Inputs:
%   W_basisC     – 1xN cell array of [P × K × L] motifs
%   cluster_idxC – 1xN cell array of cluster labels (length K)
%% whereabouts 
filePath = '/Volumes/buschman/Rodent Data/dualImaging_parkj/collectData'; 
if ispc
   filePath = 'Z:\Rodent Data\dualImaging_parkj\collectData'; 
end
load(fullfile(filePath, 'clusterW_output_DAmotifs_072125.mat'), 'W_aligned', 'cluster_idxC', 'nanpxs')

%% get Ws (remove NaN pixels and padded frames if applicable)
if size(W_aligned,1)==64*64 % this is the expected condition
    valPxs = ~ismember(1:size(W_aligned,1), nanpxs); 
else
    valPxs = true(1, size(W_aligned,1)); 
end

if size(W_aligned,3)==30 % this is the expected condition with 10 padded frames before and after the relevant frames
   valLags = 11:20; 
else
   valLags = 1:size(W_aligned,3); 
end

Ws = W_aligned(valPxs, :, valLags); % Must use "Temporally aligned motifs (W_aligned)"

%% Run t-SNE
Ytsne = motifClusterTsne(Ws); 

%% Visualize and Print
for j = 1:numel(cluster_idxC)
    motifVisualizationTsneFunc(Ytsne, cluster_idxC{1, j}, 1)
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function motifVisualizationTsneFunc(Y, cluster_idx, figSaveLogic)
% Visualize t-SNE embedding of motifs, color-coded by cluster
% Ws         – [P × nMotifs × L] matrix (temporally aligned motifs)
% cluster_idx – [nMotifs × 1] vector of cluster assignments

saveFigDir = '/Volumes/buschman/Rodent Data/dualImaging_parkj/collectFigure/visualCluster';
if ispc
    saveFigDir = 'Z:\Rodent Data\dualImaging_parkj\collectFigure\visualCluster';
end

% Identify clusters
cluster_ids = unique(cluster_idx);
K = numel(cluster_ids);

% % Compute centroids in t-SNE space
centroids = zeros(K, 2);
for c = 1:K
    centroids(c, :) = median(Y(cluster_idx == cluster_ids(c), :), 1);
end

% Normalize t-SNE coordinates to [0, 1] for consistent color mapping
Y_norm = normalize(Y, 'range');

% Identify clusters
cluster_ids = unique(cluster_idx);
K = numel(cluster_ids);

% Compute centroids (median of t-SNE coordinates)
centroids = zeros(K, 2);
for c = 1:K
    centroids(c, :) = median(Y_norm(cluster_idx == cluster_ids(c), :), 1);
end

% Convert centroids to polar angle to get hue
x_cent = centroids(:,1) - 0.5;
y_cent = centroids(:,2) - 0.5;
angles = atan2(y_cent, x_cent);  % Range: [-pi, pi]
hues = mod(angles / (2*pi), 1);  % Normalize to [0,1]

% Fixed saturation and brightness
sat = 0.65;
val = 0.95;

% Convert HSV to RGB
hsv_colors = hsv2rgb([hues, repmat(sat, K, 1), repmat(val, K, 1)]);

% Map cluster colors to each motif
cluster_color_map = containers.Map(num2cell(cluster_ids), num2cell(hsv_colors, 2));
motif_colors = zeros(size(cluster_idx,1), 3);
for j = 1:length(cluster_idx)
    motif_colors(j, :) = cluster_color_map(cluster_idx(j));
end

% Plot
figure;
scatter(Y(:,1), Y(:,2), 10, motif_colors, 'filled'); hold on;
title(sprintf('t-SNE: %d Clusters', K));
axis equal off;

% Annotate cluster centers using original Y coordinates (not normalized)
centroids_plot = zeros(K, 2);
for c = 1:K
    centroids_plot(c, :) = median(Y(cluster_idx == cluster_ids(c), :), 1);
end
for c = 1:K
    text(centroids_plot(c,1), centroids_plot(c,2), sprintf('c%d', cluster_ids(c)), ...
        'Color', 'k', 'FontSize', 10, 'HorizontalAlignment', 'center');
end

% Print (Optional)
if figSaveLogic
    timestampStr = datestr(now, 'mmddyy_HHMMSS');  % Date and time string
    figSaveName = sprintf('clusterTsne_total%dmotifs_%s', K, timestampStr);
    print(fullfile(saveFigDir, figSaveName),'-dpdf','-painters','-bestfit')
end

end












