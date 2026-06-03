function motifVisualizationUmapFunc(Y, cluster_idx, figSaveLogic, saveNameKeyword)
% Visualize UMAP embedding of motifs, color-coded by cluster
%
% Inputs:
%   Y            – [nMotifs × 2] UMAP coordinates
%   cluster_idx  – [nMotifs × 1] cluster assignments (NaN allowed; will be dropped)
%   figSaveLogic – true/false
%   saveNameKeyword – optional string/char keyword for saving

if nargin < 4
    if nargin == 3
        saveNameKeyword = '';
    else
        error("Check the number of inputs!")
    end
end

saveFigDir = '/Volumes/buschman/Rodent Data/dualImaging_parkj/collectFigure/visualCluster';
if ispc
    saveFigDir = 'Z:\Rodent Data\dualImaging_parkj\collectFigure\visualCluster';
end

% -------------------------------------------------------
% 0) Basic checks
% -------------------------------------------------------
if size(Y, 2) ~= 2
    error('Y must be [nMotifs × 2]. Current size is [%d × %d].', size(Y,1), size(Y,2));
end

if size(cluster_idx, 2) > 1
    cluster_idx = cluster_idx(:);
end

if size(Y, 1) ~= numel(cluster_idx)
    error('Number of rows in Y (%d) must match numel(cluster_idx) (%d).', ...
        size(Y,1), numel(cluster_idx));
end

% -------------------------------------------------------
% 1) Remove NaN cluster assignments
% -------------------------------------------------------
validMask = ~isnan(cluster_idx);

if ~any(validMask)
    warning('All cluster_idx values are NaN. Nothing to plot.');
    return;
end

Y = Y(validMask, :);
cluster_idx = cluster_idx(validMask);

% -------------------------------------------------------
% 2) Identify clusters
% -------------------------------------------------------
cluster_ids = unique(cluster_idx);
K = numel(cluster_ids);

% -------------------------------------------------------
% 3) Normalize UMAP coordinates for color computation
% -------------------------------------------------------
Y_norm = normalize(Y, 'range');

% -------------------------------------------------------
% 4) Compute centroids in normalized space for hue assignment
% -------------------------------------------------------
centroids = zeros(K, 2);

for c = 1:K
    centroids(c, :) = median(Y_norm(cluster_idx == cluster_ids(c), :), 1);
end

% Convert centroid angle around embedding center to hue
x_cent = centroids(:,1) - 0.5;
y_cent = centroids(:,2) - 0.5;

angles = atan2(y_cent, x_cent);      % [-pi, pi]
hues   = mod(angles / (2*pi), 1);    % [0, 1]

% Fixed saturation and value
sat = 0.65;
val = 0.95;

hsv_colors = hsv2rgb([hues, repmat(sat, K, 1), repmat(val, K, 1)]);

% -------------------------------------------------------
% 5) Assign motif colors
% -------------------------------------------------------
motif_colors = zeros(size(cluster_idx, 1), 3);

for c = 1:K
    mask = cluster_idx == cluster_ids(c);
    motif_colors(mask, :) = repmat(hsv_colors(c, :), sum(mask), 1);
end

% -------------------------------------------------------
% 6) Plot
% -------------------------------------------------------
figure;
scatter(Y(:,1), Y(:,2), 10, motif_colors, 'filled');
hold on;

title(sprintf('UMAP: %d Clusters', K));
axis equal off;

% -------------------------------------------------------
% 7) Annotate cluster centers using original UMAP coordinates
% -------------------------------------------------------
centroids_plot = zeros(K, 2);

for c = 1:K
    centroids_plot(c, :) = median(Y(cluster_idx == cluster_ids(c), :), 1);
end

for c = 1:K
    text(centroids_plot(c,1), centroids_plot(c,2), ...
        sprintf('c%d', cluster_ids(c)), ...
        'Color', 'k', ...
        'FontSize', 10, ...
        'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold', ...
        'FontAngle', 'italic');
end

% -------------------------------------------------------
% 8) Save optionally
% -------------------------------------------------------
if figSaveLogic
    if ~exist(saveFigDir, 'dir')
        mkdir(saveFigDir);
    end

    timestampStr = datestr(now, 'mmddyy_HHMMSS');

    if ~isempty(saveNameKeyword)
        figSaveName = sprintf('clusterUmap_%s_total%dclusters_%s', ...
            saveNameKeyword, K, timestampStr);
    else
        figSaveName = sprintf('clusterUmap_total%dclusters_%s', ...
            K, timestampStr);
    end

    print(fullfile(saveFigDir, figSaveName), ...
        '-dpdf', '-painters', '-bestfit');
end

end