function motifVisualizationTsneFunc_v2(Y, cluster_idx, figSaveLogic, saveNameKeyword)
% Visualize t-SNE embedding of motifs, color-coded by cluster.
%
% This version keeps the angle-based color assignment but guarantees that
% each cluster gets exactly one RGB color.
%
% Inputs:
%   Y               – [nMotifs × 2] t-SNE coordinates
%   cluster_idx     – [nMotifs × 1] cluster assignments, NaN allowed
%   figSaveLogic    – true/false
%   saveNameKeyword – optional save keyword

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
if size(Y,2) ~= 2
    error('Y must be [nMotifs × 2].');
end

cluster_idx = cluster_idx(:);

if size(Y,1) ~= numel(cluster_idx)
    error('size(Y,1) must match numel(cluster_idx).');
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
cluster_ids = unique(cluster_idx(:), 'stable');
K = numel(cluster_ids);

% -------------------------------------------------------
% 3) Normalize t-SNE coordinates for color computation
% -------------------------------------------------------
Y_norm = normalize(Y, 'range');

% -------------------------------------------------------
% 4) Compute cluster centroids
% -------------------------------------------------------
centroids_norm = nan(K, 2);
centroids_plot = nan(K, 2);

for c = 1:K
    mask = cluster_idx == cluster_ids(c);

    centroids_norm(c,:) = median(Y_norm(mask,:), 1, 'omitnan');
    centroids_plot(c,:) = median(Y(mask,:), 1, 'omitnan');
end

% -------------------------------------------------------
% 5) Angle-based hue assignment
% -------------------------------------------------------
x_cent = centroids_norm(:,1) - 0.5;
y_cent = centroids_norm(:,2) - 0.5;

angles = atan2(y_cent, x_cent);      % [-pi, pi]
hues   = mod(angles / (2*pi), 1);    % [0, 1]

% -------------------------------------------------------
% 6) Stabilize colors for near-center clusters
% -------------------------------------------------------
% For clusters close to the center, angle is poorly defined.
% So give them deterministic fallback hues.
r_cent = sqrt(x_cent.^2 + y_cent.^2);

centerThresh = 0.08;  % adjust if needed
centerMask = r_cent < centerThresh;

if any(centerMask)
    centerIdx = find(centerMask);

    % Assign evenly spaced fallback hues.
    % Offset avoids duplicating pure red too often.
    fallbackHues = mod(linspace(0, 1, numel(centerIdx)+1)' + 0.07, 1);
    fallbackHues = fallbackHues(1:end-1);

    % Optionally sort center clusters by cluster ID for reproducibility
    [~, sortOrder] = sort(cluster_ids(centerIdx));
    centerIdxSorted = centerIdx(sortOrder);

    hues(centerIdxSorted) = fallbackHues;
end

% Fixed saturation and value
sat = 0.65;
val = 0.95;

cluster_colors = hsv2rgb([hues, repmat(sat, K, 1), repmat(val, K, 1)]);

% -------------------------------------------------------
% 7) Plot each cluster separately
% -------------------------------------------------------
figure;
hold on;

pointSize = 10;

for c = 1:K
    mask = cluster_idx == cluster_ids(c);

    scatter(Y(mask,1), Y(mask,2), pointSize, ...
        repmat(cluster_colors(c,:), sum(mask), 1), ...
        'filled');
end

title(sprintf('t-SNE: %d Clusters', K));
axis equal off;

% -------------------------------------------------------
% 8) Annotate cluster centers
% -------------------------------------------------------
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
% 9) Save optionally
% -------------------------------------------------------
if figSaveLogic
    if ~exist(saveFigDir, 'dir')
        mkdir(saveFigDir);
    end

    timestampStr = datestr(now, 'mmddyy_HHMMSS');

    if ~isempty(saveNameKeyword)
        figSaveName = sprintf('clusterTsne_%s_total%dclusters_%s', ...
            saveNameKeyword, K, timestampStr);
    else
        figSaveName = sprintf('clusterTsne_total%dclusters_%s', ...
            K, timestampStr);
    end

    print(fullfile(saveFigDir, figSaveName), ...
        '-dpdf', '-painters', '-bestfit');
end

end