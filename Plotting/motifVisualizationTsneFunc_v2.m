function motifVisualizationTsneFunc_v2(Y, cluster_idx, figSaveLogic, saveNameKeyword)
% Visualize t-SNE embedding of motifs, color-coded by cluster.
%
% This version:
%   1. Assigns one RGB color per cluster.
%   2. Uses centroid angle to choose initial hue.
%   3. Detects near-overlapping hues and separates them.
%   4. Handles near-center clusters with deterministic fallback hues.
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

%% -------------------------------------------------------
%  0) Basic checks
% --------------------------------------------------------

if size(Y,2) ~= 2
    error('Y must be [nMotifs × 2].');
end

cluster_idx = cluster_idx(:);

if size(Y,1) ~= numel(cluster_idx)
    error('size(Y,1) must match numel(cluster_idx).');
end

%% -------------------------------------------------------
%  1) Remove NaN cluster assignments
% --------------------------------------------------------

validMask = ~isnan(cluster_idx);

if ~any(validMask)
    warning('All cluster_idx values are NaN. Nothing to plot.');
    return;
end

Y = Y(validMask, :);
cluster_idx = cluster_idx(validMask);

%% -------------------------------------------------------
%  2) Identify clusters
% --------------------------------------------------------

cluster_ids = unique(cluster_idx(:), 'stable');
K = numel(cluster_ids);

%% -------------------------------------------------------
%  3) Normalize t-SNE coordinates for color computation
% --------------------------------------------------------

Y_norm = normalize(Y, 'range');

%% -------------------------------------------------------
%  4) Compute cluster centroids
% --------------------------------------------------------

centroids_norm = nan(K, 2);
centroids_plot = nan(K, 2);
clusterSizes   = nan(K, 1);

for c = 1:K
    mask = cluster_idx == cluster_ids(c);

    centroids_norm(c,:) = median(Y_norm(mask,:), 1, 'omitnan');
    centroids_plot(c,:) = median(Y(mask,:), 1, 'omitnan');
    clusterSizes(c) = sum(mask);
end

%% -------------------------------------------------------
%  5) Initial angle-based hue assignment
% --------------------------------------------------------

x_cent = centroids_norm(:,1) - 0.5;
y_cent = centroids_norm(:,2) - 0.5;

angles = atan2(y_cent, x_cent);      % [-pi, pi]
hues   = mod(angles / (2*pi), 1);    % [0, 1]

r_cent = sqrt(x_cent.^2 + y_cent.^2);

%% -------------------------------------------------------
%  6) Handle near-center clusters
% --------------------------------------------------------
%
% For clusters close to the center, angle is poorly defined.
% Assign deterministic fallback hues.

centerThresh = 0.08;
centerMask = r_cent < centerThresh;

if any(centerMask)

    centerIdx = find(centerMask);

    fallbackHues = localGenerateSeparatedHues(numel(centerIdx), 0.07);

    [~, sortOrder] = sort(cluster_ids(centerIdx));
    centerIdxSorted = centerIdx(sortOrder);

    hues(centerIdxSorted) = fallbackHues(:);
end

%% -------------------------------------------------------
%  7) Separate clusters with overlapping / too-similar hues
% --------------------------------------------------------
%
% This is the key added step.
%
% If two clusters have very similar angular hue, reassign one of them to
% the available hue that is maximally separated from already-used hues.

minHueSep = 0.15;  % 0.10 hue units = 36 degrees on HSV wheel
hues = localSeparateCloseHues(hues, cluster_ids, r_cent, clusterSizes, minHueSep);

%% -------------------------------------------------------
%  8) Convert hues to RGB
% --------------------------------------------------------

sat = 0.70;
val = 0.95;

cluster_colors = hsv2rgb([hues(:), repmat(sat, K, 1), repmat(val, K, 1)]);

%% -------------------------------------------------------
%  9) Plot each cluster separately
% --------------------------------------------------------

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

%% -------------------------------------------------------
%  10) Annotate cluster centers
% --------------------------------------------------------

for c = 1:K
    text(centroids_plot(c,1), centroids_plot(c,2), ...
        sprintf('c%d', cluster_ids(c)), ...
        'Color', 'k', ...
        'FontSize', 10, ...
        'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold', ...
        'FontAngle', 'italic');
end

%% -------------------------------------------------------
%  11) Optional color audit in command window
% --------------------------------------------------------

fprintf('\nCluster color assignment:\n');
for c = 1:K
    fprintf('  c%d: hue = %.3f, size = %d, radius = %.3f\n', ...
        cluster_ids(c), hues(c), clusterSizes(c), r_cent(c));
end

%% -------------------------------------------------------
%  12) Save optionally
% --------------------------------------------------------

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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function huesOut = localSeparateCloseHues(huesIn, cluster_ids, r_cent, clusterSizes, minHueSep)
% localSeparateCloseHues
%
% Reassigns hues when clusters are too close on the HSV color wheel.
%
% Strategy:
%   - Keep the most stable/important clusters first.
%   - Stability priority:
%       1. farther from center, because angle is more meaningful
%       2. larger cluster size
%       3. lower cluster ID
%   - If a new hue is too close to an already-used hue, assign the most
%     separated candidate hue from a dense grid.

huesOut = huesIn(:);
K = numel(huesOut);

if K <= 1
    return
end

% Sort by radial distance, then cluster size.
% Larger radius and larger size get priority to keep their angle-derived hue.
sortMat = [-r_cent(:), -clusterSizes(:), cluster_ids(:)];
[~, order] = sortrows(sortMat, [1 2 3]);

usedHues = nan(0,1);
finalHues = nan(K,1);

candidateGrid = linspace(0, 1, 361)';
candidateGrid(end) = [];  % remove duplicate 1 == 0

for oo = 1:K

    c = order(oo);
    proposedHue = huesOut(c);

    if isempty(usedHues)
        finalHues(c) = proposedHue;
        usedHues(end+1,1) = proposedHue; %#ok<AGROW>
        continue
    end

    dToUsed = localCircularHueDistance(proposedHue, usedHues);

    if min(dToUsed) >= minHueSep

        finalHues(c) = proposedHue;
        usedHues(end+1,1) = proposedHue; %#ok<AGROW>

    else

        % Find candidate hue that maximizes distance from used hues.
        distMat = zeros(numel(candidateGrid), numel(usedHues));

        for uu = 1:numel(usedHues)
            distMat(:,uu) = localCircularHueDistance(candidateGrid, usedHues(uu));
        end

        minDistToUsed = min(distMat, [], 2);

        % Prefer candidates far from used hues.
        maxMinDist = max(minDistToUsed);
        bestI = find(abs(minDistToUsed - maxMinDist) < 1e-12);

        % Among equally good candidates, choose closest to proposed hue.
        dToProposed = localCircularHueDistance(candidateGrid(bestI), proposedHue);
        [~, bestLocalI] = min(dToProposed);

        newHue = candidateGrid(bestI(bestLocalI));

        finalHues(c) = newHue;
        usedHues(end+1,1) = newHue; %#ok<AGROW>

        fprintf('Hue collision fixed: cluster c%d reassigned %.3f -> %.3f\n', ...
            cluster_ids(c), proposedHue, newHue);
    end
end

huesOut = finalHues;

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = localCircularHueDistance(h1, h2)
% Circular distance on HSV hue wheel.
%
% h1, h2 in [0, 1].
% Output is in [0, 0.5].

dRaw = abs(h1(:) - h2(:)');
d = min(dRaw, 1 - dRaw);

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hues = localGenerateSeparatedHues(n, offset)
% Generate n evenly spaced hues with optional circular offset.

if nargin < 2 || isempty(offset)
    offset = 0;
end

if n <= 0
    hues = [];
    return
end

hues = mod((0:n-1)' ./ n + offset, 1);

end