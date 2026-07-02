function plotDAsC_heatmap(DAsC, blk, varargin)
% plotDAsC_heatmap
%
% Plot motif-projected dopamine traces as a motif x time heatmap.
%
% Input:
%   DAsC{1, blk} : K x frames motif-projected DA traces
%   DAsC{2, blk} : 1 x frames timestamps, optional
%
% Example:
%   plotDAsC_heatmap(DAsC, 1);
%
% Optional:
%   plotDAsC_heatmap(DAsC, 1, 'zscoreRows', true);
%   plotDAsC_heatmap(DAsC, 1, 'titleStr', 'm1045 block 1 projected DA');
%   plotDAsC_heatmap(DAsC, 1, 'useTime', true);

%% Parse inputs

p = inputParser;
addParameter(p, 'zscoreRows', false, @islogical);
addParameter(p, 'useTime', true, @islogical);
addParameter(p, 'titleStr', '', @ischar);
addParameter(p, 'climPrctile', [1 99], @isnumeric);
addParameter(p, 'showColorbar', true, @islogical);
parse(p, varargin{:});

zscoreRows = p.Results.zscoreRows;
useTime = p.Results.useTime;
titleStr = p.Results.titleStr;
climPrctile = p.Results.climPrctile;
showColorbar = p.Results.showColorbar;

%% Check data

if blk > size(DAsC, 2)
    error('Requested block %d, but DAsC only has %d blocks.', blk, size(DAsC, 2));
end

if isempty(DAsC{1, blk})
    error('DAsC{1, %d} is empty.', blk);
end

DAmat = double(DAsC{1, blk});   % K x frames

[K, nFrames] = size(DAmat);

%% Optional row-wise z-score

if zscoreRows
    DAmatPlot = rowZscore_omitnan(DAmat);
    colorLabel = 'Projected DA, row z-score';
else
    DAmatPlot = DAmat;
    colorLabel = 'Projected DA dF/F';
end

%% X-axis: timestamps if available, otherwise frame index

xVals = 1:nFrames;
xLabelStr = 'Frame';

if useTime && size(DAsC, 1) >= 2 && ~isempty(DAsC{2, blk})

    t = double(DAsC{2, blk}(:)');

    if numel(t) >= nFrames
        t = t(1:nFrames);
    end

    if numel(t) == nFrames && any(isfinite(t))
        xVals = t;
        xLabelStr = 'Time, sec';
    end
end

%% Plot

figure;
imagesc(xVals, 1:K, DAmatPlot);
axis xy;

xlabel(xLabelStr);
ylabel('Calcium motif #');

if isempty(titleStr)
    if zscoreRows
        title(sprintf('Motif-projected DA traces, block %d, row z-scored', blk));
    else
        title(sprintf('Motif-projected DA traces, block %d', blk));
    end
else
    title(titleStr);
end

%% Color scaling

finiteVals = DAmatPlot(isfinite(DAmatPlot));

if ~isempty(finiteVals) && numel(climPrctile) == 2
    cVals = prctile(finiteVals, climPrctile);

    if cVals(1) < cVals(2)
        caxis(cVals);
    end
end

%% Colorbar

if showColorbar
    cb = colorbar;
    ylabel(cb, colorLabel);
end

box off;

end


%% ========================================================================
function Z = rowZscore_omitnan(X)
% rowZscore_omitnan
%
% Row-wise z-scoring with NaN handling.
%
% Input:
%   X : rows x columns
%
% Output:
%   Z : rows x columns

Z = NaN(size(X));

for i = 1:size(X, 1)

    x = X(i, :);

    mu = mean(x, 'omitnan');
    sd = std(x, 0, 'omitnan');

    if isfinite(sd) && sd > 0
        Z(i, :) = (x - mu) ./ sd;
    end
end

end