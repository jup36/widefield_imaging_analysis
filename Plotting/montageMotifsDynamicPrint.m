function montageMotifsDynamicPrint(motif_w, nanpxs, varargin)
% montageMotifsDynamicPrint
%
% Display and optionally print spatiotemporal motif montages with scaling
% designed to highlight temporal dynamics.
%
% INPUT
%   motif_w : [P x K x L] motifs, pixels x motifs x frames/lags
%   nanpxs  : NaN pixel info used by conditionDffMat
%
% Name-value pairs:
%   'motifsPerFig' : number of motifs per figure, default = 10
%   'printLogic'   : whether to print PDF, default = true
%   'selectMotifs' : original motif IDs to display, default = all
%   'scaleMode'    : 'perMotifPercentile', 'globalPercentile', or 'fixed'
%   'climPrctile'  : percentile range for scaling, default = [1 99.5]
%   'fixedCLim'    : fixed color limits if scaleMode = 'fixed', default = [0 0.3]
%   'saveFigDir'   : output directory
%   'figPrefix'    : output filename prefix
%   'showFrameLabels' : show lag/frame labels on first row, default = false
%   'showTitle'     : show figure title, default = false

%% Default save directory
saveFigDir = '/Volumes/buschman/Rodent Data/dualImaging_parkj/collectFigure/montage_dynamic';
if ispc
    saveFigDir = 'Z:\Rodent Data\dualImaging_parkj\collectFigure\montage_dynamic';
end

%% Parse inputs
p = inputParser;
p.addParameter('motifsPerFig', 10, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('printLogic', true, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('selectMotifs', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('scaleMode', 'perMotifPercentile', @(x) ischar(x) || isstring(x));
p.addParameter('climPrctile', [1 99.5], @(x) isnumeric(x) && numel(x) == 2);
p.addParameter('fixedCLim', [0 0.3], @(x) isnumeric(x) && numel(x) == 2);
p.addParameter('saveFigDir', saveFigDir, @(x) ischar(x) || isstring(x));
p.addParameter('figPrefix', 'CAmotifs_dynamic', @(x) ischar(x) || isstring(x));
p.addParameter('showFrameLabels', false, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('showTitle', false, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));

p.parse(varargin{:});

motifsPerFig    = p.Results.motifsPerFig;
printLogic      = logical(p.Results.printLogic);
selectMotifs    = p.Results.selectMotifs;
scaleMode       = char(p.Results.scaleMode);
climPrctile     = p.Results.climPrctile;
fixedCLim       = p.Results.fixedCLim;
saveFigDir      = char(p.Results.saveFigDir);
figPrefix       = char(p.Results.figPrefix);
showFrameLabels = logical(p.Results.showFrameLabels);
showTitle       = logical(p.Results.showTitle);

if exist(saveFigDir, 'dir') ~= 7
    mkdir(saveFigDir);
end

[P, nMotifTotal, nFrames] = size(motif_w);

if isempty(selectMotifs)
    selectMotifs = 1:nMotifTotal;
else
    selectMotifs = unique(selectMotifs(:))';
    if any(selectMotifs < 1) || any(selectMotifs > nMotifTotal)
        error('selectMotifs contains indices outside the valid motif range 1:%d.', nMotifTotal);
    end
end

nMotifSel = numel(selectMotifs);
nFigures = ceil(nMotifSel / motifsPerFig);

%% Precompute global color limit if requested
if strcmpi(scaleMode, 'globalPercentile')
    allVals = motif_w(:);
    allVals = allVals(~isnan(allVals) & isfinite(allVals));

    if isempty(allVals)
        globalCLim = [0 1];
    else
        globalCLim = prctile(allVals, climPrctile);
        if globalCLim(1) == globalCLim(2)
            globalCLim = [min(allVals), max(allVals)];
        end
        if globalCLim(1) == globalCLim(2)
            globalCLim = [0 1];
        end
    end
else
    globalCLim = [];
end

%% Main loop
for figIdx = 1:nFigures

    startIdxLocal = (figIdx - 1) * motifsPerFig + 1;
    endIdxLocal   = min(figIdx * motifsPerFig, nMotifSel);

    motifIDsThisFig = selectMotifs(startIdxLocal:endIdxLocal);
    motifsThisFig   = numel(motifIDsThisFig);

    % Use a grid-friendly figure aspect ratio
    figW = 0.95;
    figH = min(0.90, max(0.25, 0.075 * motifsThisFig));

    h = figure('Color', 'k', ...
        'Units', 'normalized', ...
        'Position', [0.02 0.05 figW figH], ...
        'InvertHardcopy', 'off');

    tl = tiledlayout(motifsThisFig, nFrames, ...
        'TileSpacing', 'none', ...
        'Padding', 'none');

    if showTitle
        title(tl, sprintf('Motifs %s | scale: %s', ...
            compressMotifIDs(motifIDsThisFig), scaleMode), ...
            'Color', 'w', ...
            'FontSize', 14, ...
            'FontWeight', 'bold');
    end

    for i = 1:motifsThisFig

        motifID = motifIDsThisFig(i);

        % Convert motif vector/tensor back to 64 x 64 x frames
        if P == 64 * 64
            motif = reshape(squeeze(motif_w(:, motifID, :)), 64, 64, []);
        else
            motif = conditionDffMat(squeeze(motif_w(:, motifID, :))', nanpxs);
        end

        % Determine color limits
        switch lower(scaleMode)

            case 'permotifpercentile'
                vals = motif(:);
                vals = vals(~isnan(vals) & isfinite(vals));

                if isempty(vals)
                    clim = [0 1];
                else
                    clim = prctile(vals, climPrctile);

                    if clim(1) == clim(2)
                        clim = [min(vals), max(vals)];
                    end

                    if clim(1) == clim(2)
                        clim = [0 1];
                    end
                end

            case 'globalpercentile'
                clim = globalCLim;

            case 'fixed'
                clim = fixedCLim;

            otherwise
                error('Unknown scaleMode: %s', scaleMode);
        end

        for t = 1:nFrames

            ax = nexttile;
            imagesc(ax, motif(:, :, t), clim);

            try
                colormap(ax, magma);
            catch
                colormap(ax, hot(256));
            end

            axis(ax, 'off');
            axis(ax, 'tight');
            pbaspect(ax, [1 1 1]);

            set(ax, ...
                'Color', 'k', ...
                'XTick', [], ...
                'YTick', [], ...
                'LooseInset', [0 0 0 0], ...
                'PositionConstraint', 'outerposition');

            % Motif ID label only on first frame of each motif row
            if t == 1
                text(ax, 2, 6, sprintf('%d', motifID), ...
                    'Color', 'w', ...
                    'FontSize', 12, ...
                    'FontWeight', 'bold', ...
                    'FontAngle', 'italic');
            end

            % Optional frame/lag labels
            if showFrameLabels && i == 1
                text(ax, 32, 4, sprintf('%d', t), ...
                    'Color', [0.85 0.85 0.85], ...
                    'FontSize', 8, ...
                    'FontWeight', 'bold', ...
                    'HorizontalAlignment', 'center');
            end
        end
    end

    %% Print
    if printLogic
        timestampStr = datestr(now, 'mmddyy_HHMMSS');
        motifLabelStr = compressMotifIDs(motifIDsThisFig);

        figSaveName = sprintf('%s_%s_%s', ...
            figPrefix, motifLabelStr, timestampStr);

        print(h, fullfile(saveFigDir, figSaveName), ...
            '-dpdf', '-painters', '-bestfit');
    end
end

end


function outStr = compressMotifIDs(ids)
% compressMotifIDs  Convert motif ID vector into compact range string.
%
% Example:
%   [1 2 3 4 5 6 7 9 12 13 14 15] -> '1-7_9_12-15'

ids = unique(ids(:))';

if isempty(ids)
    outStr = '';
    return;
end

rangeParts = {};
rangeStart = ids(1);
prevVal = ids(1);

for ii = 2:numel(ids)
    if ids(ii) == prevVal + 1
        prevVal = ids(ii);
    else
        rangeParts{end+1} = localRangeToStr(rangeStart, prevVal); %#ok<AGROW>
        rangeStart = ids(ii);
        prevVal = ids(ii);
    end
end

rangeParts{end+1} = localRangeToStr(rangeStart, prevVal);
outStr = strjoin(rangeParts, '_');

end


function s = localRangeToStr(a, b)

if a == b
    s = sprintf('%d', a);
else
    s = sprintf('%d-%d', a, b);
end

end