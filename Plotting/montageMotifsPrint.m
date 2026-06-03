function montageMotifsPrint(motif_w, nanpxs, varargin)
% montageMotifsPrint  Display and optionally print DA motif montages.
%
% Usage:
%   montageMotifsPrint(motif_w, nanpxs)
%   montageMotifsPrint(motif_w, nanpxs, 'motifsPerFig', 3)
%   montageMotifsPrint(motif_w, nanpxs, 'printLogic', false)
%   montageMotifsPrint(motif_w, nanpxs, 'selectMotifs', [1 3 6 7 9 12 13 14 15])
%   montageMotifsPrint(motif_w, nanpxs, 'prctileRange', [1 99.5])
%
% Inputs:
%   motif_w      : [P x K x L] array of spatiotemporal motifs
%                  pixels x motifs x frames/lags
%   nanpxs       : NaN pixel information used by conditionDffMat
%
% Name-value pairs:
%   'motifsPerFig' : number of motifs per figure, default = 1
%   'printLogic'   : logical scalar, whether to print to PDF, default = true
%   'selectMotifs' : vector of original motif IDs to display, default = all
%   'prctileRange' : percentile range for per-motif scaling, default = [1 99.5]
%   'scaleMode'    : 'perMotif', 'global', or 'none', default = 'perMotif'
%   'displayRange' : display range after scaling, default = [0 1]
%
% Notes:
%   - Per-motif scaling rescales each motif across all lags using the requested
%     percentile range.
%   - This prevents one high-amplitude motif from saturating the whole montage.
%   - Printed filenames preserve original motif IDs.

saveFigDir = '/Volumes/buschman/Rodent Data/dualImaging_parkj/collectFigure/montage';
if ispc
    saveFigDir = 'Z:\Rodent Data\dualImaging_parkj\collectFigure\montage';
end

if exist(saveFigDir, 'dir') ~= 7
    mkdir(saveFigDir);
end

% ---------------- Parse inputs ----------------
p = inputParser;
p.addParameter('motifsPerFig', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('printLogic', true, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('selectMotifs', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('prctileRange', [1 99.5], @(x) isnumeric(x) && numel(x) == 2 && x(1) < x(2));
p.addParameter('scaleMode', 'perMotif', @(x) ischar(x) || isstring(x));
p.addParameter('displayRange', [0 1], @(x) isnumeric(x) && numel(x) == 2 && x(1) < x(2));
p.addParameter('colormapName', 'magma', @(x) ischar(x) || isstring(x) || isnumeric(x));
p.addParameter('figNamePrefix', 'Motifs', @(x) ischar(x) || isstring(x));

p.parse(varargin{:});

motifsPerFig = p.Results.motifsPerFig;
printLogic   = logical(p.Results.printLogic);
selectMotifs = p.Results.selectMotifs;
prctileRange = p.Results.prctileRange;
scaleMode    = char(p.Results.scaleMode);
displayRange = p.Results.displayRange;
colormapName = p.Results.colormapName;
figNamePrefix = char(p.Results.figNamePrefix);

[P, nMotifTotal, nFrames] = size(motif_w);

if isempty(selectMotifs)
    selectMotifs = 1:nMotifTotal;
else
    selectMotifs = unique(selectMotifs(:))';
    if any(selectMotifs < 1) || any(selectMotifs > nMotifTotal)
        error('selectMotifs contains indices outside valid motif range 1:%d.', nMotifTotal);
    end
end

% Restrict displayed motifs, but preserve original IDs in selectMotifs
motif_w_sel = motif_w(:, selectMotifs, :);
nMotifSel = numel(selectMotifs);
nFigures = ceil(nMotifSel / motifsPerFig);

% Optional global scaling range
switch lower(scaleMode)
    case 'global'
        vals = motif_w_sel(:);
        vals = vals(isfinite(vals));
        vals = vals(vals ~= 0);

        if isempty(vals)
            globalLow = 0;
            globalHigh = 1;
        else
            globalLow = prctile(vals, prctileRange(1));
            globalHigh = prctile(vals, prctileRange(2));

            if globalHigh <= globalLow
                globalLow = min(vals);
                globalHigh = max(vals);
            end

            if globalHigh <= globalLow
                globalLow = 0;
                globalHigh = 1;
            end
        end

    case {'permotif', 'none'}
        globalLow = [];
        globalHigh = [];

    otherwise
        error('scaleMode must be ''perMotif'', ''global'', or ''none''.');
end

% ---------------- Main loop ----------------
for figIdx = 1:nFigures

    % Indices within selected motif list
    startIdxLocal = (figIdx - 1) * motifsPerFig + 1;
    endIdxLocal   = min(figIdx * motifsPerFig, nMotifSel);
    motifsThisFig = endIdxLocal - startIdxLocal + 1;

    % Original motif IDs shown in this figure
    motifIDsThisFig = selectMotifs(startIdxLocal:endIdxLocal);

    allMotifs = zeros(64, 64, 1, motifsThisFig * nFrames);

    for i = 1:motifsThisFig

        motifIdxLocal = startIdxLocal + i - 1;

        if P == 64 * 64
            motif = reshape(squeeze(motif_w_sel(:, motifIdxLocal, :)), 64, 64, []);
        else
            motif = conditionDffMat(squeeze(motif_w_sel(:, motifIdxLocal, :))', nanpxs);
        end

        % ---------------- Scaling ----------------
        switch lower(scaleMode)

            case 'permotif'
                motif = localScaleByPercentile(motif, prctileRange);

            case 'global'
                motif = localScaleByFixedRange(motif, globalLow, globalHigh);

            case 'none'
                % Leave motif as-is.

        end

        motif(~isfinite(motif)) = 0;

        allMotifs(:, :, :, (i-1)*nFrames + (1:nFrames)) = motif;
    end

    % Show montage
    h = figure('Color', 'k');
    montage(allMotifs, ...
        'Size', [motifsThisFig, nFrames], ...
        'DisplayRange', displayRange);

    localApplyColormap(colormapName);
    set(gca, 'Color', 'k');
    set(gcf, 'InvertHardcopy', 'off');

    title(sprintf('Motifs %s', compressMotifIDs(motifIDsThisFig)), ...
        'Color', 'w', ...
        'Interpreter', 'none');

    % Print montage
    if printLogic
        timestampStr = datestr(now, 'mmddyy_HHMMSS');
        motifLabelStr = compressMotifIDs(motifIDsThisFig);
        figSaveName = sprintf('%s_%s_%s', figNamePrefix, motifLabelStr, timestampStr);
        
        % Print montage
        set(h, 'InvertHardcopy', 'off');  % preserve black background
        print(h, fullfile(saveFigDir, figSaveName), ...
            '-painters', '-bestfit', '-dpdf');

        % If exportgraphics is unavailable in your MATLAB version, use:
        % print(h, fullfile(saveFigDir, figSaveName), '-dpdf', '-bestfit');
    end
end

end


function motifScaled = localScaleByPercentile(motif, prctileRange)

vals = motif(:);
vals = vals(isfinite(vals));
vals = vals(vals ~= 0);

if isempty(vals)
    motifScaled = zeros(size(motif), 'like', motif);
    return
end

lo = prctile(vals, prctileRange(1));
hi = prctile(vals, prctileRange(2));

if hi <= lo
    lo = min(vals);
    hi = max(vals);
end

if hi <= lo
    motifScaled = zeros(size(motif), 'like', motif);
    return
end

motifScaled = (motif - lo) ./ (hi - lo);
motifScaled(motifScaled < 0) = 0;
motifScaled(motifScaled > 1) = 1;

end


function motifScaled = localScaleByFixedRange(motif, lo, hi)

if hi <= lo
    motifScaled = zeros(size(motif), 'like', motif);
    return
end

motifScaled = (motif - lo) ./ (hi - lo);
motifScaled(motifScaled < 0) = 0;
motifScaled(motifScaled > 1) = 1;

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

function localApplyColormap(colormapName)
% localApplyColormap
%
% Applies a colormap from either:
%   - a string/char name, e.g. 'magma', 'parula', 'hot', 'turbo'
%   - an explicit N × 3 colormap matrix

if isnumeric(colormapName)
    colormap(colormapName);
    return
end

cmapName = char(colormapName);

try
    colormap(feval(cmapName, 256));
catch
    try
        colormap(cmapName);
    catch
        warning('Could not apply colormap "%s". Falling back to parula.', cmapName);
        colormap(parula);
    end
end

end