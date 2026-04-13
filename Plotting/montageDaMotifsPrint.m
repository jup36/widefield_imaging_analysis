function montageDaMotifsPrint(motif_w, nanpxs, varargin)
% montageDaMotifsPrint  Display and optionally print DA motif montages.
%
% Usage:
%   montageDaMotifsPrint(motif_w, nanpxs)
%   montageDaMotifsPrint(motif_w, nanpxs, 'motifsPerFig', 3)
%   montageDaMotifsPrint(motif_w, nanpxs, 'printLogic', false)
%   montageDaMotifsPrint(motif_w, nanpxs, 'selectMotifs', [1 3 6 7 9 12 13 14 15])
%
% Inputs:
%   motif_w      : [P x K x L] array of spatiotemporal motifs
%                  (pixels x motifs x frames)
%   nanpxs       : NaN pixel information used by conditionDffMat
%
% Name-value pairs:
%   'motifsPerFig' : number of motifs per figure (default = 1)
%   'printLogic'   : logical scalar, whether to print to PDF (default = true)
%   'selectMotifs' : vector of original motif IDs to display (default = all)
%
% Notes:
%   - Printed filenames preserve the ORIGINAL motif IDs.
%   - For discontinuous motif selections, filenames compress motif ranges:
%       e.g. [1 2 3 4 5 6 7 9 12 13 14 15] -> '1-7_9_12-15'

saveFigDir = '/Volumes/buschman/Rodent Data/dualImaging_parkj/collectFigure/montage';
if ispc
    saveFigDir = 'Z:\Rodent Data\dualImaging_parkj\collectFigure\montage';
end

% ---------------- Parse inputs ----------------
p = inputParser;
p.addParameter('motifsPerFig', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('printLogic', true, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('selectMotifs', [], @(x) isempty(x) || isnumeric(x));

p.parse(varargin{:});

motifsPerFig = p.Results.motifsPerFig;
printLogic   = logical(p.Results.printLogic);
selectMotifs = p.Results.selectMotifs;

[P, nMotifTotal, nFrames] = size(motif_w);

if isempty(selectMotifs)
    selectMotifs = 1:nMotifTotal;
else
    selectMotifs = unique(selectMotifs(:))';  % row vector, sorted, unique
    if any(selectMotifs < 1) || any(selectMotifs > nMotifTotal)
        error('selectMotifs contains indices outside the valid motif range 1:%d.', nMotifTotal);
    end
end

% Restrict displayed motifs, but preserve original IDs in selectMotifs
motif_w_sel = motif_w(:, selectMotifs, :);
nMotifSel = numel(selectMotifs);
nFigures = ceil(nMotifSel / motifsPerFig);

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

        allMotifs(:, :, :, (i-1)*nFrames + (1:nFrames)) = motif;
    end

    % Show montage
    figure;
    montage(allMotifs, 'Size', [motifsThisFig, nFrames], 'DisplayRange', [0 0.3]);
    colormap magma;
    title(sprintf('Motifs %s', compressMotifIDs(motifIDsThisFig)));

    % Print montage
    if printLogic
        timestampStr = datestr(now, 'mmddyy_HHMMSS');
        motifLabelStr = compressMotifIDs(motifIDsThisFig);
        figSaveName = sprintf('DAmotifs_%s_%s', motifLabelStr, timestampStr);
        print(fullfile(saveFigDir, figSaveName), '-dpdf', '-painters', '-bestfit');
    end
end

end


function outStr = compressMotifIDs(ids)
% compressMotifIDs  Convert motif ID vector into compact range string.
%
% Example:
%   [1 2 3 4 5 6 7 9 12 13 14 15] -> '1-7_9_12-15'
%   [3] -> '3'
%   [3 5 6 8] -> '3_5-6_8'

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