function montageMotifsPrintAdvanced(motif_w, nanpxs, varargin)
% MontageMotifsPrintAdvanced  Display and optionally print motif montages.
%
% Usage:
%   MontageMotifsPrintAdvanced(motif_w, nanpxs)
%   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'motifsPerFig', 3)
%   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'printLogic', false)
%   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'selectMotifs', [1 3 6 7 9 12 13 14 15])
%   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'prctileRange', [1 99.5])
%   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'smoothLogic', true, 'gaussianSigma', 1)
%   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'nFrames', 7)
%   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'nFrames', 7, 'frameSelectionScope', 'perMotif')
%   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'upsampleFactor', 6, 'interpMethod', 'bicubic')
%
% Inputs:
%   motif_w      : [P x K x L] array of spatiotemporal motifs
%                  pixels x motifs x frames/lags
%   nanpxs       : NaN pixel information used by conditionDffMat. Also used
%                  here to derive the dorsal-cortex mask/boundary.
%
% Name-value pairs:
%   'motifsPerFig'  : number of motifs per figure, default = 1
%   'printLogic'    : logical scalar, whether to print to PDF, default = true
%   'selectMotifs'  : vector of original motif IDs to display, default = all
%   'prctileRange'  : percentile range for per-motif scaling, default = [1 99.5]
%   'scaleMode'     : 'perMotif', 'global', or 'none', default = 'perMotif'
%   'displayRange'  : display range after scaling, default = [0 1]
%   'colormapName'  : colormap name or N x 3 matrix, default = 'magma'
%   'figNamePrefix' : output figure filename prefix, default = 'Motifs'
%   'smoothLogic'   : apply Gaussian smoothing frame-by-frame, default = false
%   'gaussianSigma' : sigma for Gaussian smoothing, default = 1
%
%   'nFrames'             : number of frames to display per motif, out of the
%                           L available. Default = [] (show all L frames).
%                           When nFrames < L, the dropped frames are chosen
%                           by an energy-based contiguous-window search (see
%                           "Frame-dropping logic" below) -- frames are never
%                           dropped out of the middle of the sequence.
%   'frameSelectionScope' : 'global' (default) or 'perMotif'.
%                           'global'   - a single common window of frames is
%                                        chosen using the pooled (summed)
%                                        energy across all displayed motifs,
%                                        so every row/tile in the montage
%                                        shares the same underlying frame
%                                        indices (recommended: keeps a
%                                        common time axis across motifs).
%                           'perMotif' - each motif independently keeps its
%                                        own best window. Frame indices may
%                                        then differ row-to-row; the kept
%                                        frame range is annotated on each row.
%
%   'cortexBoundaryLogic' : draw the dorsal-cortex boundary in gray,
%                           default = true (only applies when nanpxs actually
%                           crops the image, i.e. P ~= 64*64).
%   'boundaryColor'        : RGB triplet for the boundary line, default = [0.5 0.5 0.5]
%   'boundaryLineWidth'    : line width for the boundary, default = 1.5
%   'boundaryNumHarmonics' : number of low-frequency Fourier harmonics kept
%                           when smoothing the cortex boundary curve,
%                           default = 15. The raw pixel boundary is traced
%                           once, resampled to uniform arc-length spacing,
%                           and reconstructed from only these harmonics --
%                           fewer harmonics = smoother/rounder outline
%                           (small anatomical notches may be smoothed away);
%                           more harmonics = closer to the raw pixel shape.
%   'boundaryResamplePoints' : number of uniformly arc-length-spaced points
%                           the raw boundary is resampled to before Fourier
%                           smoothing, default = 400. Should comfortably
%                           exceed 2x boundaryNumHarmonics.
%
%   'colGapFrac' : gap between adjacent frame columns, as a fraction of the
%                  figure width, default = 0.004 (very tight). Set
%                  independently from 'rowGapFrac' (tiles are laid out with
%                  manually positioned axes rather than tiledlayout, since
%                  tiledlayout's 'TileSpacing' cannot differ by direction).
%   'rowGapFrac' : gap between adjacent motif rows, as a fraction of the
%                  figure height, default = 0.018.
%
%   'upsampleFactor' : spatial upsampling factor applied to each 64x64 frame
%                      before display/printing, default = 4 (i.e. 256x256).
%                      Interpolation is mask-aware (see Notes) to avoid
%                      bleeding intensity across the cortex boundary.
%   'interpMethod'   : interpolation method passed to imresize for the
%                      upsampling step, default = 'bicubic'.
%
% Notes:
%   - Gaussian smoothing and percentile/global scaling are computed at the
%     native 64x64 resolution (as before); spatial upsampling happens last,
%     purely for display/print rendering.
%   - Background is now white. Pixels outside the dorsal-cortex mask are
%     blended toward white (with a soft, anti-aliased edge) and baked
%     directly into an opaque RGB image, rather than relying on continuous
%     AlphaData transparency -- the 'painters' renderer used for -dpdf
%     printing does not reliably support partial alpha and can crash on it.
%   - Per-motif scaling rescales each motif across all lags using the
%     requested percentile range.
%   - Printed filenames preserve original motif IDs.
%
% Frame-dropping logic:
%   Given L available frames and a request to keep nFrames <= L, the
%   function computes, for each frame, an "energy" value (sum of squared
%   pixel values across the cortex, ignoring NaNs). It then slides a window
%   of length nFrames across the 1..L sequence and keeps whichever
%   contiguous window has the greatest total energy. Because the window is
%   contiguous, its complement (the dropped frames) is always a prefix
%   and/or suffix of the sequence -- frames are never dropped from the
%   middle, and frame order is preserved.

%% ------------------------------------------------------------------------
%  Save directory
% -------------------------------------------------------------------------

saveFigDir = '/Volumes/buschman/Rodent Data/dualImaging_parkj/collectFigure/montage';

if ispc
    saveFigDir = 'Z:\Rodent Data\dualImaging_parkj\collectFigure\montage';
end

if exist(saveFigDir, 'dir') ~= 7
    mkdir(saveFigDir);
end

%% ------------------------------------------------------------------------
%  Parse inputs
% -------------------------------------------------------------------------

p = inputParser;

p.addParameter('motifsPerFig', 1, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1);

p.addParameter('printLogic', true, ...
    @(x) islogical(x) || (isnumeric(x) && isscalar(x)));

p.addParameter('selectMotifs', [], ...
    @(x) isempty(x) || isnumeric(x));

p.addParameter('prctileRange', [1 99.5], ...
    @(x) isnumeric(x) && numel(x) == 2 && x(1) < x(2));

p.addParameter('scaleMode', 'perMotif', ...
    @(x) ischar(x) || isstring(x));

p.addParameter('displayRange', [0 1], ...
    @(x) isnumeric(x) && numel(x) == 2 && x(1) < x(2));

p.addParameter('colormapName', 'magma', ...
    @(x) ischar(x) || isstring(x) || isnumeric(x));

p.addParameter('figNamePrefix', 'Motifs', ...
    @(x) ischar(x) || isstring(x));

p.addParameter('smoothLogic', false, ...
    @(x) islogical(x) || (isnumeric(x) && isscalar(x)));

p.addParameter('gaussianSigma', 1, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);

% New: frame-count control + drop logic
p.addParameter('nFrames', [], ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 1 && mod(x,1) == 0));

p.addParameter('frameSelectionScope', 'global', ...
    @(x) ischar(x) || isstring(x));

% New: cortex boundary overlay
p.addParameter('cortexBoundaryLogic', true, ...
    @(x) islogical(x) || (isnumeric(x) && isscalar(x)));

p.addParameter('boundaryColor', [0.5 0.5 0.5], ...
    @(x) isnumeric(x) && numel(x) == 3);

p.addParameter('boundaryLineWidth', 1.5, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);

% Boundary smoothing via Fourier descriptors (see localGetCortexMaskAndBoundary):
% the raw pixel boundary is traced once, resampled to uniform arc-length
% spacing, and reconstructed from only its low-frequency components,
% guaranteeing a smooth, seamlessly closed curve.
% New: sub-pixel boundary smoothness controls (see localGetCortexMaskAndBoundary)
p.addParameter('boundaryNumHarmonics', 15, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1);

p.addParameter('boundaryResamplePoints', 400, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 16);

% New: independent horizontal/vertical tile spacing (fraction of figure
% width/height given to the gap between adjacent tiles). Unlike
% tiledlayout's 'TileSpacing', these are set independently per axis.
p.addParameter('colGapFrac', 0.004, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 0 && x < 1);

p.addParameter('rowGapFrac', 0.018, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 0 && x < 1);

% New: spatial upsampling
p.addParameter('upsampleFactor', 4, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1);

p.addParameter('interpMethod', 'bicubic', ...
    @(x) ischar(x) || isstring(x));

p.parse(varargin{:});

motifsPerFig    = p.Results.motifsPerFig;
printLogic      = logical(p.Results.printLogic);
selectMotifs    = p.Results.selectMotifs;
prctileRange    = p.Results.prctileRange;
scaleMode       = char(p.Results.scaleMode);
displayRange    = p.Results.displayRange;
colormapName    = p.Results.colormapName;
figNamePrefix   = char(p.Results.figNamePrefix);
smoothLogic     = logical(p.Results.smoothLogic);
gaussianSigma   = p.Results.gaussianSigma;

nFramesRequest      = p.Results.nFrames;
frameSelectionScope = lower(char(p.Results.frameSelectionScope));

cortexBoundaryLogic = logical(p.Results.cortexBoundaryLogic);
boundaryColor       = p.Results.boundaryColor;
boundaryLineWidth   = p.Results.boundaryLineWidth;
boundaryNumHarmonics   = p.Results.boundaryNumHarmonics;
boundaryResamplePoints = p.Results.boundaryResamplePoints;
colGapFrac          = p.Results.colGapFrac;
rowGapFrac          = p.Results.rowGapFrac;

upsampleFactor = p.Results.upsampleFactor;
interpMethod   = char(p.Results.interpMethod);

if ~ismember(frameSelectionScope, {'global', 'permotif'})
    error('frameSelectionScope must be ''global'' or ''perMotif''.');
end

%% ------------------------------------------------------------------------
%  Validate motif tensor
% -------------------------------------------------------------------------

if ~isnumeric(motif_w) || ndims(motif_w) ~= 3
    error('motif_w must be a numeric [P x K x L] array.');
end

[P, nMotifTotal, nFramesTotal] = size(motif_w);

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

%% ------------------------------------------------------------------------
%  Cortex mask + boundary (used for transparency and the gray outline)
% -------------------------------------------------------------------------

[cortexMask, cortexBoundaryXY] = localGetCortexMaskAndBoundary( ...
    nanpxs, cortexBoundaryLogic, boundaryNumHarmonics, boundaryResamplePoints);

if isempty(cortexBoundaryXY)
    cortexBoundaryLogic = false;
end

% Precompute the high-resolution alpha mask ONCE, by rasterizing the exact
% same smooth boundary curve used for the plotted line (via poly2mask)
% rather than separately bicubic-upsampling the original blocky binary
% mask. Deriving the mask and the line from two different smoothing
% methods (Fourier descriptors for one, raster blurring for the other) has
% no guarantee of agreeing pixel-for-pixel, which is what kept showing up
% as a residual fringe/staircase along the edge. Using the identical curve
% for both eliminates that mismatch entirely. This is also computed once
% here rather than per motif/frame, since it doesn't depend on the data.
hiresMaskShared = localRasterizeSmoothMask( ...
    cortexMask, cortexBoundaryXY, upsampleFactor);

%% ------------------------------------------------------------------------
%  Frame selection (drop only from the edges, keep highest-energy window)
% -------------------------------------------------------------------------

if isempty(nFramesRequest)
    nFramesKeep = nFramesTotal;
else
    nFramesKeep = nFramesRequest;
end

if nFramesKeep > nFramesTotal
    error('nFrames (%d) cannot exceed the number of available frames (%d).', ...
        nFramesKeep, nFramesTotal);
end

[motif_w_kept, keptFrameLabels] = localSelectFrameWindow( ...
    motif_w_sel, nFramesKeep, frameSelectionScope);

nMotifSel_check = size(motif_w_kept, 2); %#ok<NASGU>
nFramesShow     = size(motif_w_kept, 3);

nFigures = ceil(nMotifSel / motifsPerFig);

%% ------------------------------------------------------------------------
%  Optional global scaling range (computed on the frames actually shown)
% -------------------------------------------------------------------------

switch lower(scaleMode)

    case 'global'

        vals = motif_w_kept(:);
        vals = vals(isfinite(vals));
        vals = vals(vals ~= 0);

        if isempty(vals)
            globalLow = 0;
            globalHigh = 1;
        else
            globalLow  = prctile(vals, prctileRange(1));
            globalHigh = prctile(vals, prctileRange(2));

            if globalHigh <= globalLow
                globalLow  = min(vals);
                globalHigh = max(vals);
            end

            if globalHigh <= globalLow
                globalLow  = 0;
                globalHigh = 1;
            end
        end

    case {'permotif', 'none'}

        globalLow = [];
        globalHigh = [];

    otherwise

        error('scaleMode must be ''perMotif'', ''global'', or ''none''.');
end

%% ------------------------------------------------------------------------
%  Resolve colormap once (used to manually bake RGB below -- see note in
%  the render loop about why we don't rely on continuous AlphaData)
% -------------------------------------------------------------------------

cmap = localResolveColormap(colormapName);

%% ------------------------------------------------------------------------
%  Main montage loop
% -------------------------------------------------------------------------

for figIdx = 1:nFigures

    % Indices within selected motif list
    startIdxLocal = (figIdx - 1) * motifsPerFig + 1;
    endIdxLocal   = min(figIdx * motifsPerFig, nMotifSel);
    motifsThisFig = endIdxLocal - startIdxLocal + 1;

    % Original motif IDs shown in this figure
    motifIDsThisFig = selectMotifs(startIdxLocal:endIdxLocal);

    h = figure('Color', 'w');

    % Layout margins (figure-normalized units) reserved for the title,
    % per-column frame headers, and per-row motif labels.
    leftMarginFrac   = 0.06;
    topMarginFrac    = 0.075;
    bottomMarginFrac = 0.01;
    rightMarginFrac  = 0.01;

    nCols = nFramesShow;
    nRows = motifsThisFig;

    tileW = (1 - leftMarginFrac - rightMarginFrac - (nCols - 1) * colGapFrac) / nCols;
    tileH = (1 - topMarginFrac  - bottomMarginFrac - (nRows - 1) * rowGapFrac) / nRows;

    % Each tile calls axis(ax,'image','off') to preserve the (square)
    % image's aspect ratio. If a tile's allocated box isn't itself square,
    % MATLAB centers the square image within a letterboxed sub-region of
    % that box -- padding that colGapFrac has no control over, since it
    % isn't axes spacing at all. Fix this at the source: size the figure's
    % physical pixels so tileW/tileH normalized fractions correspond to an
    % actually-square physical tile, and axis('image') never needs to pad.
    %
    % IMPORTANT: the required figWidthPx/figHeightPx RATIO is what matters
    % for squareness, not their absolute size. Scaling by a fixed
    % targetTilePx let the absolute figure size grow with row/column count
    % (e.g. 8 rows could require a figure >2000px tall) -- if that exceeds
    % the screen, MATLAB/the OS window manager clips or rescales the
    % window to fit, which can distort the ratio non-uniformly and
    % reintroduce letterboxing (this is exactly why it worked for a single
    % row but broke for many). Instead, fix the ratio and cap the larger
    % absolute dimension to a size that reliably fits on screen.
    maxFigDimPx = 1200;
    aspectRatio = tileH / tileW; % required figWidthPx / figHeightPx

    if aspectRatio >= 1
        figWidthPx  = maxFigDimPx;
        figHeightPx = maxFigDimPx / aspectRatio;
    else
        figHeightPx = maxFigDimPx;
        figWidthPx  = maxFigDimPx * aspectRatio;
    end

    set(h, 'Units', 'pixels', 'Position', [100, 100, figWidthPx, figHeightPx]);

    for i = 1:motifsThisFig

        motifIdxLocal = startIdxLocal + i - 1;

        %% ----------------------------------------------------------------
        %  Reconstruct motif into image stack
        % -----------------------------------------------------------------

        if P == 64 * 64

            motif = reshape(squeeze(motif_w_kept(:, motifIdxLocal, :)), ...
                64, 64, []);

        else

            motif = conditionDffMat( ...
                squeeze(motif_w_kept(:, motifIdxLocal, :))', ...
                nanpxs);
        end

        %% ----------------------------------------------------------------
        %  Optional Gaussian smoothing (native resolution, before scaling)
        % -----------------------------------------------------------------

        if smoothLogic
            motif = applyImgaussfilt(motif, 'sigma', gaussianSigma);
        end

        %% ----------------------------------------------------------------
        %  Scaling
        % -----------------------------------------------------------------

        switch lower(scaleMode)

            case 'permotif'

                motif = localScaleByPercentile(motif, prctileRange);

            case 'global'

                motif = localScaleByFixedRange(motif, globalLow, globalHigh);

            case 'none'

                % Leave motif as-is.
        end

        %% ----------------------------------------------------------------
        %  Render each frame as its own tile (mask-aware upsampling +
        %  transparent background + gray cortex boundary overlay)
        % -----------------------------------------------------------------

        for fIdx = 1:nFramesShow

            frameData = motif(:, :, fIdx);
            frameData(~isfinite(frameData)) = 0;

            hiresData = localUpsampleFrameData(frameData, upsampleFactor, interpMethod);
            hiresMask = hiresMaskShared;

            % Bake the soft (anti-aliased) mask directly into an opaque RGB
            % image by blending the colormapped data with a white
            % background, rather than using continuous AlphaData. The
            % 'painters' renderer (used below for -dpdf printing) does not
            % reliably support partial/continuous transparency -- passing
            % it a non-binary AlphaData matrix can crash MATLAB. Producing
            % a fully opaque true-color image sidesteps that entirely while
            % still giving the same smooth, anti-aliased edge.
            rgbTile = localComposeRGBWithWhiteBackground( ...
                hiresData, hiresMask, displayRange, cmap);

            tileX = leftMarginFrac + (fIdx - 1) * (tileW + colGapFrac);
            tileY = 1 - topMarginFrac - i * tileH - (i - 1) * rowGapFrac;

            ax = axes('Parent', h, 'Position', [tileX, tileY, tileW, tileH]); %#ok<LAXES>
            image(ax, rgbTile);
            axis(ax, 'image', 'off');
            set(ax, 'Color', 'w');
            hold(ax, 'on');

            if cortexBoundaryLogic
                for bIdx = 1:numel(cortexBoundaryXY)
                    bxy = cortexBoundaryXY{bIdx} * upsampleFactor;
                    % plot() draws point 1 -> N but does not automatically
                    % connect N back to 1 -- explicitly close the loop by
                    % repeating the first point at the end, otherwise a gap
                    % appears at the seam (wherever the raw boundary trace
                    % happened to start).
                    bxyClosed = [bxy; bxy(1, :)];
                    plot(ax, bxyClosed(:,1), bxyClosed(:,2), '-', ...
                        'Color', boundaryColor, ...
                        'LineWidth', boundaryLineWidth);
                end
            end

            hold(ax, 'off');

            % Column header (original frame index) on the first row only,
            % and only when meaningful (i.e. a single shared frame axis).
            if i == 1 && strcmp(frameSelectionScope, 'global')
                title(ax, sprintf('f%d', keptFrameLabels.global(fIdx)), ...
                    'Color', 'k', 'FontSize', 8, 'FontWeight', 'normal');
            end

            % Row label (motif ID [+ frame range if perMotif]) on first column
            if fIdx == 1
                if strcmp(frameSelectionScope, 'permotif')
                    rowLabel = sprintf('M%d (f%d-%d)', ...
                        motifIDsThisFig(i), ...
                        keptFrameLabels.perMotif(motifIdxLocal, 1), ...
                        keptFrameLabels.perMotif(motifIdxLocal, end));
                else
                    rowLabel = sprintf('M%d', motifIDsThisFig(i));
                end
                ylabel(ax, rowLabel, 'Color', 'k', 'FontSize', 8, ...
                    'Rotation', 0, 'HorizontalAlignment', 'right', ...
                    'VerticalAlignment', 'middle', 'Visible', 'on');
                % axis('off') hides the ylabel too, so force it visible
                ax.YLabel.Visible = 'on';
            end
        end
    end

    if smoothLogic
        smoothLabel = sprintf(' | Gaussian \\sigma = %.2g', gaussianSigma);
    else
        smoothLabel = '';
    end

    if strcmp(frameSelectionScope, 'global') && nFramesShow < nFramesTotal
        frameLabel = sprintf(' | frames %d-%d of %d', ...
            keptFrameLabels.global(1), keptFrameLabels.global(end), nFramesTotal);
    elseif strcmp(frameSelectionScope, 'permotif') && nFramesShow < nFramesTotal
        frameLabel = sprintf(' | %d/%d frames (per-motif window)', ...
            nFramesShow, nFramesTotal);
    else
        frameLabel = '';
    end

    sgtitle(h, sprintf('Motifs %s%s%s', ...
        compressMotifIDs(motifIDsThisFig), frameLabel, smoothLabel), ...
        'Color', 'k', 'Interpreter', 'tex');

    %% --------------------------------------------------------------------
    %  Print montage
    % ---------------------------------------------------------------------

    if printLogic

        timestampStr  = datestr(now, 'mmddyy_HHMMSS');
        motifLabelStr = compressMotifIDs(motifIDsThisFig);

        if smoothLogic
            smoothFileStr = sprintf('_gaussSigma%.2g', gaussianSigma);
            smoothFileStr = strrep(smoothFileStr, '.', 'p');
        else
            smoothFileStr = '';
        end

        figSaveName = sprintf('%s_%s%s_%s', ...
            figNamePrefix, ...
            motifLabelStr, ...
            smoothFileStr, ...
            timestampStr);

        set(h, 'InvertHardcopy', 'off');  % preserve exact on-screen colors

        % Each tile is now a full true-color RGB bitmap (see
        % localComposeRGBWithWhiteBackground), not a small indexed image --
        % 'painters' is a pure vector renderer and is not well-suited to
        % embedding large raster content in a PDF; it has been a source of
        % crashes/instability for exactly this kind of mixed raster+vector
        % (image tiles + boundary lines) figure. 'opengl' handles this
        % mixed content far more robustly.
        print(h, fullfile(saveFigDir, figSaveName), ...
            '-painters', '-bestfit', '-dpdf');
    end
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function motifScaled = localScaleByFixedRange(motif, lo, hi)

if hi <= lo
    motifScaled = zeros(size(motif), 'like', motif);
    return
end

motifScaled = (motif - lo) ./ (hi - lo);
motifScaled(motifScaled < 0) = 0;
motifScaled(motifScaled > 1) = 1;

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = localRangeToStr(a, b)

if a == b
    s = sprintf('%d', a);
else
    s = sprintf('%d-%d', a, b);
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cmap = localResolveColormap(colormapName)
% localResolveColormap
%
% Resolves a colormap from either:
%   - a string/char name, e.g. 'magma', 'parula', 'hot', 'turbo'
%   - an explicit N x 3 colormap matrix

if isnumeric(colormapName)
    cmap = colormapName;
    return
end

cmapName = char(colormapName);

try
    cmap = feval(cmapName, 256);
catch
    try
        cmap = eval(cmapName); %#ok<EVLDIR>
    catch
        warning('Could not resolve colormap "%s". Falling back to parula.', cmapName);
        cmap = parula(256);
    end
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rgbImage = localComposeRGBWithWhiteBackground(data, alphaMask, displayRange, cmap)
% localComposeRGBWithWhiteBackground
%
% Manually maps `data` through `cmap` (using displayRange as the color
% axis limits, matching what imagesc/clim would normally do) and blends
% the result with a white background using alphaMask (a continuous [0,1]
% field), producing a fully opaque H x W x 3 RGB image.
%
% This exists specifically to AVOID passing a continuous (non-binary)
% AlphaData matrix to imagesc: the 'painters' renderer (used elsewhere in
% this function for -dpdf printing) does not reliably support partial
% transparency, and doing so can crash MATLAB. Baking the blend into plain
% RGB values sidesteps alpha compositing entirely while still producing
% the same smooth, anti-aliased edge.

lo = displayRange(1);
hi = displayRange(2);

normVal = (data - lo) ./ (hi - lo);
normVal = min(max(normVal, 0), 1);

nColors = size(cmap, 1);
colorIdx = round(normVal * (nColors - 1)) + 1;
colorIdx = min(max(colorIdx, 1), nColors);

R = reshape(cmap(colorIdx(:), 1), size(data));
G = reshape(cmap(colorIdx(:), 2), size(data));
B = reshape(cmap(colorIdx(:), 3), size(data));

a = alphaMask;

Rout = R .* a + 1 .* (1 - a);
Gout = G .* a + 1 .* (1 - a);
Bout = B .* a + 1 .* (1 - a);

rgbImage = cat(3, Rout, Gout, Bout);

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cortexMask, boundaryXY] = localGetCortexMaskAndBoundary( ...
    nanpxs, wantBoundary, numHarmonics, resampleN)
% localGetCortexMaskAndBoundary
%
% Derives a 64x64 logical dorsal-cortex mask directly from nanpxs -- the
% linear indices (or logical mask) of the non-cortex ("NaN") pixels within
% the full 64x64 = 4096 pixel grid -- and its boundary as a cell array of
% [x y] coordinate lists (native 64x64 pixel units, one cell per connected
% boundary component), smoothed via truncated Fourier descriptors.
%
% IMPORTANT: whether a real cortex mask exists is entirely a function of
% nanpxs, NOT of P (the first dimension of motif_w). motif_w can be stored
% either as [nValidPixels x K x L] (reconstructed via conditionDffMat) or
% as an already-reshaped [4096 x K x L] full-grid array -- either way, if
% nanpxs is supplied it still marks the true dorsal-cortex boundary within
% that 64x64 grid and should be used.
%
% If nanpxs is empty, there is no mask information available and the whole
% 64x64 frame is treated as valid (boundaryXY is returned empty).
%
% Boundary smoothness: rather than blurring a rasterized mask (which only
% ever anti-aliases individual pixel-step edges and struggles to remove the
% macroscopic staircase left by the native 64x64 resolution), the actual
% pixel boundary is traced once with bwboundaries, resampled to
% resampleN uniformly arc-length-spaced points (needed for a well-posed
% Fourier truncation), represented as a complex sequence z = x + 1i*y, and
% reconstructed from only its lowest numHarmonics frequency components via
% FFT/IFFT. Truncating high frequencies of a closed, periodic curve
% guarantees a result that is both smooth AND exactly seamlessly closed --
% there is no "sigma in the wrong units" pitfall here, since the smoothing
% is applied directly to the curve's shape, independent of any raster
% resolution or supersampling choice.

boundaryXY = {};

if isempty(nanpxs)
    cortexMask = true(64, 64);
    return
end

nanFlagVec = false(64 * 64, 1);

if islogical(nanpxs)
    nanFlagVec(:) = nanpxs(:);
else
    nanFlagVec(nanpxs(:)) = true;
end

cortexMask = reshape(~nanFlagVec, 64, 64);

if ~wantBoundary
    return
end

if exist('bwboundaries', 'file') ~= 2
    warning('bwboundaries (Image Processing Toolbox) not found; skipping cortex boundary overlay.');
    return
end

rawBoundaries = bwboundaries(cortexMask, 'noholes');
boundaryXY = cell(size(rawBoundaries));

for k = 1:numel(rawBoundaries)
    % bwboundaries returns [row col] = [y x]; convert to [x y]
    rawXY = [rawBoundaries{k}(:,2), rawBoundaries{k}(:,1)];

    resampledXY = localResampleClosedCurve(rawXY, resampleN);
    boundaryXY{k} = localSmoothClosedContourFourier(resampledXY, numHarmonics);
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function resampledXY = localResampleClosedCurve(xy, nResample)
% localResampleClosedCurve
%
% Resamples a closed 2D polygon (M x 2, [x y]) to nResample points spaced
% uniformly by arc length around the loop. Uniform spacing is required for
% a clean/well-posed Fourier-descriptor truncation afterward.

% Ensure explicitly closed (first point repeated at the end) for arc-length
% accumulation, then drop the duplicate after resampling.
if norm(xy(1,:) - xy(end,:)) > 1e-9
    xy = [xy; xy(1,:)];
end

% Guard against zero-length (duplicate/repeated) points -- bwboundaries can
% occasionally emit a repeated point where the traced path touches a thin
% or degenerate pixel connection. A zero-length segment here would give
% interp1 a non-strictly-increasing cumulative-distance vector, which can
% produce a small localized artifact in the curve after Fourier smoothing.
segLenRaw = sqrt(sum(diff(xy).^2, 2));
keepPoint = [true; segLenRaw > 1e-9];
xy = xy(keepPoint, :);

segLen  = sqrt(sum(diff(xy).^2, 2));
cumDist = [0; cumsum(segLen)];
totalLen = cumDist(end);

if totalLen == 0
    resampledXY = repmat(xy(1,:), nResample, 1);
    return
end

targetDist = linspace(0, totalLen, nResample + 1);
targetDist(end) = []; % drop duplicate closing point

xResampled = interp1(cumDist, xy(:,1), targetDist, 'linear');
yResampled = interp1(cumDist, xy(:,2), targetDist, 'linear');

resampledXY = [xResampled(:), yResampled(:)];

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function smoothXY = localSmoothClosedContourFourier(xy, numHarmonics)
% localSmoothClosedContourFourier
%
% Smooths a closed, uniformly arc-length-resampled 2D curve (N x 2, [x y])
% by representing it as a complex sequence z = x + 1i*y, truncating its
% discrete Fourier transform to the lowest numHarmonics positive and
% negative frequencies (plus the DC term), and inverting. This is the
% classic "Fourier descriptor" smoothing approach: it is guaranteed to
% return an exactly closed curve (the representation is inherently
% periodic) that is smooth by construction, with numHarmonics as the sole,
% resolution-independent smoothness knob (fewer harmonics = smoother/
% rounder; more harmonics = closer to the original traced shape).

N = size(xy, 1);
z = complex(xy(:,1), xy(:,2));

Z = fft(z);

numHarmonics = min(numHarmonics, floor((N - 1) / 2));

keepMask = false(N, 1);
keepMask(1) = true; % DC term
keepMask(2:(numHarmonics + 1)) = true;               % low positive frequencies
keepMask((N - numHarmonics + 1):N) = true;           % mirrored negative frequencies

Z(~keepMask) = 0;

zSmooth = ifft(Z);

smoothXY = [real(zSmooth), imag(zSmooth)];

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hiresData = localUpsampleFrameData(frameData, upsampleFactor, interpMethod)
% localUpsampleFrameData
%
% Spatially upsamples a single 64x64 frame's data for high-resolution
% display. The data is interpolated as-is, with no masking applied here --
% masking/fading to white is handled entirely by the shared alpha mask
% (see localRasterizeSmoothMask), computed once from the same smooth
% boundary curve used for the plotted line. Zeroing data outside the mask
% before interpolation was tried and rejected: combined with a separately
% smoothed alpha, it causes visible color fringing (a dark tinge from the
% hard zero shows through wherever alpha is only partially faded).

if upsampleFactor == 1
    hiresData = frameData;
    return
end

targetSize = size(frameData) * upsampleFactor;
hiresData = imresize(frameData, targetSize, interpMethod);

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hiresMask = localRasterizeSmoothMask(cortexMask, boundaryXY, upsampleFactor)
% localRasterizeSmoothMask
%
% Builds the high-resolution alpha mask by rasterizing the SAME smooth
% boundary curve(s) used for the plotted gray line (via poly2mask), rather
% than independently bicubic-upsampling the original native-resolution
% binary mask. Using two different smoothing methods for the line (Fourier
% descriptors) and the mask (raster blurring) gives no guarantee they agree
% pixel-for-pixel, which produced a visible fringe/staircase along the
% edge in earlier attempts. Deriving both from the identical curve
% eliminates that mismatch structurally. A light final Gaussian pass
% anti-aliases the (now much finer, supersampled-grid-scale) rasterization
% step, which is a much smaller and less objectionable artifact than the
% native-pixel-grid blockiness this replaces.
%
% Falls back to a plain bicubic-upsampled binary mask if no boundary curve
% is available (e.g. cortexBoundaryLogic was false, or nanpxs was empty).

targetSize = size(cortexMask) * upsampleFactor;

if isempty(boundaryXY) || exist('poly2mask', 'file') ~= 2
    hiresMask = imresize(double(cortexMask), targetSize, 'bicubic');
    hiresMask = min(max(hiresMask, 0), 1);
    return
end

maskAccum = false(targetSize);

for k = 1:numel(boundaryXY)
    xy = boundaryXY{k} * upsampleFactor;
    bw = poly2mask(xy(:,1), xy(:,2), targetSize(1), targetSize(2));
    maskAccum = maskAccum | bw;
end

hiresMask = double(maskAccum);

if exist('imgaussfilt', 'file') == 2
    hiresMask = imgaussfilt(hiresMask, 1.0); % mild anti-aliasing at supersampled-grid scale
end

hiresMask = min(max(hiresMask, 0), 1);

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [motif_w_kept, keptFrameLabels] = localSelectFrameWindow(motif_w_sel, nFramesKeep, scope)
% localSelectFrameWindow
%
% Selects, for each motif (or globally across all motifs), the contiguous
% window of nFramesKeep frames with the highest total energy, out of the
% nFramesTotal frames available. Energy for a frame is the sum of squared
% pixel values across all pixels/motifs being pooled, ignoring NaNs.
%
% Because the chosen window is always contiguous, its complement (the
% dropped frames) is always a prefix and/or suffix of the sequence -- i.e.
% frames are dropped only from the edges, never from the middle, and frame
% order is preserved.
%
% Outputs:
%   motif_w_kept    : [P x nMotifSel x nFramesKeep] tensor with the
%                     dropped frames removed.
%   keptFrameLabels : struct with fields:
%                       .global   - 1 x nFramesKeep vector of original frame
%                                   indices (valid/used when scope=='global')
%                       .perMotif - nMotifSel x nFramesKeep matrix of
%                                   original frame indices per motif (valid/
%                                   used when scope=='perMotif')

[P, nMotifSel, nFramesTotal] = size(motif_w_sel);

keptFrameLabels = struct('global', [], 'perMotif', []);

if nFramesKeep == nFramesTotal
    motif_w_kept = motif_w_sel;
    keptFrameLabels.global   = 1:nFramesTotal;
    keptFrameLabels.perMotif = repmat(1:nFramesTotal, nMotifSel, 1);
    return
end

% Energy per motif per frame: [nMotifSel x nFramesTotal]
sq = motif_w_sel .^ 2;
sq(~isfinite(sq)) = 0;
energyMotifFrame = squeeze(sum(sq, 1));       % nMotifSel x nFramesTotal
if nMotifSel == 1
    energyMotifFrame = reshape(energyMotifFrame, 1, nFramesTotal);
end

switch scope

    case 'global'

        aggregateEnergy = sum(energyMotifFrame, 1);         % 1 x nFramesTotal
        winStart = localBestWindowStart(aggregateEnergy, nFramesKeep);
        keepIdx  = winStart:(winStart + nFramesKeep - 1);

        motif_w_kept = motif_w_sel(:, :, keepIdx);
        keptFrameLabels.global   = keepIdx;
        keptFrameLabels.perMotif = repmat(keepIdx, nMotifSel, 1);

    case 'permotif'

        motif_w_kept = zeros(P, nMotifSel, nFramesKeep, 'like', motif_w_sel);
        perMotifIdx  = zeros(nMotifSel, nFramesKeep);

        for m = 1:nMotifSel
            winStart = localBestWindowStart(energyMotifFrame(m, :), nFramesKeep);
            keepIdx  = winStart:(winStart + nFramesKeep - 1);

            motif_w_kept(:, m, :) = motif_w_sel(:, m, keepIdx);
            perMotifIdx(m, :) = keepIdx;
        end

        keptFrameLabels.perMotif = perMotifIdx;

    otherwise

        error('Unknown frameSelectionScope: %s', scope);
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function winStart = localBestWindowStart(energyVec, winLen)
% localBestWindowStart
%
% Slides a window of length winLen across energyVec (1 x L) and returns the
% start index of the window with the maximum total energy. Ties are broken
% in favor of the earliest (smallest-index) window.

L = numel(energyVec);
nWindows = L - winLen + 1;

winSums = zeros(1, nWindows);
for s = 1:nWindows
    winSums(s) = sum(energyVec(s:(s + winLen - 1)));
end

[~, winStart] = max(winSums);

end


% function montageMotifsPrintAdvanced(motif_w, nanpxs, varargin)
% % MontageMotifsPrintAdvanced  Display and optionally print motif montages.
% %
% % Usage:
% %   MontageMotifsPrintAdvanced(motif_w, nanpxs)
% %   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'motifsPerFig', 3)
% %   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'printLogic', false)
% %   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'selectMotifs', [1 3 6 7 9 12 13 14 15])
% %   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'prctileRange', [1 99.5])
% %   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'smoothLogic', true, 'gaussianSigma', 1)
% %   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'nFrames', 7)
% %   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'nFrames', 7, 'frameSelectionScope', 'perMotif')
% %   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'upsampleFactor', 6, 'interpMethod', 'bicubic')
% %
% % Inputs:
% %   motif_w      : [P x K x L] array of spatiotemporal motifs
% %                  pixels x motifs x frames/lags
% %   nanpxs       : NaN pixel information used by conditionDffMat. Also used
% %                  here to derive the dorsal-cortex mask/boundary.
% %
% % Name-value pairs:
% %   'motifsPerFig'  : number of motifs per figure, default = 1
% %   'printLogic'    : logical scalar, whether to print to PDF, default = true
% %   'selectMotifs'  : vector of original motif IDs to display, default = all
% %   'prctileRange'  : percentile range for per-motif scaling, default = [1 99.5]
% %   'scaleMode'     : 'perMotif', 'global', or 'none', default = 'perMotif'
% %   'displayRange'  : display range after scaling, default = [0 1]
% %   'colormapName'  : colormap name or N x 3 matrix, default = 'magma'
% %   'figNamePrefix' : output figure filename prefix, default = 'Motifs'
% %   'smoothLogic'   : apply Gaussian smoothing frame-by-frame, default = false
% %   'gaussianSigma' : sigma for Gaussian smoothing, default = 1
% %
% %   'nFrames'             : number of frames to display per motif, out of the
% %                           L available. Default = [] (show all L frames).
% %                           When nFrames < L, the dropped frames are chosen
% %                           by an energy-based contiguous-window search (see
% %                           "Frame-dropping logic" below) -- frames are never
% %                           dropped out of the middle of the sequence.
% %   'frameSelectionScope' : 'global' (default) or 'perMotif'.
% %                           'global'   - a single common window of frames is
% %                                        chosen using the pooled (summed)
% %                                        energy across all displayed motifs,
% %                                        so every row/tile in the montage
% %                                        shares the same underlying frame
% %                                        indices (recommended: keeps a
% %                                        common time axis across motifs).
% %                           'perMotif' - each motif independently keeps its
% %                                        own best window. Frame indices may
% %                                        then differ row-to-row; the kept
% %                                        frame range is annotated on each row.
% %
% %   'cortexBoundaryLogic' : draw the dorsal-cortex boundary in gray,
% %                           default = true (only applies when nanpxs actually
% %                           crops the image, i.e. P ~= 64*64).
% %   'boundaryColor'        : RGB triplet for the boundary line, default = [0.5 0.5 0.5]
% %   'boundaryLineWidth'    : line width for the boundary, default = 1.5
% %   'boundaryNumHarmonics' : number of low-frequency Fourier harmonics kept
% %                           when smoothing the cortex boundary curve,
% %                           default = 15. The raw pixel boundary is traced
% %                           once, resampled to uniform arc-length spacing,
% %                           and reconstructed from only these harmonics --
% %                           fewer harmonics = smoother/rounder outline
% %                           (small anatomical notches may be smoothed away);
% %                           more harmonics = closer to the raw pixel shape.
% %   'boundaryResamplePoints' : number of uniformly arc-length-spaced points
% %                           the raw boundary is resampled to before Fourier
% %                           smoothing, default = 400. Should comfortably
% %                           exceed 2x boundaryNumHarmonics.
% %
% %   'colGapFrac' : gap between adjacent frame columns, as a fraction of the
% %                  figure width, default = 0.004 (very tight). Set
% %                  independently from 'rowGapFrac' (tiles are laid out with
% %                  manually positioned axes rather than tiledlayout, since
% %                  tiledlayout's 'TileSpacing' cannot differ by direction).
% %   'rowGapFrac' : gap between adjacent motif rows, as a fraction of the
% %                  figure height, default = 0.018.
% %
% %   'upsampleFactor' : spatial upsampling factor applied to each 64x64 frame
% %                      before display/printing, default = 4 (i.e. 256x256).
% %                      Interpolation is mask-aware (see Notes) to avoid
% %                      bleeding intensity across the cortex boundary.
% %   'interpMethod'   : interpolation method passed to imresize for the
% %                      upsampling step, default = 'bicubic'.
% %
% % Notes:
% %   - Gaussian smoothing and percentile/global scaling are computed at the
% %     native 64x64 resolution (as before); spatial upsampling happens last,
% %     purely for display/print rendering.
% %   - Background is now white. Pixels outside the dorsal-cortex mask are
% %     blended toward white (with a soft, anti-aliased edge) and baked
% %     directly into an opaque RGB image, rather than relying on continuous
% %     AlphaData transparency -- the 'painters' renderer used for -dpdf
% %     printing does not reliably support partial alpha and can crash on it.
% %   - Per-motif scaling rescales each motif across all lags using the
% %     requested percentile range.
% %   - Printed filenames preserve original motif IDs.
% %
% % Frame-dropping logic:
% %   Given L available frames and a request to keep nFrames <= L, the
% %   function computes, for each frame, an "energy" value (sum of squared
% %   pixel values across the cortex, ignoring NaNs). It then slides a window
% %   of length nFrames across the 1..L sequence and keeps whichever
% %   contiguous window has the greatest total energy. Because the window is
% %   contiguous, its complement (the dropped frames) is always a prefix
% %   and/or suffix of the sequence -- frames are never dropped from the
% %   middle, and frame order is preserved.
% 
% %% ------------------------------------------------------------------------
% %  Save directory
% % -------------------------------------------------------------------------
% 
% saveFigDir = '/Volumes/buschman/Rodent Data/dualImaging_parkj/collectFigure/montage';
% 
% if ispc
%     saveFigDir = 'Z:\Rodent Data\dualImaging_parkj\collectFigure\montage';
% end
% 
% if exist(saveFigDir, 'dir') ~= 7
%     mkdir(saveFigDir);
% end
% 
% %% ------------------------------------------------------------------------
% %  Parse inputs
% % -------------------------------------------------------------------------
% 
% p = inputParser;
% 
% p.addParameter('motifsPerFig', 1, ...
%     @(x) isnumeric(x) && isscalar(x) && x >= 1);
% 
% p.addParameter('printLogic', true, ...
%     @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
% 
% p.addParameter('selectMotifs', [], ...
%     @(x) isempty(x) || isnumeric(x));
% 
% p.addParameter('prctileRange', [1 99.5], ...
%     @(x) isnumeric(x) && numel(x) == 2 && x(1) < x(2));
% 
% p.addParameter('scaleMode', 'perMotif', ...
%     @(x) ischar(x) || isstring(x));
% 
% p.addParameter('displayRange', [0 1], ...
%     @(x) isnumeric(x) && numel(x) == 2 && x(1) < x(2));
% 
% p.addParameter('colormapName', 'magma', ...
%     @(x) ischar(x) || isstring(x) || isnumeric(x));
% 
% p.addParameter('figNamePrefix', 'Motifs', ...
%     @(x) ischar(x) || isstring(x));
% 
% p.addParameter('smoothLogic', false, ...
%     @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
% 
% p.addParameter('gaussianSigma', 1, ...
%     @(x) isnumeric(x) && isscalar(x) && x > 0);
% 
% % New: frame-count control + drop logic
% p.addParameter('nFrames', [], ...
%     @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 1 && mod(x,1) == 0));
% 
% p.addParameter('frameSelectionScope', 'global', ...
%     @(x) ischar(x) || isstring(x));
% 
% % New: cortex boundary overlay
% p.addParameter('cortexBoundaryLogic', true, ...
%     @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
% 
% p.addParameter('boundaryColor', [0.5 0.5 0.5], ...
%     @(x) isnumeric(x) && numel(x) == 3);
% 
% p.addParameter('boundaryLineWidth', 1.5, ...
%     @(x) isnumeric(x) && isscalar(x) && x > 0);
% 
% % Boundary smoothing via Fourier descriptors (see localGetCortexMaskAndBoundary):
% % the raw pixel boundary is traced once, resampled to uniform arc-length
% % spacing, and reconstructed from only its low-frequency components,
% % guaranteeing a smooth, seamlessly closed curve.
% % New: sub-pixel boundary smoothness controls (see localGetCortexMaskAndBoundary)
% p.addParameter('boundaryNumHarmonics', 15, ...
%     @(x) isnumeric(x) && isscalar(x) && x >= 1);
% 
% p.addParameter('boundaryResamplePoints', 400, ...
%     @(x) isnumeric(x) && isscalar(x) && x >= 16);
% 
% % New: independent horizontal/vertical tile spacing (fraction of figure
% % width/height given to the gap between adjacent tiles). Unlike
% % tiledlayout's 'TileSpacing', these are set independently per axis.
% p.addParameter('colGapFrac', 0.004, ...
%     @(x) isnumeric(x) && isscalar(x) && x >= 0 && x < 1);
% 
% p.addParameter('rowGapFrac', 0.018, ...
%     @(x) isnumeric(x) && isscalar(x) && x >= 0 && x < 1);
% 
% % New: spatial upsampling
% p.addParameter('upsampleFactor', 4, ...
%     @(x) isnumeric(x) && isscalar(x) && x >= 1);
% 
% p.addParameter('interpMethod', 'bicubic', ...
%     @(x) ischar(x) || isstring(x));
% 
% p.parse(varargin{:});
% 
% motifsPerFig    = p.Results.motifsPerFig;
% printLogic      = logical(p.Results.printLogic);
% selectMotifs    = p.Results.selectMotifs;
% prctileRange    = p.Results.prctileRange;
% scaleMode       = char(p.Results.scaleMode);
% displayRange    = p.Results.displayRange;
% colormapName    = p.Results.colormapName;
% figNamePrefix   = char(p.Results.figNamePrefix);
% smoothLogic     = logical(p.Results.smoothLogic);
% gaussianSigma   = p.Results.gaussianSigma;
% 
% nFramesRequest      = p.Results.nFrames;
% frameSelectionScope = lower(char(p.Results.frameSelectionScope));
% 
% cortexBoundaryLogic = logical(p.Results.cortexBoundaryLogic);
% boundaryColor       = p.Results.boundaryColor;
% boundaryLineWidth   = p.Results.boundaryLineWidth;
% boundaryNumHarmonics   = p.Results.boundaryNumHarmonics;
% boundaryResamplePoints = p.Results.boundaryResamplePoints;
% colGapFrac          = p.Results.colGapFrac;
% rowGapFrac          = p.Results.rowGapFrac;
% 
% upsampleFactor = p.Results.upsampleFactor;
% interpMethod   = char(p.Results.interpMethod);
% 
% if ~ismember(frameSelectionScope, {'global', 'permotif'})
%     error('frameSelectionScope must be ''global'' or ''perMotif''.');
% end
% 
% %% ------------------------------------------------------------------------
% %  Validate motif tensor
% % -------------------------------------------------------------------------
% 
% if ~isnumeric(motif_w) || ndims(motif_w) ~= 3
%     error('motif_w must be a numeric [P x K x L] array.');
% end
% 
% [P, nMotifTotal, nFramesTotal] = size(motif_w);
% 
% if isempty(selectMotifs)
%     selectMotifs = 1:nMotifTotal;
% else
%     selectMotifs = unique(selectMotifs(:))';
% 
%     if any(selectMotifs < 1) || any(selectMotifs > nMotifTotal)
%         error('selectMotifs contains indices outside valid motif range 1:%d.', nMotifTotal);
%     end
% end
% 
% % Restrict displayed motifs, but preserve original IDs in selectMotifs
% motif_w_sel = motif_w(:, selectMotifs, :);
% 
% nMotifSel = numel(selectMotifs);
% 
% %% ------------------------------------------------------------------------
% %  Cortex mask + boundary (used for transparency and the gray outline)
% % -------------------------------------------------------------------------
% 
% [cortexMask, cortexBoundaryXY] = localGetCortexMaskAndBoundary( ...
%     nanpxs, cortexBoundaryLogic, boundaryNumHarmonics, boundaryResamplePoints);
% 
% if isempty(cortexBoundaryXY)
%     cortexBoundaryLogic = false;
% end
% 
% % Precompute the high-resolution alpha mask ONCE, by rasterizing the exact
% % same smooth boundary curve used for the plotted line (via poly2mask)
% % rather than separately bicubic-upsampling the original blocky binary
% % mask. Deriving the mask and the line from two different smoothing
% % methods (Fourier descriptors for one, raster blurring for the other) has
% % no guarantee of agreeing pixel-for-pixel, which is what kept showing up
% % as a residual fringe/staircase along the edge. Using the identical curve
% % for both eliminates that mismatch entirely. This is also computed once
% % here rather than per motif/frame, since it doesn't depend on the data.
% hiresMaskShared = localRasterizeSmoothMask( ...
%     cortexMask, cortexBoundaryXY, upsampleFactor);
% 
% %% ------------------------------------------------------------------------
% %  Frame selection (drop only from the edges, keep highest-energy window)
% % -------------------------------------------------------------------------
% 
% if isempty(nFramesRequest)
%     nFramesKeep = nFramesTotal;
% else
%     nFramesKeep = nFramesRequest;
% end
% 
% if nFramesKeep > nFramesTotal
%     error('nFrames (%d) cannot exceed the number of available frames (%d).', ...
%         nFramesKeep, nFramesTotal);
% end
% 
% [motif_w_kept, keptFrameLabels] = localSelectFrameWindow( ...
%     motif_w_sel, nFramesKeep, frameSelectionScope);
% 
% nMotifSel_check = size(motif_w_kept, 2); %#ok<NASGU>
% nFramesShow     = size(motif_w_kept, 3);
% 
% nFigures = ceil(nMotifSel / motifsPerFig);
% 
% %% ------------------------------------------------------------------------
% %  Optional global scaling range (computed on the frames actually shown)
% % -------------------------------------------------------------------------
% 
% switch lower(scaleMode)
% 
%     case 'global'
% 
%         vals = motif_w_kept(:);
%         vals = vals(isfinite(vals));
%         vals = vals(vals ~= 0);
% 
%         if isempty(vals)
%             globalLow = 0;
%             globalHigh = 1;
%         else
%             globalLow  = prctile(vals, prctileRange(1));
%             globalHigh = prctile(vals, prctileRange(2));
% 
%             if globalHigh <= globalLow
%                 globalLow  = min(vals);
%                 globalHigh = max(vals);
%             end
% 
%             if globalHigh <= globalLow
%                 globalLow  = 0;
%                 globalHigh = 1;
%             end
%         end
% 
%     case {'permotif', 'none'}
% 
%         globalLow = [];
%         globalHigh = [];
% 
%     otherwise
% 
%         error('scaleMode must be ''perMotif'', ''global'', or ''none''.');
% end
% 
% %% ------------------------------------------------------------------------
% %  Resolve colormap once (used to manually bake RGB below -- see note in
% %  the render loop about why we don't rely on continuous AlphaData)
% % -------------------------------------------------------------------------
% 
% cmap = localResolveColormap(colormapName);
% 
% %% ------------------------------------------------------------------------
% %  Main montage loop
% % -------------------------------------------------------------------------
% 
% for figIdx = 1:nFigures
% 
%     % Indices within selected motif list
%     startIdxLocal = (figIdx - 1) * motifsPerFig + 1;
%     endIdxLocal   = min(figIdx * motifsPerFig, nMotifSel);
%     motifsThisFig = endIdxLocal - startIdxLocal + 1;
% 
%     % Original motif IDs shown in this figure
%     motifIDsThisFig = selectMotifs(startIdxLocal:endIdxLocal);
% 
%     h = figure('Color', 'w');
% 
%     % Layout margins (figure-normalized units) reserved for the title,
%     % per-column frame headers, and per-row motif labels.
%     leftMarginFrac   = 0.06;
%     topMarginFrac    = 0.075;
%     bottomMarginFrac = 0.01;
%     rightMarginFrac  = 0.01;
% 
%     nCols = nFramesShow;
%     nRows = motifsThisFig;
% 
%     tileW = (1 - leftMarginFrac - rightMarginFrac - (nCols - 1) * colGapFrac) / nCols;
%     tileH = (1 - topMarginFrac  - bottomMarginFrac - (nRows - 1) * rowGapFrac) / nRows;
% 
%     % Each tile calls axis(ax,'image','off') to preserve the (square)
%     % image's aspect ratio. If a tile's allocated box isn't itself square,
%     % MATLAB centers the square image within a letterboxed sub-region of
%     % that box -- padding that colGapFrac has no control over, since it
%     % isn't axes spacing at all. Fix this at the source: size the figure's
%     % physical pixels so tileW/tileH normalized fractions correspond to an
%     % actually-square physical tile, and axis('image') never needs to pad.
%     targetTilePx = 220;
%     figWidthPx  = targetTilePx / tileW;
%     figHeightPx = targetTilePx / tileH;
%     set(h, 'Units', 'pixels', 'Position', [100, 100, figWidthPx, figHeightPx]);
% 
%     for i = 1:motifsThisFig
% 
%         motifIdxLocal = startIdxLocal + i - 1;
% 
%         %% ----------------------------------------------------------------
%         %  Reconstruct motif into image stack
%         % -----------------------------------------------------------------
% 
%         if P == 64 * 64
% 
%             motif = reshape(squeeze(motif_w_kept(:, motifIdxLocal, :)), ...
%                 64, 64, []);
% 
%         else
% 
%             motif = conditionDffMat( ...
%                 squeeze(motif_w_kept(:, motifIdxLocal, :))', ...
%                 nanpxs);
%         end
% 
%         %% ----------------------------------------------------------------
%         %  Optional Gaussian smoothing (native resolution, before scaling)
%         % -----------------------------------------------------------------
% 
%         if smoothLogic
%             motif = applyImgaussfilt(motif, 'sigma', gaussianSigma);
%         end
% 
%         %% ----------------------------------------------------------------
%         %  Scaling
%         % -----------------------------------------------------------------
% 
%         switch lower(scaleMode)
% 
%             case 'permotif'
% 
%                 motif = localScaleByPercentile(motif, prctileRange);
% 
%             case 'global'
% 
%                 motif = localScaleByFixedRange(motif, globalLow, globalHigh);
% 
%             case 'none'
% 
%                 % Leave motif as-is.
%         end
% 
%         %% ----------------------------------------------------------------
%         %  Render each frame as its own tile (mask-aware upsampling +
%         %  transparent background + gray cortex boundary overlay)
%         % -----------------------------------------------------------------
% 
%         for fIdx = 1:nFramesShow
% 
%             frameData = motif(:, :, fIdx);
%             frameData(~isfinite(frameData)) = 0;
% 
%             hiresData = localUpsampleFrameData(frameData, upsampleFactor, interpMethod);
%             hiresMask = hiresMaskShared;
% 
%             % Bake the soft (anti-aliased) mask directly into an opaque RGB
%             % image by blending the colormapped data with a white
%             % background, rather than using continuous AlphaData. The
%             % 'painters' renderer (used below for -dpdf printing) does not
%             % reliably support partial/continuous transparency -- passing
%             % it a non-binary AlphaData matrix can crash MATLAB. Producing
%             % a fully opaque true-color image sidesteps that entirely while
%             % still giving the same smooth, anti-aliased edge.
%             rgbTile = localComposeRGBWithWhiteBackground( ...
%                 hiresData, hiresMask, displayRange, cmap);
% 
%             tileX = leftMarginFrac + (fIdx - 1) * (tileW + colGapFrac);
%             tileY = 1 - topMarginFrac - i * tileH - (i - 1) * rowGapFrac;
% 
%             ax = axes('Parent', h, 'Position', [tileX, tileY, tileW, tileH]); %#ok<LAXES>
%             image(ax, rgbTile);
%             axis(ax, 'image', 'off');
%             set(ax, 'Color', 'w');
%             hold(ax, 'on');
% 
%             if cortexBoundaryLogic
%                 for bIdx = 1:numel(cortexBoundaryXY)
%                     bxy = cortexBoundaryXY{bIdx} * upsampleFactor;
%                     % plot() draws point 1 -> N but does not automatically
%                     % connect N back to 1 -- explicitly close the loop by
%                     % repeating the first point at the end, otherwise a gap
%                     % appears at the seam (wherever the raw boundary trace
%                     % happened to start).
%                     bxyClosed = [bxy; bxy(1, :)];
%                     plot(ax, bxyClosed(:,1), bxyClosed(:,2), '-', ...
%                         'Color', boundaryColor, ...
%                         'LineWidth', boundaryLineWidth);
%                 end
%             end
% 
%             hold(ax, 'off');
% 
%             % Column header (original frame index) on the first row only,
%             % and only when meaningful (i.e. a single shared frame axis).
%             if i == 1 && strcmp(frameSelectionScope, 'global')
%                 title(ax, sprintf('f%d', keptFrameLabels.global(fIdx)), ...
%                     'Color', 'k', 'FontSize', 8, 'FontWeight', 'normal');
%             end
% 
%             % Row label (motif ID [+ frame range if perMotif]) on first column
%             if fIdx == 1
%                 if strcmp(frameSelectionScope, 'permotif')
%                     rowLabel = sprintf('M%d (f%d-%d)', ...
%                         motifIDsThisFig(i), ...
%                         keptFrameLabels.perMotif(motifIdxLocal, 1), ...
%                         keptFrameLabels.perMotif(motifIdxLocal, end));
%                 else
%                     rowLabel = sprintf('M%d', motifIDsThisFig(i));
%                 end
%                 ylabel(ax, rowLabel, 'Color', 'k', 'FontSize', 8, ...
%                     'Rotation', 0, 'HorizontalAlignment', 'right', ...
%                     'VerticalAlignment', 'middle', 'Visible', 'on');
%                 % axis('off') hides the ylabel too, so force it visible
%                 ax.YLabel.Visible = 'on';
%             end
%         end
%     end
% 
%     if smoothLogic
%         smoothLabel = sprintf(' | Gaussian \\sigma = %.2g', gaussianSigma);
%     else
%         smoothLabel = '';
%     end
% 
%     if strcmp(frameSelectionScope, 'global') && nFramesShow < nFramesTotal
%         frameLabel = sprintf(' | frames %d-%d of %d', ...
%             keptFrameLabels.global(1), keptFrameLabels.global(end), nFramesTotal);
%     elseif strcmp(frameSelectionScope, 'permotif') && nFramesShow < nFramesTotal
%         frameLabel = sprintf(' | %d/%d frames (per-motif window)', ...
%             nFramesShow, nFramesTotal);
%     else
%         frameLabel = '';
%     end
% 
%     sgtitle(h, sprintf('Motifs %s%s%s', ...
%         compressMotifIDs(motifIDsThisFig), frameLabel, smoothLabel), ...
%         'Color', 'k', 'Interpreter', 'tex');
% 
%     %% --------------------------------------------------------------------
%     %  Print montage
%     % ---------------------------------------------------------------------
% 
%     if printLogic
% 
%         timestampStr  = datestr(now, 'mmddyy_HHMMSS');
%         motifLabelStr = compressMotifIDs(motifIDsThisFig);
% 
%         if smoothLogic
%             smoothFileStr = sprintf('_gaussSigma%.2g', gaussianSigma);
%             smoothFileStr = strrep(smoothFileStr, '.', 'p');
%         else
%             smoothFileStr = '';
%         end
% 
%         figSaveName = sprintf('%s_%s%s_%s', ...
%             figNamePrefix, ...
%             motifLabelStr, ...
%             smoothFileStr, ...
%             timestampStr);
% 
%         set(h, 'InvertHardcopy', 'off');  % preserve exact on-screen colors
% 
%         % Each tile is now a full true-color RGB bitmap (see
%         % localComposeRGBWithWhiteBackground), not a small indexed image --
%         % 'painters' is a pure vector renderer and is not well-suited to
%         % embedding large raster content in a PDF; it has been a source of
%         % crashes/instability for exactly this kind of mixed raster+vector
%         % (image tiles + boundary lines) figure. 'opengl' handles this
%         % mixed content far more robustly.
%         print(h, fullfile(saveFigDir, figSaveName), ...
%             '-opengl', '-bestfit', '-dpdf');
%     end
% end
% 
% end
% 
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function motifScaled = localScaleByPercentile(motif, prctileRange)
% 
% vals = motif(:);
% vals = vals(isfinite(vals));
% vals = vals(vals ~= 0);
% 
% if isempty(vals)
%     motifScaled = zeros(size(motif), 'like', motif);
%     return
% end
% 
% lo = prctile(vals, prctileRange(1));
% hi = prctile(vals, prctileRange(2));
% 
% if hi <= lo
%     lo = min(vals);
%     hi = max(vals);
% end
% 
% if hi <= lo
%     motifScaled = zeros(size(motif), 'like', motif);
%     return
% end
% 
% motifScaled = (motif - lo) ./ (hi - lo);
% motifScaled(motifScaled < 0) = 0;
% motifScaled(motifScaled > 1) = 1;
% 
% end
% 
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function motifScaled = localScaleByFixedRange(motif, lo, hi)
% 
% if hi <= lo
%     motifScaled = zeros(size(motif), 'like', motif);
%     return
% end
% 
% motifScaled = (motif - lo) ./ (hi - lo);
% motifScaled(motifScaled < 0) = 0;
% motifScaled(motifScaled > 1) = 1;
% 
% end
% 
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function outStr = compressMotifIDs(ids)
% % compressMotifIDs  Convert motif ID vector into compact range string.
% %
% % Example:
% %   [1 2 3 4 5 6 7 9 12 13 14 15] -> '1-7_9_12-15'
% 
% ids = unique(ids(:))';
% 
% if isempty(ids)
%     outStr = '';
%     return;
% end
% 
% rangeParts = {};
% rangeStart = ids(1);
% prevVal = ids(1);
% 
% for ii = 2:numel(ids)
% 
%     if ids(ii) == prevVal + 1
% 
%         prevVal = ids(ii);
% 
%     else
% 
%         rangeParts{end+1} = localRangeToStr(rangeStart, prevVal); %#ok<AGROW>
%         rangeStart = ids(ii);
%         prevVal = ids(ii);
%     end
% end
% 
% rangeParts{end+1} = localRangeToStr(rangeStart, prevVal);
% outStr = strjoin(rangeParts, '_');
% 
% end
% 
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function s = localRangeToStr(a, b)
% 
% if a == b
%     s = sprintf('%d', a);
% else
%     s = sprintf('%d-%d', a, b);
% end
% 
% end
% 
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function cmap = localResolveColormap(colormapName)
% % localResolveColormap
% %
% % Resolves a colormap from either:
% %   - a string/char name, e.g. 'magma', 'parula', 'hot', 'turbo'
% %   - an explicit N x 3 colormap matrix
% 
% if isnumeric(colormapName)
%     cmap = colormapName;
%     return
% end
% 
% cmapName = char(colormapName);
% 
% try
%     cmap = feval(cmapName, 256);
% catch
%     try
%         cmap = eval(cmapName); %#ok<EVLDIR>
%     catch
%         warning('Could not resolve colormap "%s". Falling back to parula.', cmapName);
%         cmap = parula(256);
%     end
% end
% 
% end
% 
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function rgbImage = localComposeRGBWithWhiteBackground(data, alphaMask, displayRange, cmap)
% % localComposeRGBWithWhiteBackground
% %
% % Manually maps `data` through `cmap` (using displayRange as the color
% % axis limits, matching what imagesc/clim would normally do) and blends
% % the result with a white background using alphaMask (a continuous [0,1]
% % field), producing a fully opaque H x W x 3 RGB image.
% %
% % This exists specifically to AVOID passing a continuous (non-binary)
% % AlphaData matrix to imagesc: the 'painters' renderer (used elsewhere in
% % this function for -dpdf printing) does not reliably support partial
% % transparency, and doing so can crash MATLAB. Baking the blend into plain
% % RGB values sidesteps alpha compositing entirely while still producing
% % the same smooth, anti-aliased edge.
% 
% lo = displayRange(1);
% hi = displayRange(2);
% 
% normVal = (data - lo) ./ (hi - lo);
% normVal = min(max(normVal, 0), 1);
% 
% nColors = size(cmap, 1);
% colorIdx = round(normVal * (nColors - 1)) + 1;
% colorIdx = min(max(colorIdx, 1), nColors);
% 
% R = reshape(cmap(colorIdx(:), 1), size(data));
% G = reshape(cmap(colorIdx(:), 2), size(data));
% B = reshape(cmap(colorIdx(:), 3), size(data));
% 
% a = alphaMask;
% 
% Rout = R .* a + 1 .* (1 - a);
% Gout = G .* a + 1 .* (1 - a);
% Bout = B .* a + 1 .* (1 - a);
% 
% rgbImage = cat(3, Rout, Gout, Bout);
% 
% end
% 
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [cortexMask, boundaryXY] = localGetCortexMaskAndBoundary( ...
%     nanpxs, wantBoundary, numHarmonics, resampleN)
% % localGetCortexMaskAndBoundary
% %
% % Derives a 64x64 logical dorsal-cortex mask directly from nanpxs -- the
% % linear indices (or logical mask) of the non-cortex ("NaN") pixels within
% % the full 64x64 = 4096 pixel grid -- and its boundary as a cell array of
% % [x y] coordinate lists (native 64x64 pixel units, one cell per connected
% % boundary component), smoothed via truncated Fourier descriptors.
% %
% % IMPORTANT: whether a real cortex mask exists is entirely a function of
% % nanpxs, NOT of P (the first dimension of motif_w). motif_w can be stored
% % either as [nValidPixels x K x L] (reconstructed via conditionDffMat) or
% % as an already-reshaped [4096 x K x L] full-grid array -- either way, if
% % nanpxs is supplied it still marks the true dorsal-cortex boundary within
% % that 64x64 grid and should be used.
% %
% % If nanpxs is empty, there is no mask information available and the whole
% % 64x64 frame is treated as valid (boundaryXY is returned empty).
% %
% % Boundary smoothness: rather than blurring a rasterized mask (which only
% % ever anti-aliases individual pixel-step edges and struggles to remove the
% % macroscopic staircase left by the native 64x64 resolution), the actual
% % pixel boundary is traced once with bwboundaries, resampled to
% % resampleN uniformly arc-length-spaced points (needed for a well-posed
% % Fourier truncation), represented as a complex sequence z = x + 1i*y, and
% % reconstructed from only its lowest numHarmonics frequency components via
% % FFT/IFFT. Truncating high frequencies of a closed, periodic curve
% % guarantees a result that is both smooth AND exactly seamlessly closed --
% % there is no "sigma in the wrong units" pitfall here, since the smoothing
% % is applied directly to the curve's shape, independent of any raster
% % resolution or supersampling choice.
% 
% boundaryXY = {};
% 
% if isempty(nanpxs)
%     cortexMask = true(64, 64);
%     return
% end
% 
% nanFlagVec = false(64 * 64, 1);
% 
% if islogical(nanpxs)
%     nanFlagVec(:) = nanpxs(:);
% else
%     nanFlagVec(nanpxs(:)) = true;
% end
% 
% cortexMask = reshape(~nanFlagVec, 64, 64);
% 
% if ~wantBoundary
%     return
% end
% 
% if exist('bwboundaries', 'file') ~= 2
%     warning('bwboundaries (Image Processing Toolbox) not found; skipping cortex boundary overlay.');
%     return
% end
% 
% rawBoundaries = bwboundaries(cortexMask, 'noholes');
% boundaryXY = cell(size(rawBoundaries));
% 
% for k = 1:numel(rawBoundaries)
%     % bwboundaries returns [row col] = [y x]; convert to [x y]
%     rawXY = [rawBoundaries{k}(:,2), rawBoundaries{k}(:,1)];
% 
%     resampledXY = localResampleClosedCurve(rawXY, resampleN);
%     boundaryXY{k} = localSmoothClosedContourFourier(resampledXY, numHarmonics);
% end
% 
% end
% 
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function resampledXY = localResampleClosedCurve(xy, nResample)
% % localResampleClosedCurve
% %
% % Resamples a closed 2D polygon (M x 2, [x y]) to nResample points spaced
% % uniformly by arc length around the loop. Uniform spacing is required for
% % a clean/well-posed Fourier-descriptor truncation afterward.
% 
% % Ensure explicitly closed (first point repeated at the end) for arc-length
% % accumulation, then drop the duplicate after resampling.
% if norm(xy(1,:) - xy(end,:)) > 1e-9
%     xy = [xy; xy(1,:)];
% end
% 
% % Guard against zero-length (duplicate/repeated) points -- bwboundaries can
% % occasionally emit a repeated point where the traced path touches a thin
% % or degenerate pixel connection. A zero-length segment here would give
% % interp1 a non-strictly-increasing cumulative-distance vector, which can
% % produce a small localized artifact in the curve after Fourier smoothing.
% segLenRaw = sqrt(sum(diff(xy).^2, 2));
% keepPoint = [true; segLenRaw > 1e-9];
% xy = xy(keepPoint, :);
% 
% segLen  = sqrt(sum(diff(xy).^2, 2));
% cumDist = [0; cumsum(segLen)];
% totalLen = cumDist(end);
% 
% if totalLen == 0
%     resampledXY = repmat(xy(1,:), nResample, 1);
%     return
% end
% 
% targetDist = linspace(0, totalLen, nResample + 1);
% targetDist(end) = []; % drop duplicate closing point
% 
% xResampled = interp1(cumDist, xy(:,1), targetDist, 'linear');
% yResampled = interp1(cumDist, xy(:,2), targetDist, 'linear');
% 
% resampledXY = [xResampled(:), yResampled(:)];
% 
% end
% 
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function smoothXY = localSmoothClosedContourFourier(xy, numHarmonics)
% % localSmoothClosedContourFourier
% %
% % Smooths a closed, uniformly arc-length-resampled 2D curve (N x 2, [x y])
% % by representing it as a complex sequence z = x + 1i*y, truncating its
% % discrete Fourier transform to the lowest numHarmonics positive and
% % negative frequencies (plus the DC term), and inverting. This is the
% % classic "Fourier descriptor" smoothing approach: it is guaranteed to
% % return an exactly closed curve (the representation is inherently
% % periodic) that is smooth by construction, with numHarmonics as the sole,
% % resolution-independent smoothness knob (fewer harmonics = smoother/
% % rounder; more harmonics = closer to the original traced shape).
% 
% N = size(xy, 1);
% z = complex(xy(:,1), xy(:,2));
% 
% Z = fft(z);
% 
% numHarmonics = min(numHarmonics, floor((N - 1) / 2));
% 
% keepMask = false(N, 1);
% keepMask(1) = true; % DC term
% keepMask(2:(numHarmonics + 1)) = true;               % low positive frequencies
% keepMask((N - numHarmonics + 1):N) = true;           % mirrored negative frequencies
% 
% Z(~keepMask) = 0;
% 
% zSmooth = ifft(Z);
% 
% smoothXY = [real(zSmooth), imag(zSmooth)];
% 
% end
% 
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function hiresData = localUpsampleFrameData(frameData, upsampleFactor, interpMethod)
% % localUpsampleFrameData
% %
% % Spatially upsamples a single 64x64 frame's data for high-resolution
% % display. The data is interpolated as-is, with no masking applied here --
% % masking/fading to white is handled entirely by the shared alpha mask
% % (see localRasterizeSmoothMask), computed once from the same smooth
% % boundary curve used for the plotted line. Zeroing data outside the mask
% % before interpolation was tried and rejected: combined with a separately
% % smoothed alpha, it causes visible color fringing (a dark tinge from the
% % hard zero shows through wherever alpha is only partially faded).
% 
% if upsampleFactor == 1
%     hiresData = frameData;
%     return
% end
% 
% targetSize = size(frameData) * upsampleFactor;
% hiresData = imresize(frameData, targetSize, interpMethod);
% 
% end
% 
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function hiresMask = localRasterizeSmoothMask(cortexMask, boundaryXY, upsampleFactor)
% % localRasterizeSmoothMask
% %
% % Builds the high-resolution alpha mask by rasterizing the SAME smooth
% % boundary curve(s) used for the plotted gray line (via poly2mask), rather
% % than independently bicubic-upsampling the original native-resolution
% % binary mask. Using two different smoothing methods for the line (Fourier
% % descriptors) and the mask (raster blurring) gives no guarantee they agree
% % pixel-for-pixel, which produced a visible fringe/staircase along the
% % edge in earlier attempts. Deriving both from the identical curve
% % eliminates that mismatch structurally. A light final Gaussian pass
% % anti-aliases the (now much finer, supersampled-grid-scale) rasterization
% % step, which is a much smaller and less objectionable artifact than the
% % native-pixel-grid blockiness this replaces.
% %
% % Falls back to a plain bicubic-upsampled binary mask if no boundary curve
% % is available (e.g. cortexBoundaryLogic was false, or nanpxs was empty).
% 
% targetSize = size(cortexMask) * upsampleFactor;
% 
% if isempty(boundaryXY) || exist('poly2mask', 'file') ~= 2
%     hiresMask = imresize(double(cortexMask), targetSize, 'bicubic');
%     hiresMask = min(max(hiresMask, 0), 1);
%     return
% end
% 
% maskAccum = false(targetSize);
% 
% for k = 1:numel(boundaryXY)
%     xy = boundaryXY{k} * upsampleFactor;
%     bw = poly2mask(xy(:,1), xy(:,2), targetSize(1), targetSize(2));
%     maskAccum = maskAccum | bw;
% end
% 
% hiresMask = double(maskAccum);
% 
% if exist('imgaussfilt', 'file') == 2
%     hiresMask = imgaussfilt(hiresMask, 1.0); % mild anti-aliasing at supersampled-grid scale
% end
% 
% hiresMask = min(max(hiresMask, 0), 1);
% 
% end
% 
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [motif_w_kept, keptFrameLabels] = localSelectFrameWindow(motif_w_sel, nFramesKeep, scope)
% % localSelectFrameWindow
% %
% % Selects, for each motif (or globally across all motifs), the contiguous
% % window of nFramesKeep frames with the highest total energy, out of the
% % nFramesTotal frames available. Energy for a frame is the sum of squared
% % pixel values across all pixels/motifs being pooled, ignoring NaNs.
% %
% % Because the chosen window is always contiguous, its complement (the
% % dropped frames) is always a prefix and/or suffix of the sequence -- i.e.
% % frames are dropped only from the edges, never from the middle, and frame
% % order is preserved.
% %
% % Outputs:
% %   motif_w_kept    : [P x nMotifSel x nFramesKeep] tensor with the
% %                     dropped frames removed.
% %   keptFrameLabels : struct with fields:
% %                       .global   - 1 x nFramesKeep vector of original frame
% %                                   indices (valid/used when scope=='global')
% %                       .perMotif - nMotifSel x nFramesKeep matrix of
% %                                   original frame indices per motif (valid/
% %                                   used when scope=='perMotif')
% 
% [P, nMotifSel, nFramesTotal] = size(motif_w_sel);
% 
% keptFrameLabels = struct('global', [], 'perMotif', []);
% 
% if nFramesKeep == nFramesTotal
%     motif_w_kept = motif_w_sel;
%     keptFrameLabels.global   = 1:nFramesTotal;
%     keptFrameLabels.perMotif = repmat(1:nFramesTotal, nMotifSel, 1);
%     return
% end
% 
% % Energy per motif per frame: [nMotifSel x nFramesTotal]
% sq = motif_w_sel .^ 2;
% sq(~isfinite(sq)) = 0;
% energyMotifFrame = squeeze(sum(sq, 1));       % nMotifSel x nFramesTotal
% if nMotifSel == 1
%     energyMotifFrame = reshape(energyMotifFrame, 1, nFramesTotal);
% end
% 
% switch scope
% 
%     case 'global'
% 
%         aggregateEnergy = sum(energyMotifFrame, 1);         % 1 x nFramesTotal
%         winStart = localBestWindowStart(aggregateEnergy, nFramesKeep);
%         keepIdx  = winStart:(winStart + nFramesKeep - 1);
% 
%         motif_w_kept = motif_w_sel(:, :, keepIdx);
%         keptFrameLabels.global   = keepIdx;
%         keptFrameLabels.perMotif = repmat(keepIdx, nMotifSel, 1);
% 
%     case 'permotif'
% 
%         motif_w_kept = zeros(P, nMotifSel, nFramesKeep, 'like', motif_w_sel);
%         perMotifIdx  = zeros(nMotifSel, nFramesKeep);
% 
%         for m = 1:nMotifSel
%             winStart = localBestWindowStart(energyMotifFrame(m, :), nFramesKeep);
%             keepIdx  = winStart:(winStart + nFramesKeep - 1);
% 
%             motif_w_kept(:, m, :) = motif_w_sel(:, m, keepIdx);
%             perMotifIdx(m, :) = keepIdx;
%         end
% 
%         keptFrameLabels.perMotif = perMotifIdx;
% 
%     otherwise
% 
%         error('Unknown frameSelectionScope: %s', scope);
% end
% 
% end
% 
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function winStart = localBestWindowStart(energyVec, winLen)
% % localBestWindowStart
% %
% % Slides a window of length winLen across energyVec (1 x L) and returns the
% % start index of the window with the maximum total energy. Ties are broken
% % in favor of the earliest (smallest-index) window.
% 
% L = numel(energyVec);
% nWindows = L - winLen + 1;
% 
% winSums = zeros(1, nWindows);
% for s = 1:nWindows
%     winSums(s) = sum(energyVec(s:(s + winLen - 1)));
% end
% 
% [~, winStart] = max(winSums);
% 
% end
% 
% 
% 
% % function montageMotifsPrintAdvanced(motif_w, nanpxs, varargin)
% % % MontageMotifsPrintAdvanced  Display and optionally print motif montages.
% % %
% % % Usage:
% % %   MontageMotifsPrintAdvanced(motif_w, nanpxs)
% % %   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'motifsPerFig', 3)
% % %   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'printLogic', false)
% % %   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'selectMotifs', [1 3 6 7 9 12 13 14 15])
% % %   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'prctileRange', [1 99.5])
% % %   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'smoothLogic', true, 'gaussianSigma', 1)
% % %   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'nFrames', 7)
% % %   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'nFrames', 7, 'frameSelectionScope', 'perMotif')
% % %   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'upsampleFactor', 6, 'interpMethod', 'bicubic')
% % %
% % % Inputs:
% % %   motif_w      : [P x K x L] array of spatiotemporal motifs
% % %                  pixels x motifs x frames/lags
% % %   nanpxs       : NaN pixel information used by conditionDffMat. Also used
% % %                  here to derive the dorsal-cortex mask/boundary.
% % %
% % % Name-value pairs:
% % %   'motifsPerFig'  : number of motifs per figure, default = 1
% % %   'printLogic'    : logical scalar, whether to print to PDF, default = true
% % %   'selectMotifs'  : vector of original motif IDs to display, default = all
% % %   'prctileRange'  : percentile range for per-motif scaling, default = [1 99.5]
% % %   'scaleMode'     : 'perMotif', 'global', or 'none', default = 'perMotif'
% % %   'displayRange'  : display range after scaling, default = [0 1]
% % %   'colormapName'  : colormap name or N x 3 matrix, default = 'magma'
% % %   'figNamePrefix' : output figure filename prefix, default = 'Motifs'
% % %   'smoothLogic'   : apply Gaussian smoothing frame-by-frame, default = false
% % %   'gaussianSigma' : sigma for Gaussian smoothing, default = 1
% % %
% % %   'nFrames'             : number of frames to display per motif, out of the
% % %                           L available. Default = [] (show all L frames).
% % %                           When nFrames < L, the dropped frames are chosen
% % %                           by an energy-based contiguous-window search (see
% % %                           "Frame-dropping logic" below) -- frames are never
% % %                           dropped out of the middle of the sequence.
% % %   'frameSelectionScope' : 'global' (default) or 'perMotif'.
% % %                           'global'   - a single common window of frames is
% % %                                        chosen using the pooled (summed)
% % %                                        energy across all displayed motifs,
% % %                                        so every row/tile in the montage
% % %                                        shares the same underlying frame
% % %                                        indices (recommended: keeps a
% % %                                        common time axis across motifs).
% % %                           'perMotif' - each motif independently keeps its
% % %                                        own best window. Frame indices may
% % %                                        then differ row-to-row; the kept
% % %                                        frame range is annotated on each row.
% % %
% % %   'cortexBoundaryLogic' : draw the dorsal-cortex boundary in gray,
% % %                           default = true (only applies when nanpxs actually
% % %                           crops the image, i.e. P ~= 64*64).
% % %   'boundaryColor'        : RGB triplet for the boundary line, default = [0.5 0.5 0.5]
% % %   'boundaryLineWidth'    : line width for the boundary, default = 1.5
% % %   'boundaryNumHarmonics' : number of low-frequency Fourier harmonics kept
% % %                           when smoothing the cortex boundary curve,
% % %                           default = 15. The raw pixel boundary is traced
% % %                           once, resampled to uniform arc-length spacing,
% % %                           and reconstructed from only these harmonics --
% % %                           fewer harmonics = smoother/rounder outline
% % %                           (small anatomical notches may be smoothed away);
% % %                           more harmonics = closer to the raw pixel shape.
% % %   'boundaryResamplePoints' : number of uniformly arc-length-spaced points
% % %                           the raw boundary is resampled to before Fourier
% % %                           smoothing, default = 400. Should comfortably
% % %                           exceed 2x boundaryNumHarmonics.
% % %
% % %   'colGapFrac' : gap between adjacent frame columns, as a fraction of the
% % %                  figure width, default = 0.004 (very tight). Set
% % %                  independently from 'rowGapFrac' (tiles are laid out with
% % %                  manually positioned axes rather than tiledlayout, since
% % %                  tiledlayout's 'TileSpacing' cannot differ by direction).
% % %   'rowGapFrac' : gap between adjacent motif rows, as a fraction of the
% % %                  figure height, default = 0.018.
% % %
% % %   'upsampleFactor' : spatial upsampling factor applied to each 64x64 frame
% % %                      before display/printing, default = 4 (i.e. 256x256).
% % %                      Interpolation is mask-aware (see Notes) to avoid
% % %                      bleeding intensity across the cortex boundary.
% % %   'interpMethod'   : interpolation method passed to imresize for the
% % %                      upsampling step, default = 'bicubic'.
% % %
% % % Notes:
% % %   - Gaussian smoothing and percentile/global scaling are computed at the
% % %     native 64x64 resolution (as before); spatial upsampling happens last,
% % %     purely for display/print rendering.
% % %   - Background is now white. Pixels outside the dorsal-cortex mask are
% % %     blended toward white (with a soft, anti-aliased edge) and baked
% % %     directly into an opaque RGB image, rather than relying on continuous
% % %     AlphaData transparency -- the 'painters' renderer used for -dpdf
% % %     printing does not reliably support partial alpha and can crash on it.
% % %   - Per-motif scaling rescales each motif across all lags using the
% % %     requested percentile range.
% % %   - Printed filenames preserve original motif IDs.
% % %
% % % Frame-dropping logic:
% % %   Given L available frames and a request to keep nFrames <= L, the
% % %   function computes, for each frame, an "energy" value (sum of squared
% % %   pixel values across the cortex, ignoring NaNs). It then slides a window
% % %   of length nFrames across the 1..L sequence and keeps whichever
% % %   contiguous window has the greatest total energy. Because the window is
% % %   contiguous, its complement (the dropped frames) is always a prefix
% % %   and/or suffix of the sequence -- frames are never dropped from the
% % %   middle, and frame order is preserved.
% % 
% % %% ------------------------------------------------------------------------
% % %  Save directory
% % % -------------------------------------------------------------------------
% % 
% % saveFigDir = '/Volumes/buschman/Rodent Data/dualImaging_parkj/collectFigure/montage';
% % 
% % if ispc
% %     saveFigDir = 'Z:\Rodent Data\dualImaging_parkj\collectFigure\montage';
% % end
% % 
% % if exist(saveFigDir, 'dir') ~= 7
% %     mkdir(saveFigDir);
% % end
% % 
% % %% ------------------------------------------------------------------------
% % %  Parse inputs
% % % -------------------------------------------------------------------------
% % 
% % p = inputParser;
% % 
% % p.addParameter('motifsPerFig', 1, ...
% %     @(x) isnumeric(x) && isscalar(x) && x >= 1);
% % 
% % p.addParameter('printLogic', true, ...
% %     @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
% % 
% % p.addParameter('selectMotifs', [], ...
% %     @(x) isempty(x) || isnumeric(x));
% % 
% % p.addParameter('prctileRange', [1 99.5], ...
% %     @(x) isnumeric(x) && numel(x) == 2 && x(1) < x(2));
% % 
% % p.addParameter('scaleMode', 'perMotif', ...
% %     @(x) ischar(x) || isstring(x));
% % 
% % p.addParameter('displayRange', [0 1], ...
% %     @(x) isnumeric(x) && numel(x) == 2 && x(1) < x(2));
% % 
% % p.addParameter('colormapName', 'magma', ...
% %     @(x) ischar(x) || isstring(x) || isnumeric(x));
% % 
% % p.addParameter('figNamePrefix', 'Motifs', ...
% %     @(x) ischar(x) || isstring(x));
% % 
% % p.addParameter('smoothLogic', false, ...
% %     @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
% % 
% % p.addParameter('gaussianSigma', 1, ...
% %     @(x) isnumeric(x) && isscalar(x) && x > 0);
% % 
% % % New: frame-count control + drop logic
% % p.addParameter('nFrames', [], ...
% %     @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 1 && mod(x,1) == 0));
% % 
% % p.addParameter('frameSelectionScope', 'global', ...
% %     @(x) ischar(x) || isstring(x));
% % 
% % % New: cortex boundary overlay
% % p.addParameter('cortexBoundaryLogic', true, ...
% %     @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
% % 
% % p.addParameter('boundaryColor', [0.5 0.5 0.5], ...
% %     @(x) isnumeric(x) && numel(x) == 3);
% % 
% % p.addParameter('boundaryLineWidth', 1.5, ...
% %     @(x) isnumeric(x) && isscalar(x) && x > 0);
% % 
% % % Boundary smoothing via Fourier descriptors (see localGetCortexMaskAndBoundary):
% % % the raw pixel boundary is traced once, resampled to uniform arc-length
% % % spacing, and reconstructed from only its low-frequency components,
% % % guaranteeing a smooth, seamlessly closed curve.
% % % New: sub-pixel boundary smoothness controls (see localGetCortexMaskAndBoundary)
% % p.addParameter('boundaryNumHarmonics', 15, ...
% %     @(x) isnumeric(x) && isscalar(x) && x >= 1);
% % 
% % p.addParameter('boundaryResamplePoints', 400, ...
% %     @(x) isnumeric(x) && isscalar(x) && x >= 16);
% % 
% % % New: independent horizontal/vertical tile spacing (fraction of figure
% % % width/height given to the gap between adjacent tiles). Unlike
% % % tiledlayout's 'TileSpacing', these are set independently per axis.
% % p.addParameter('colGapFrac', 0.004, ...
% %     @(x) isnumeric(x) && isscalar(x) && x >= 0 && x < 1);
% % 
% % p.addParameter('rowGapFrac', 0.018, ...
% %     @(x) isnumeric(x) && isscalar(x) && x >= 0 && x < 1);
% % 
% % % New: spatial upsampling
% % p.addParameter('upsampleFactor', 4, ...
% %     @(x) isnumeric(x) && isscalar(x) && x >= 1);
% % 
% % p.addParameter('interpMethod', 'bicubic', ...
% %     @(x) ischar(x) || isstring(x));
% % 
% % p.parse(varargin{:});
% % 
% % motifsPerFig    = p.Results.motifsPerFig;
% % printLogic      = logical(p.Results.printLogic);
% % selectMotifs    = p.Results.selectMotifs;
% % prctileRange    = p.Results.prctileRange;
% % scaleMode       = char(p.Results.scaleMode);
% % displayRange    = p.Results.displayRange;
% % colormapName    = p.Results.colormapName;
% % figNamePrefix   = char(p.Results.figNamePrefix);
% % smoothLogic     = logical(p.Results.smoothLogic);
% % gaussianSigma   = p.Results.gaussianSigma;
% % 
% % nFramesRequest      = p.Results.nFrames;
% % frameSelectionScope = lower(char(p.Results.frameSelectionScope));
% % 
% % cortexBoundaryLogic = logical(p.Results.cortexBoundaryLogic);
% % boundaryColor       = p.Results.boundaryColor;
% % boundaryLineWidth   = p.Results.boundaryLineWidth;
% % boundaryNumHarmonics   = p.Results.boundaryNumHarmonics;
% % boundaryResamplePoints = p.Results.boundaryResamplePoints;
% % colGapFrac          = p.Results.colGapFrac;
% % rowGapFrac          = p.Results.rowGapFrac;
% % 
% % upsampleFactor = p.Results.upsampleFactor;
% % interpMethod   = char(p.Results.interpMethod);
% % 
% % if ~ismember(frameSelectionScope, {'global', 'permotif'})
% %     error('frameSelectionScope must be ''global'' or ''perMotif''.');
% % end
% % 
% % %% ------------------------------------------------------------------------
% % %  Validate motif tensor
% % % -------------------------------------------------------------------------
% % 
% % if ~isnumeric(motif_w) || ndims(motif_w) ~= 3
% %     error('motif_w must be a numeric [P x K x L] array.');
% % end
% % 
% % [P, nMotifTotal, nFramesTotal] = size(motif_w);
% % 
% % if isempty(selectMotifs)
% %     selectMotifs = 1:nMotifTotal;
% % else
% %     selectMotifs = unique(selectMotifs(:))';
% % 
% %     if any(selectMotifs < 1) || any(selectMotifs > nMotifTotal)
% %         error('selectMotifs contains indices outside valid motif range 1:%d.', nMotifTotal);
% %     end
% % end
% % 
% % % Restrict displayed motifs, but preserve original IDs in selectMotifs
% % motif_w_sel = motif_w(:, selectMotifs, :);
% % 
% % nMotifSel = numel(selectMotifs);
% % 
% % %% ------------------------------------------------------------------------
% % %  Cortex mask + boundary (used for transparency and the gray outline)
% % % -------------------------------------------------------------------------
% % 
% % [cortexMask, cortexBoundaryXY] = localGetCortexMaskAndBoundary( ...
% %     nanpxs, cortexBoundaryLogic, boundaryNumHarmonics, boundaryResamplePoints);
% % 
% % if isempty(cortexBoundaryXY)
% %     cortexBoundaryLogic = false;
% % end
% % 
% % % Precompute the high-resolution alpha mask ONCE, by rasterizing the exact
% % % same smooth boundary curve used for the plotted line (via poly2mask)
% % % rather than separately bicubic-upsampling the original blocky binary
% % % mask. Deriving the mask and the line from two different smoothing
% % % methods (Fourier descriptors for one, raster blurring for the other) has
% % % no guarantee of agreeing pixel-for-pixel, which is what kept showing up
% % % as a residual fringe/staircase along the edge. Using the identical curve
% % % for both eliminates that mismatch entirely. This is also computed once
% % % here rather than per motif/frame, since it doesn't depend on the data.
% % hiresMaskShared = localRasterizeSmoothMask( ...
% %     cortexMask, cortexBoundaryXY, upsampleFactor);
% % 
% % %% ------------------------------------------------------------------------
% % %  Frame selection (drop only from the edges, keep highest-energy window)
% % % -------------------------------------------------------------------------
% % 
% % if isempty(nFramesRequest)
% %     nFramesKeep = nFramesTotal;
% % else
% %     nFramesKeep = nFramesRequest;
% % end
% % 
% % if nFramesKeep > nFramesTotal
% %     error('nFrames (%d) cannot exceed the number of available frames (%d).', ...
% %         nFramesKeep, nFramesTotal);
% % end
% % 
% % [motif_w_kept, keptFrameLabels] = localSelectFrameWindow( ...
% %     motif_w_sel, nFramesKeep, frameSelectionScope);
% % 
% % nMotifSel_check = size(motif_w_kept, 2); %#ok<NASGU>
% % nFramesShow     = size(motif_w_kept, 3);
% % 
% % nFigures = ceil(nMotifSel / motifsPerFig);
% % 
% % %% ------------------------------------------------------------------------
% % %  Optional global scaling range (computed on the frames actually shown)
% % % -------------------------------------------------------------------------
% % 
% % switch lower(scaleMode)
% % 
% %     case 'global'
% % 
% %         vals = motif_w_kept(:);
% %         vals = vals(isfinite(vals));
% %         vals = vals(vals ~= 0);
% % 
% %         if isempty(vals)
% %             globalLow = 0;
% %             globalHigh = 1;
% %         else
% %             globalLow  = prctile(vals, prctileRange(1));
% %             globalHigh = prctile(vals, prctileRange(2));
% % 
% %             if globalHigh <= globalLow
% %                 globalLow  = min(vals);
% %                 globalHigh = max(vals);
% %             end
% % 
% %             if globalHigh <= globalLow
% %                 globalLow  = 0;
% %                 globalHigh = 1;
% %             end
% %         end
% % 
% %     case {'permotif', 'none'}
% % 
% %         globalLow = [];
% %         globalHigh = [];
% % 
% %     otherwise
% % 
% %         error('scaleMode must be ''perMotif'', ''global'', or ''none''.');
% % end
% % 
% % %% ------------------------------------------------------------------------
% % %  Resolve colormap once (used to manually bake RGB below -- see note in
% % %  the render loop about why we don't rely on continuous AlphaData)
% % % -------------------------------------------------------------------------
% % 
% % cmap = localResolveColormap(colormapName);
% % 
% % %% ------------------------------------------------------------------------
% % %  Main montage loop
% % % -------------------------------------------------------------------------
% % 
% % for figIdx = 1:nFigures
% % 
% %     % Indices within selected motif list
% %     startIdxLocal = (figIdx - 1) * motifsPerFig + 1;
% %     endIdxLocal   = min(figIdx * motifsPerFig, nMotifSel);
% %     motifsThisFig = endIdxLocal - startIdxLocal + 1;
% % 
% %     % Original motif IDs shown in this figure
% %     motifIDsThisFig = selectMotifs(startIdxLocal:endIdxLocal);
% % 
% %     h = figure('Color', 'w');
% % 
% %     % Layout margins (figure-normalized units) reserved for the title,
% %     % per-column frame headers, and per-row motif labels.
% %     leftMarginFrac   = 0.06;
% %     topMarginFrac    = 0.075;
% %     bottomMarginFrac = 0.01;
% %     rightMarginFrac  = 0.01;
% % 
% %     nCols = nFramesShow;
% %     nRows = motifsThisFig;
% % 
% %     tileW = (1 - leftMarginFrac - rightMarginFrac - (nCols - 1) * colGapFrac) / nCols;
% %     tileH = (1 - topMarginFrac  - bottomMarginFrac - (nRows - 1) * rowGapFrac) / nRows;
% % 
% %     % Each tile calls axis(ax,'image','off') to preserve the (square)
% %     % image's aspect ratio. If a tile's allocated box isn't itself square,
% %     % MATLAB centers the square image within a letterboxed sub-region of
% %     % that box -- padding that colGapFrac has no control over, since it
% %     % isn't axes spacing at all. Fix this at the source: size the figure's
% %     % physical pixels so tileW/tileH normalized fractions correspond to an
% %     % actually-square physical tile, and axis('image') never needs to pad.
% %     targetTilePx = 220;
% %     figWidthPx  = targetTilePx / tileW;
% %     figHeightPx = targetTilePx / tileH;
% %     set(h, 'Units', 'pixels', 'Position', [100, 100, figWidthPx, figHeightPx]);
% % 
% %     for i = 1:motifsThisFig
% % 
% %         motifIdxLocal = startIdxLocal + i - 1;
% % 
% %         %% ----------------------------------------------------------------
% %         %  Reconstruct motif into image stack
% %         % -----------------------------------------------------------------
% % 
% %         if P == 64 * 64
% % 
% %             motif = reshape(squeeze(motif_w_kept(:, motifIdxLocal, :)), ...
% %                 64, 64, []);
% % 
% %         else
% % 
% %             motif = conditionDffMat( ...
% %                 squeeze(motif_w_kept(:, motifIdxLocal, :))', ...
% %                 nanpxs);
% %         end
% % 
% %         %% ----------------------------------------------------------------
% %         %  Optional Gaussian smoothing (native resolution, before scaling)
% %         % -----------------------------------------------------------------
% % 
% %         if smoothLogic
% %             motif = applyImgaussfilt(motif, 'sigma', gaussianSigma);
% %         end
% % 
% %         %% ----------------------------------------------------------------
% %         %  Scaling
% %         % -----------------------------------------------------------------
% % 
% %         switch lower(scaleMode)
% % 
% %             case 'permotif'
% % 
% %                 motif = localScaleByPercentile(motif, prctileRange);
% % 
% %             case 'global'
% % 
% %                 motif = localScaleByFixedRange(motif, globalLow, globalHigh);
% % 
% %             case 'none'
% % 
% %                 % Leave motif as-is.
% %         end
% % 
% %         %% ----------------------------------------------------------------
% %         %  Render each frame as its own tile (mask-aware upsampling +
% %         %  transparent background + gray cortex boundary overlay)
% %         % -----------------------------------------------------------------
% % 
% %         for fIdx = 1:nFramesShow
% % 
% %             frameData = motif(:, :, fIdx);
% %             frameData(~isfinite(frameData)) = 0;
% % 
% %             hiresData = localUpsampleFrameData(frameData, upsampleFactor, interpMethod);
% %             hiresMask = hiresMaskShared;
% % 
% %             % Bake the soft (anti-aliased) mask directly into an opaque RGB
% %             % image by blending the colormapped data with a white
% %             % background, rather than using continuous AlphaData. The
% %             % 'painters' renderer (used below for -dpdf printing) does not
% %             % reliably support partial/continuous transparency -- passing
% %             % it a non-binary AlphaData matrix can crash MATLAB. Producing
% %             % a fully opaque true-color image sidesteps that entirely while
% %             % still giving the same smooth, anti-aliased edge.
% %             rgbTile = localComposeRGBWithWhiteBackground( ...
% %                 hiresData, hiresMask, displayRange, cmap);
% % 
% %             tileX = leftMarginFrac + (fIdx - 1) * (tileW + colGapFrac);
% %             tileY = 1 - topMarginFrac - i * tileH - (i - 1) * rowGapFrac;
% % 
% %             ax = axes('Parent', h, 'Position', [tileX, tileY, tileW, tileH]); %#ok<LAXES>
% %             image(ax, rgbTile);
% %             axis(ax, 'image', 'off');
% %             set(ax, 'Color', 'w');
% %             hold(ax, 'on');
% % 
% %             if cortexBoundaryLogic
% %                 for bIdx = 1:numel(cortexBoundaryXY)
% %                     bxy = cortexBoundaryXY{bIdx} * upsampleFactor;
% %                     % plot() draws point 1 -> N but does not automatically
% %                     % connect N back to 1 -- explicitly close the loop by
% %                     % repeating the first point at the end, otherwise a gap
% %                     % appears at the seam (wherever the raw boundary trace
% %                     % happened to start).
% %                     bxyClosed = [bxy; bxy(1, :)];
% %                     plot(ax, bxyClosed(:,1), bxyClosed(:,2), '-', ...
% %                         'Color', boundaryColor, ...
% %                         'LineWidth', boundaryLineWidth);
% %                 end
% %             end
% % 
% %             hold(ax, 'off');
% % 
% %             % Column header (original frame index) on the first row only,
% %             % and only when meaningful (i.e. a single shared frame axis).
% %             if i == 1 && strcmp(frameSelectionScope, 'global')
% %                 title(ax, sprintf('f%d', keptFrameLabels.global(fIdx)), ...
% %                     'Color', 'k', 'FontSize', 8, 'FontWeight', 'normal');
% %             end
% % 
% %             % Row label (motif ID [+ frame range if perMotif]) on first column
% %             if fIdx == 1
% %                 if strcmp(frameSelectionScope, 'permotif')
% %                     rowLabel = sprintf('M%d (f%d-%d)', ...
% %                         motifIDsThisFig(i), ...
% %                         keptFrameLabels.perMotif(motifIdxLocal, 1), ...
% %                         keptFrameLabels.perMotif(motifIdxLocal, end));
% %                 else
% %                     rowLabel = sprintf('M%d', motifIDsThisFig(i));
% %                 end
% %                 ylabel(ax, rowLabel, 'Color', 'k', 'FontSize', 8, ...
% %                     'Rotation', 0, 'HorizontalAlignment', 'right', ...
% %                     'VerticalAlignment', 'middle', 'Visible', 'on');
% %                 % axis('off') hides the ylabel too, so force it visible
% %                 ax.YLabel.Visible = 'on';
% %             end
% %         end
% %     end
% % 
% %     if smoothLogic
% %         smoothLabel = sprintf(' | Gaussian \\sigma = %.2g', gaussianSigma);
% %     else
% %         smoothLabel = '';
% %     end
% % 
% %     if strcmp(frameSelectionScope, 'global') && nFramesShow < nFramesTotal
% %         frameLabel = sprintf(' | frames %d-%d of %d', ...
% %             keptFrameLabels.global(1), keptFrameLabels.global(end), nFramesTotal);
% %     elseif strcmp(frameSelectionScope, 'permotif') && nFramesShow < nFramesTotal
% %         frameLabel = sprintf(' | %d/%d frames (per-motif window)', ...
% %             nFramesShow, nFramesTotal);
% %     else
% %         frameLabel = '';
% %     end
% % 
% %     sgtitle(h, sprintf('Motifs %s%s%s', ...
% %         compressMotifIDs(motifIDsThisFig), frameLabel, smoothLabel), ...
% %         'Color', 'k', 'Interpreter', 'tex');
% % 
% %     %% --------------------------------------------------------------------
% %     %  Print montage
% %     % ---------------------------------------------------------------------
% % 
% %     if printLogic
% % 
% %         timestampStr  = datestr(now, 'mmddyy_HHMMSS');
% %         motifLabelStr = compressMotifIDs(motifIDsThisFig);
% % 
% %         if smoothLogic
% %             smoothFileStr = sprintf('_gaussSigma%.2g', gaussianSigma);
% %             smoothFileStr = strrep(smoothFileStr, '.', 'p');
% %         else
% %             smoothFileStr = '';
% %         end
% % 
% %         figSaveName = sprintf('%s_%s%s_%s', ...
% %             figNamePrefix, ...
% %             motifLabelStr, ...
% %             smoothFileStr, ...
% %             timestampStr);
% % 
% %         set(h, 'InvertHardcopy', 'off');  % preserve exact on-screen colors
% % 
% %         % Each tile is now a full true-color RGB bitmap (see
% %         % localComposeRGBWithWhiteBackground), not a small indexed image --
% %         % 'painters' is a pure vector renderer and is not well-suited to
% %         % embedding large raster content in a PDF; it has been a source of
% %         % crashes/instability for exactly this kind of mixed raster+vector
% %         % (image tiles + boundary lines) figure. 'opengl' handles this
% %         % mixed content far more robustly.
% %         print(h, fullfile(saveFigDir, figSaveName), ...
% %             '-opengl', '-bestfit', '-dpdf');
% %     end
% % end
% % 
% % end
% % 
% % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % function motifScaled = localScaleByPercentile(motif, prctileRange)
% % 
% % vals = motif(:);
% % vals = vals(isfinite(vals));
% % vals = vals(vals ~= 0);
% % 
% % if isempty(vals)
% %     motifScaled = zeros(size(motif), 'like', motif);
% %     return
% % end
% % 
% % lo = prctile(vals, prctileRange(1));
% % hi = prctile(vals, prctileRange(2));
% % 
% % if hi <= lo
% %     lo = min(vals);
% %     hi = max(vals);
% % end
% % 
% % if hi <= lo
% %     motifScaled = zeros(size(motif), 'like', motif);
% %     return
% % end
% % 
% % motifScaled = (motif - lo) ./ (hi - lo);
% % motifScaled(motifScaled < 0) = 0;
% % motifScaled(motifScaled > 1) = 1;
% % 
% % end
% % 
% % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % function motifScaled = localScaleByFixedRange(motif, lo, hi)
% % 
% % if hi <= lo
% %     motifScaled = zeros(size(motif), 'like', motif);
% %     return
% % end
% % 
% % motifScaled = (motif - lo) ./ (hi - lo);
% % motifScaled(motifScaled < 0) = 0;
% % motifScaled(motifScaled > 1) = 1;
% % 
% % end
% % 
% % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % function outStr = compressMotifIDs(ids)
% % % compressMotifIDs  Convert motif ID vector into compact range string.
% % %
% % % Example:
% % %   [1 2 3 4 5 6 7 9 12 13 14 15] -> '1-7_9_12-15'
% % 
% % ids = unique(ids(:))';
% % 
% % if isempty(ids)
% %     outStr = '';
% %     return;
% % end
% % 
% % rangeParts = {};
% % rangeStart = ids(1);
% % prevVal = ids(1);
% % 
% % for ii = 2:numel(ids)
% % 
% %     if ids(ii) == prevVal + 1
% % 
% %         prevVal = ids(ii);
% % 
% %     else
% % 
% %         rangeParts{end+1} = localRangeToStr(rangeStart, prevVal); %#ok<AGROW>
% %         rangeStart = ids(ii);
% %         prevVal = ids(ii);
% %     end
% % end
% % 
% % rangeParts{end+1} = localRangeToStr(rangeStart, prevVal);
% % outStr = strjoin(rangeParts, '_');
% % 
% % end
% % 
% % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % function s = localRangeToStr(a, b)
% % 
% % if a == b
% %     s = sprintf('%d', a);
% % else
% %     s = sprintf('%d-%d', a, b);
% % end
% % 
% % end
% % 
% % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % function cmap = localResolveColormap(colormapName)
% % % localResolveColormap
% % %
% % % Resolves a colormap from either:
% % %   - a string/char name, e.g. 'magma', 'parula', 'hot', 'turbo'
% % %   - an explicit N x 3 colormap matrix
% % 
% % if isnumeric(colormapName)
% %     cmap = colormapName;
% %     return
% % end
% % 
% % cmapName = char(colormapName);
% % 
% % try
% %     cmap = feval(cmapName, 256);
% % catch
% %     try
% %         cmap = eval(cmapName); %#ok<EVLDIR>
% %     catch
% %         warning('Could not resolve colormap "%s". Falling back to parula.', cmapName);
% %         cmap = parula(256);
% %     end
% % end
% % 
% % end
% % 
% % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % function rgbImage = localComposeRGBWithWhiteBackground(data, alphaMask, displayRange, cmap)
% % % localComposeRGBWithWhiteBackground
% % %
% % % Manually maps `data` through `cmap` (using displayRange as the color
% % % axis limits, matching what imagesc/clim would normally do) and blends
% % % the result with a white background using alphaMask (a continuous [0,1]
% % % field), producing a fully opaque H x W x 3 RGB image.
% % %
% % % This exists specifically to AVOID passing a continuous (non-binary)
% % % AlphaData matrix to imagesc: the 'painters' renderer (used elsewhere in
% % % this function for -dpdf printing) does not reliably support partial
% % % transparency, and doing so can crash MATLAB. Baking the blend into plain
% % % RGB values sidesteps alpha compositing entirely while still producing
% % % the same smooth, anti-aliased edge.
% % 
% % lo = displayRange(1);
% % hi = displayRange(2);
% % 
% % normVal = (data - lo) ./ (hi - lo);
% % normVal = min(max(normVal, 0), 1);
% % 
% % nColors = size(cmap, 1);
% % colorIdx = round(normVal * (nColors - 1)) + 1;
% % colorIdx = min(max(colorIdx, 1), nColors);
% % 
% % R = reshape(cmap(colorIdx(:), 1), size(data));
% % G = reshape(cmap(colorIdx(:), 2), size(data));
% % B = reshape(cmap(colorIdx(:), 3), size(data));
% % 
% % a = alphaMask;
% % 
% % Rout = R .* a + 1 .* (1 - a);
% % Gout = G .* a + 1 .* (1 - a);
% % Bout = B .* a + 1 .* (1 - a);
% % 
% % rgbImage = cat(3, Rout, Gout, Bout);
% % 
% % end
% % 
% % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % function [cortexMask, boundaryXY] = localGetCortexMaskAndBoundary( ...
% %     nanpxs, wantBoundary, numHarmonics, resampleN)
% % % localGetCortexMaskAndBoundary
% % %
% % % Derives a 64x64 logical dorsal-cortex mask directly from nanpxs -- the
% % % linear indices (or logical mask) of the non-cortex ("NaN") pixels within
% % % the full 64x64 = 4096 pixel grid -- and its boundary as a cell array of
% % % [x y] coordinate lists (native 64x64 pixel units, one cell per connected
% % % boundary component), smoothed via truncated Fourier descriptors.
% % %
% % % IMPORTANT: whether a real cortex mask exists is entirely a function of
% % % nanpxs, NOT of P (the first dimension of motif_w). motif_w can be stored
% % % either as [nValidPixels x K x L] (reconstructed via conditionDffMat) or
% % % as an already-reshaped [4096 x K x L] full-grid array -- either way, if
% % % nanpxs is supplied it still marks the true dorsal-cortex boundary within
% % % that 64x64 grid and should be used.
% % %
% % % If nanpxs is empty, there is no mask information available and the whole
% % % 64x64 frame is treated as valid (boundaryXY is returned empty).
% % %
% % % Boundary smoothness: rather than blurring a rasterized mask (which only
% % % ever anti-aliases individual pixel-step edges and struggles to remove the
% % % macroscopic staircase left by the native 64x64 resolution), the actual
% % % pixel boundary is traced once with bwboundaries, resampled to
% % % resampleN uniformly arc-length-spaced points (needed for a well-posed
% % % Fourier truncation), represented as a complex sequence z = x + 1i*y, and
% % % reconstructed from only its lowest numHarmonics frequency components via
% % % FFT/IFFT. Truncating high frequencies of a closed, periodic curve
% % % guarantees a result that is both smooth AND exactly seamlessly closed --
% % % there is no "sigma in the wrong units" pitfall here, since the smoothing
% % % is applied directly to the curve's shape, independent of any raster
% % % resolution or supersampling choice.
% % 
% % boundaryXY = {};
% % 
% % if isempty(nanpxs)
% %     cortexMask = true(64, 64);
% %     return
% % end
% % 
% % nanFlagVec = false(64 * 64, 1);
% % 
% % if islogical(nanpxs)
% %     nanFlagVec(:) = nanpxs(:);
% % else
% %     nanFlagVec(nanpxs(:)) = true;
% % end
% % 
% % cortexMask = reshape(~nanFlagVec, 64, 64);
% % 
% % if ~wantBoundary
% %     return
% % end
% % 
% % if exist('bwboundaries', 'file') ~= 2
% %     warning('bwboundaries (Image Processing Toolbox) not found; skipping cortex boundary overlay.');
% %     return
% % end
% % 
% % rawBoundaries = bwboundaries(cortexMask, 'noholes');
% % boundaryXY = cell(size(rawBoundaries));
% % 
% % for k = 1:numel(rawBoundaries)
% %     % bwboundaries returns [row col] = [y x]; convert to [x y]
% %     rawXY = [rawBoundaries{k}(:,2), rawBoundaries{k}(:,1)];
% % 
% %     resampledXY = localResampleClosedCurve(rawXY, resampleN);
% %     boundaryXY{k} = localSmoothClosedContourFourier(resampledXY, numHarmonics);
% % end
% % 
% % end
% % 
% % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % function resampledXY = localResampleClosedCurve(xy, nResample)
% % % localResampleClosedCurve
% % %
% % % Resamples a closed 2D polygon (M x 2, [x y]) to nResample points spaced
% % % uniformly by arc length around the loop. Uniform spacing is required for
% % % a clean/well-posed Fourier-descriptor truncation afterward.
% % 
% % % Ensure explicitly closed (first point repeated at the end) for arc-length
% % % accumulation, then drop the duplicate after resampling.
% % if norm(xy(1,:) - xy(end,:)) > 1e-9
% %     xy = [xy; xy(1,:)];
% % end
% % 
% % % Guard against zero-length (duplicate/repeated) points -- bwboundaries can
% % % occasionally emit a repeated point where the traced path touches a thin
% % % or degenerate pixel connection. A zero-length segment here would give
% % % interp1 a non-strictly-increasing cumulative-distance vector, which can
% % % produce a small localized artifact in the curve after Fourier smoothing.
% % segLenRaw = sqrt(sum(diff(xy).^2, 2));
% % keepPoint = [true; segLenRaw > 1e-9];
% % xy = xy(keepPoint, :);
% % 
% % segLen  = sqrt(sum(diff(xy).^2, 2));
% % cumDist = [0; cumsum(segLen)];
% % totalLen = cumDist(end);
% % 
% % if totalLen == 0
% %     resampledXY = repmat(xy(1,:), nResample, 1);
% %     return
% % end
% % 
% % targetDist = linspace(0, totalLen, nResample + 1);
% % targetDist(end) = []; % drop duplicate closing point
% % 
% % xResampled = interp1(cumDist, xy(:,1), targetDist, 'linear');
% % yResampled = interp1(cumDist, xy(:,2), targetDist, 'linear');
% % 
% % resampledXY = [xResampled(:), yResampled(:)];
% % 
% % end
% % 
% % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % function smoothXY = localSmoothClosedContourFourier(xy, numHarmonics)
% % % localSmoothClosedContourFourier
% % %
% % % Smooths a closed, uniformly arc-length-resampled 2D curve (N x 2, [x y])
% % % by representing it as a complex sequence z = x + 1i*y, truncating its
% % % discrete Fourier transform to the lowest numHarmonics positive and
% % % negative frequencies (plus the DC term), and inverting. This is the
% % % classic "Fourier descriptor" smoothing approach: it is guaranteed to
% % % return an exactly closed curve (the representation is inherently
% % % periodic) that is smooth by construction, with numHarmonics as the sole,
% % % resolution-independent smoothness knob (fewer harmonics = smoother/
% % % rounder; more harmonics = closer to the original traced shape).
% % 
% % N = size(xy, 1);
% % z = complex(xy(:,1), xy(:,2));
% % 
% % Z = fft(z);
% % 
% % numHarmonics = min(numHarmonics, floor((N - 1) / 2));
% % 
% % keepMask = false(N, 1);
% % keepMask(1) = true; % DC term
% % keepMask(2:(numHarmonics + 1)) = true;               % low positive frequencies
% % keepMask((N - numHarmonics + 1):N) = true;           % mirrored negative frequencies
% % 
% % Z(~keepMask) = 0;
% % 
% % zSmooth = ifft(Z);
% % 
% % smoothXY = [real(zSmooth), imag(zSmooth)];
% % 
% % end
% % 
% % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % function hiresData = localUpsampleFrameData(frameData, upsampleFactor, interpMethod)
% % % localUpsampleFrameData
% % %
% % % Spatially upsamples a single 64x64 frame's data for high-resolution
% % % display. The data is interpolated as-is, with no masking applied here --
% % % masking/fading to white is handled entirely by the shared alpha mask
% % % (see localRasterizeSmoothMask), computed once from the same smooth
% % % boundary curve used for the plotted line. Zeroing data outside the mask
% % % before interpolation was tried and rejected: combined with a separately
% % % smoothed alpha, it causes visible color fringing (a dark tinge from the
% % % hard zero shows through wherever alpha is only partially faded).
% % 
% % if upsampleFactor == 1
% %     hiresData = frameData;
% %     return
% % end
% % 
% % targetSize = size(frameData) * upsampleFactor;
% % hiresData = imresize(frameData, targetSize, interpMethod);
% % 
% % end
% % 
% % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % function hiresMask = localRasterizeSmoothMask(cortexMask, boundaryXY, upsampleFactor)
% % % localRasterizeSmoothMask
% % %
% % % Builds the high-resolution alpha mask by rasterizing the SAME smooth
% % % boundary curve(s) used for the plotted gray line (via poly2mask), rather
% % % than independently bicubic-upsampling the original native-resolution
% % % binary mask. Using two different smoothing methods for the line (Fourier
% % % descriptors) and the mask (raster blurring) gives no guarantee they agree
% % % pixel-for-pixel, which produced a visible fringe/staircase along the
% % % edge in earlier attempts. Deriving both from the identical curve
% % % eliminates that mismatch structurally. A light final Gaussian pass
% % % anti-aliases the (now much finer, supersampled-grid-scale) rasterization
% % % step, which is a much smaller and less objectionable artifact than the
% % % native-pixel-grid blockiness this replaces.
% % %
% % % Falls back to a plain bicubic-upsampled binary mask if no boundary curve
% % % is available (e.g. cortexBoundaryLogic was false, or nanpxs was empty).
% % 
% % targetSize = size(cortexMask) * upsampleFactor;
% % 
% % if isempty(boundaryXY) || exist('poly2mask', 'file') ~= 2
% %     hiresMask = imresize(double(cortexMask), targetSize, 'bicubic');
% %     hiresMask = min(max(hiresMask, 0), 1);
% %     return
% % end
% % 
% % maskAccum = false(targetSize);
% % 
% % for k = 1:numel(boundaryXY)
% %     xy = boundaryXY{k} * upsampleFactor;
% %     bw = poly2mask(xy(:,1), xy(:,2), targetSize(1), targetSize(2));
% %     maskAccum = maskAccum | bw;
% % end
% % 
% % hiresMask = double(maskAccum);
% % 
% % if exist('imgaussfilt', 'file') == 2
% %     hiresMask = imgaussfilt(hiresMask, 1.0); % mild anti-aliasing at supersampled-grid scale
% % end
% % 
% % hiresMask = min(max(hiresMask, 0), 1);
% % 
% % end
% % 
% % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % function [motif_w_kept, keptFrameLabels] = localSelectFrameWindow(motif_w_sel, nFramesKeep, scope)
% % % localSelectFrameWindow
% % %
% % % Selects, for each motif (or globally across all motifs), the contiguous
% % % window of nFramesKeep frames with the highest total energy, out of the
% % % nFramesTotal frames available. Energy for a frame is the sum of squared
% % % pixel values across all pixels/motifs being pooled, ignoring NaNs.
% % %
% % % Because the chosen window is always contiguous, its complement (the
% % % dropped frames) is always a prefix and/or suffix of the sequence -- i.e.
% % % frames are dropped only from the edges, never from the middle, and frame
% % % order is preserved.
% % %
% % % Outputs:
% % %   motif_w_kept    : [P x nMotifSel x nFramesKeep] tensor with the
% % %                     dropped frames removed.
% % %   keptFrameLabels : struct with fields:
% % %                       .global   - 1 x nFramesKeep vector of original frame
% % %                                   indices (valid/used when scope=='global')
% % %                       .perMotif - nMotifSel x nFramesKeep matrix of
% % %                                   original frame indices per motif (valid/
% % %                                   used when scope=='perMotif')
% % 
% % [P, nMotifSel, nFramesTotal] = size(motif_w_sel);
% % 
% % keptFrameLabels = struct('global', [], 'perMotif', []);
% % 
% % if nFramesKeep == nFramesTotal
% %     motif_w_kept = motif_w_sel;
% %     keptFrameLabels.global   = 1:nFramesTotal;
% %     keptFrameLabels.perMotif = repmat(1:nFramesTotal, nMotifSel, 1);
% %     return
% % end
% % 
% % % Energy per motif per frame: [nMotifSel x nFramesTotal]
% % sq = motif_w_sel .^ 2;
% % sq(~isfinite(sq)) = 0;
% % energyMotifFrame = squeeze(sum(sq, 1));       % nMotifSel x nFramesTotal
% % if nMotifSel == 1
% %     energyMotifFrame = reshape(energyMotifFrame, 1, nFramesTotal);
% % end
% % 
% % switch scope
% % 
% %     case 'global'
% % 
% %         aggregateEnergy = sum(energyMotifFrame, 1);         % 1 x nFramesTotal
% %         winStart = localBestWindowStart(aggregateEnergy, nFramesKeep);
% %         keepIdx  = winStart:(winStart + nFramesKeep - 1);
% % 
% %         motif_w_kept = motif_w_sel(:, :, keepIdx);
% %         keptFrameLabels.global   = keepIdx;
% %         keptFrameLabels.perMotif = repmat(keepIdx, nMotifSel, 1);
% % 
% %     case 'permotif'
% % 
% %         motif_w_kept = zeros(P, nMotifSel, nFramesKeep, 'like', motif_w_sel);
% %         perMotifIdx  = zeros(nMotifSel, nFramesKeep);
% % 
% %         for m = 1:nMotifSel
% %             winStart = localBestWindowStart(energyMotifFrame(m, :), nFramesKeep);
% %             keepIdx  = winStart:(winStart + nFramesKeep - 1);
% % 
% %             motif_w_kept(:, m, :) = motif_w_sel(:, m, keepIdx);
% %             perMotifIdx(m, :) = keepIdx;
% %         end
% % 
% %         keptFrameLabels.perMotif = perMotifIdx;
% % 
% %     otherwise
% % 
% %         error('Unknown frameSelectionScope: %s', scope);
% % end
% % 
% % end
% % 
% % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % function winStart = localBestWindowStart(energyVec, winLen)
% % % localBestWindowStart
% % %
% % % Slides a window of length winLen across energyVec (1 x L) and returns the
% % % start index of the window with the maximum total energy. Ties are broken
% % % in favor of the earliest (smallest-index) window.
% % 
% % L = numel(energyVec);
% % nWindows = L - winLen + 1;
% % 
% % winSums = zeros(1, nWindows);
% % for s = 1:nWindows
% %     winSums(s) = sum(energyVec(s:(s + winLen - 1)));
% % end
% % 
% % [~, winStart] = max(winSums);
% % 
% % end
% % 
% % 
% % % function montageMotifsPrintAdvanced(motif_w, nanpxs, varargin)
% % % % MontageMotifsPrintAdvanced  Display and optionally print motif montages.
% % % %
% % % % Usage:
% % % %   MontageMotifsPrintAdvanced(motif_w, nanpxs)
% % % %   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'motifsPerFig', 3)
% % % %   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'printLogic', false)
% % % %   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'selectMotifs', [1 3 6 7 9 12 13 14 15])
% % % %   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'prctileRange', [1 99.5])
% % % %   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'smoothLogic', true, 'gaussianSigma', 1)
% % % %   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'nFrames', 7)
% % % %   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'nFrames', 7, 'frameSelectionScope', 'perMotif')
% % % %   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'upsampleFactor', 6, 'interpMethod', 'bicubic')
% % % %
% % % % Inputs:
% % % %   motif_w      : [P x K x L] array of spatiotemporal motifs
% % % %                  pixels x motifs x frames/lags
% % % %   nanpxs       : NaN pixel information used by conditionDffMat. Also used
% % % %                  here to derive the dorsal-cortex mask/boundary.
% % % %
% % % % Name-value pairs:
% % % %   'motifsPerFig'  : number of motifs per figure, default = 1
% % % %   'printLogic'    : logical scalar, whether to print to PDF, default = true
% % % %   'selectMotifs'  : vector of original motif IDs to display, default = all
% % % %   'prctileRange'  : percentile range for per-motif scaling, default = [1 99.5]
% % % %   'scaleMode'     : 'perMotif', 'global', or 'none', default = 'perMotif'
% % % %   'displayRange'  : display range after scaling, default = [0 1]
% % % %   'colormapName'  : colormap name or N x 3 matrix, default = 'magma'
% % % %   'figNamePrefix' : output figure filename prefix, default = 'Motifs'
% % % %   'smoothLogic'   : apply Gaussian smoothing frame-by-frame, default = false
% % % %   'gaussianSigma' : sigma for Gaussian smoothing, default = 1
% % % %
% % % %   'nFrames'             : number of frames to display per motif, out of the
% % % %                           L available. Default = [] (show all L frames).
% % % %                           When nFrames < L, the dropped frames are chosen
% % % %                           by an energy-based contiguous-window search (see
% % % %                           "Frame-dropping logic" below) -- frames are never
% % % %                           dropped out of the middle of the sequence.
% % % %   'frameSelectionScope' : 'global' (default) or 'perMotif'.
% % % %                           'global'   - a single common window of frames is
% % % %                                        chosen using the pooled (summed)
% % % %                                        energy across all displayed motifs,
% % % %                                        so every row/tile in the montage
% % % %                                        shares the same underlying frame
% % % %                                        indices (recommended: keeps a
% % % %                                        common time axis across motifs).
% % % %                           'perMotif' - each motif independently keeps its
% % % %                                        own best window. Frame indices may
% % % %                                        then differ row-to-row; the kept
% % % %                                        frame range is annotated on each row.
% % % %
% % % %   'cortexBoundaryLogic' : draw the dorsal-cortex boundary in gray,
% % % %                           default = true (only applies when nanpxs actually
% % % %                           crops the image, i.e. P ~= 64*64).
% % % %   'boundaryColor'        : RGB triplet for the boundary line, default = [0.5 0.5 0.5]
% % % %   'boundaryLineWidth'    : line width for the boundary, default = 1.5
% % % %   'boundaryNumHarmonics' : number of low-frequency Fourier harmonics kept
% % % %                           when smoothing the cortex boundary curve,
% % % %                           default = 15. The raw pixel boundary is traced
% % % %                           once, resampled to uniform arc-length spacing,
% % % %                           and reconstructed from only these harmonics --
% % % %                           fewer harmonics = smoother/rounder outline
% % % %                           (small anatomical notches may be smoothed away);
% % % %                           more harmonics = closer to the raw pixel shape.
% % % %   'boundaryResamplePoints' : number of uniformly arc-length-spaced points
% % % %                           the raw boundary is resampled to before Fourier
% % % %                           smoothing, default = 400. Should comfortably
% % % %                           exceed 2x boundaryNumHarmonics.
% % % %
% % % %   'colGapFrac' : gap between adjacent frame columns, as a fraction of the
% % % %                  figure width, default = 0.004 (very tight). Set
% % % %                  independently from 'rowGapFrac' (tiles are laid out with
% % % %                  manually positioned axes rather than tiledlayout, since
% % % %                  tiledlayout's 'TileSpacing' cannot differ by direction).
% % % %   'rowGapFrac' : gap between adjacent motif rows, as a fraction of the
% % % %                  figure height, default = 0.018.
% % % %
% % % %   'upsampleFactor' : spatial upsampling factor applied to each 64x64 frame
% % % %                      before display/printing, default = 4 (i.e. 256x256).
% % % %                      Interpolation is mask-aware (see Notes) to avoid
% % % %                      bleeding intensity across the cortex boundary.
% % % %   'interpMethod'   : interpolation method passed to imresize for the
% % % %                      upsampling step, default = 'bicubic'.
% % % %
% % % % Notes:
% % % %   - Gaussian smoothing and percentile/global scaling are computed at the
% % % %     native 64x64 resolution (as before); spatial upsampling happens last,
% % % %     purely for display/print rendering.
% % % %   - Background is now white. Pixels outside the dorsal-cortex mask are
% % % %     blended toward white (with a soft, anti-aliased edge) and baked
% % % %     directly into an opaque RGB image, rather than relying on continuous
% % % %     AlphaData transparency -- the 'painters' renderer used for -dpdf
% % % %     printing does not reliably support partial alpha and can crash on it.
% % % %   - Per-motif scaling rescales each motif across all lags using the
% % % %     requested percentile range.
% % % %   - Printed filenames preserve original motif IDs.
% % % %
% % % % Frame-dropping logic:
% % % %   Given L available frames and a request to keep nFrames <= L, the
% % % %   function computes, for each frame, an "energy" value (sum of squared
% % % %   pixel values across the cortex, ignoring NaNs). It then slides a window
% % % %   of length nFrames across the 1..L sequence and keeps whichever
% % % %   contiguous window has the greatest total energy. Because the window is
% % % %   contiguous, its complement (the dropped frames) is always a prefix
% % % %   and/or suffix of the sequence -- frames are never dropped from the
% % % %   middle, and frame order is preserved.
% % % 
% % % %% ------------------------------------------------------------------------
% % % %  Save directory
% % % % -------------------------------------------------------------------------
% % % 
% % % saveFigDir = '/Volumes/buschman/Rodent Data/dualImaging_parkj/collectFigure/montage';
% % % 
% % % if ispc
% % %     saveFigDir = 'Z:\Rodent Data\dualImaging_parkj\collectFigure\montage';
% % % end
% % % 
% % % if exist(saveFigDir, 'dir') ~= 7
% % %     mkdir(saveFigDir);
% % % end
% % % 
% % % %% ------------------------------------------------------------------------
% % % %  Parse inputs
% % % % -------------------------------------------------------------------------
% % % 
% % % p = inputParser;
% % % 
% % % p.addParameter('motifsPerFig', 1, ...
% % %     @(x) isnumeric(x) && isscalar(x) && x >= 1);
% % % 
% % % p.addParameter('printLogic', true, ...
% % %     @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
% % % 
% % % p.addParameter('selectMotifs', [], ...
% % %     @(x) isempty(x) || isnumeric(x));
% % % 
% % % p.addParameter('prctileRange', [1 99.5], ...
% % %     @(x) isnumeric(x) && numel(x) == 2 && x(1) < x(2));
% % % 
% % % p.addParameter('scaleMode', 'perMotif', ...
% % %     @(x) ischar(x) || isstring(x));
% % % 
% % % p.addParameter('displayRange', [0 1], ...
% % %     @(x) isnumeric(x) && numel(x) == 2 && x(1) < x(2));
% % % 
% % % p.addParameter('colormapName', 'magma', ...
% % %     @(x) ischar(x) || isstring(x) || isnumeric(x));
% % % 
% % % p.addParameter('figNamePrefix', 'Motifs', ...
% % %     @(x) ischar(x) || isstring(x));
% % % 
% % % p.addParameter('smoothLogic', false, ...
% % %     @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
% % % 
% % % p.addParameter('gaussianSigma', 1, ...
% % %     @(x) isnumeric(x) && isscalar(x) && x > 0);
% % % 
% % % % New: frame-count control + drop logic
% % % p.addParameter('nFrames', [], ...
% % %     @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 1 && mod(x,1) == 0));
% % % 
% % % p.addParameter('frameSelectionScope', 'global', ...
% % %     @(x) ischar(x) || isstring(x));
% % % 
% % % % New: cortex boundary overlay
% % % p.addParameter('cortexBoundaryLogic', true, ...
% % %     @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
% % % 
% % % p.addParameter('boundaryColor', [0.5 0.5 0.5], ...
% % %     @(x) isnumeric(x) && numel(x) == 3);
% % % 
% % % p.addParameter('boundaryLineWidth', 1.5, ...
% % %     @(x) isnumeric(x) && isscalar(x) && x > 0);
% % % 
% % % % Boundary smoothing via Fourier descriptors (see localGetCortexMaskAndBoundary):
% % % % the raw pixel boundary is traced once, resampled to uniform arc-length
% % % % spacing, and reconstructed from only its low-frequency components,
% % % % guaranteeing a smooth, seamlessly closed curve.
% % % % New: sub-pixel boundary smoothness controls (see localGetCortexMaskAndBoundary)
% % % p.addParameter('boundaryNumHarmonics', 15, ...
% % %     @(x) isnumeric(x) && isscalar(x) && x >= 1);
% % % 
% % % p.addParameter('boundaryResamplePoints', 400, ...
% % %     @(x) isnumeric(x) && isscalar(x) && x >= 16);
% % % 
% % % % New: independent horizontal/vertical tile spacing (fraction of figure
% % % % width/height given to the gap between adjacent tiles). Unlike
% % % % tiledlayout's 'TileSpacing', these are set independently per axis.
% % % p.addParameter('colGapFrac', 0.004, ...
% % %     @(x) isnumeric(x) && isscalar(x) && x >= 0 && x < 1);
% % % 
% % % p.addParameter('rowGapFrac', 0.018, ...
% % %     @(x) isnumeric(x) && isscalar(x) && x >= 0 && x < 1);
% % % 
% % % % New: spatial upsampling
% % % p.addParameter('upsampleFactor', 4, ...
% % %     @(x) isnumeric(x) && isscalar(x) && x >= 1);
% % % 
% % % p.addParameter('interpMethod', 'bicubic', ...
% % %     @(x) ischar(x) || isstring(x));
% % % 
% % % p.parse(varargin{:});
% % % 
% % % motifsPerFig    = p.Results.motifsPerFig;
% % % printLogic      = logical(p.Results.printLogic);
% % % selectMotifs    = p.Results.selectMotifs;
% % % prctileRange    = p.Results.prctileRange;
% % % scaleMode       = char(p.Results.scaleMode);
% % % displayRange    = p.Results.displayRange;
% % % colormapName    = p.Results.colormapName;
% % % figNamePrefix   = char(p.Results.figNamePrefix);
% % % smoothLogic     = logical(p.Results.smoothLogic);
% % % gaussianSigma   = p.Results.gaussianSigma;
% % % 
% % % nFramesRequest      = p.Results.nFrames;
% % % frameSelectionScope = lower(char(p.Results.frameSelectionScope));
% % % 
% % % cortexBoundaryLogic = logical(p.Results.cortexBoundaryLogic);
% % % boundaryColor       = p.Results.boundaryColor;
% % % boundaryLineWidth   = p.Results.boundaryLineWidth;
% % % boundaryNumHarmonics   = p.Results.boundaryNumHarmonics;
% % % boundaryResamplePoints = p.Results.boundaryResamplePoints;
% % % colGapFrac          = p.Results.colGapFrac;
% % % rowGapFrac          = p.Results.rowGapFrac;
% % % 
% % % upsampleFactor = p.Results.upsampleFactor;
% % % interpMethod   = char(p.Results.interpMethod);
% % % 
% % % if ~ismember(frameSelectionScope, {'global', 'permotif'})
% % %     error('frameSelectionScope must be ''global'' or ''perMotif''.');
% % % end
% % % 
% % % %% ------------------------------------------------------------------------
% % % %  Validate motif tensor
% % % % -------------------------------------------------------------------------
% % % 
% % % if ~isnumeric(motif_w) || ndims(motif_w) ~= 3
% % %     error('motif_w must be a numeric [P x K x L] array.');
% % % end
% % % 
% % % [P, nMotifTotal, nFramesTotal] = size(motif_w);
% % % 
% % % if isempty(selectMotifs)
% % %     selectMotifs = 1:nMotifTotal;
% % % else
% % %     selectMotifs = unique(selectMotifs(:))';
% % % 
% % %     if any(selectMotifs < 1) || any(selectMotifs > nMotifTotal)
% % %         error('selectMotifs contains indices outside valid motif range 1:%d.', nMotifTotal);
% % %     end
% % % end
% % % 
% % % % Restrict displayed motifs, but preserve original IDs in selectMotifs
% % % motif_w_sel = motif_w(:, selectMotifs, :);
% % % 
% % % nMotifSel = numel(selectMotifs);
% % % 
% % % %% ------------------------------------------------------------------------
% % % %  Cortex mask + boundary (used for transparency and the gray outline)
% % % % -------------------------------------------------------------------------
% % % 
% % % [cortexMask, cortexBoundaryXY] = localGetCortexMaskAndBoundary( ...
% % %     nanpxs, cortexBoundaryLogic, boundaryNumHarmonics, boundaryResamplePoints);
% % % 
% % % if isempty(cortexBoundaryXY)
% % %     cortexBoundaryLogic = false;
% % % end
% % % 
% % % % Precompute the high-resolution alpha mask ONCE, by rasterizing the exact
% % % % same smooth boundary curve used for the plotted line (via poly2mask)
% % % % rather than separately bicubic-upsampling the original blocky binary
% % % % mask. Deriving the mask and the line from two different smoothing
% % % % methods (Fourier descriptors for one, raster blurring for the other) has
% % % % no guarantee of agreeing pixel-for-pixel, which is what kept showing up
% % % % as a residual fringe/staircase along the edge. Using the identical curve
% % % % for both eliminates that mismatch entirely. This is also computed once
% % % % here rather than per motif/frame, since it doesn't depend on the data.
% % % hiresMaskShared = localRasterizeSmoothMask( ...
% % %     cortexMask, cortexBoundaryXY, upsampleFactor);
% % % 
% % % %% ------------------------------------------------------------------------
% % % %  Frame selection (drop only from the edges, keep highest-energy window)
% % % % -------------------------------------------------------------------------
% % % 
% % % if isempty(nFramesRequest)
% % %     nFramesKeep = nFramesTotal;
% % % else
% % %     nFramesKeep = nFramesRequest;
% % % end
% % % 
% % % if nFramesKeep > nFramesTotal
% % %     error('nFrames (%d) cannot exceed the number of available frames (%d).', ...
% % %         nFramesKeep, nFramesTotal);
% % % end
% % % 
% % % [motif_w_kept, keptFrameLabels] = localSelectFrameWindow( ...
% % %     motif_w_sel, nFramesKeep, frameSelectionScope);
% % % 
% % % nMotifSel_check = size(motif_w_kept, 2); %#ok<NASGU>
% % % nFramesShow     = size(motif_w_kept, 3);
% % % 
% % % nFigures = ceil(nMotifSel / motifsPerFig);
% % % 
% % % %% ------------------------------------------------------------------------
% % % %  Optional global scaling range (computed on the frames actually shown)
% % % % -------------------------------------------------------------------------
% % % 
% % % switch lower(scaleMode)
% % % 
% % %     case 'global'
% % % 
% % %         vals = motif_w_kept(:);
% % %         vals = vals(isfinite(vals));
% % %         vals = vals(vals ~= 0);
% % % 
% % %         if isempty(vals)
% % %             globalLow = 0;
% % %             globalHigh = 1;
% % %         else
% % %             globalLow  = prctile(vals, prctileRange(1));
% % %             globalHigh = prctile(vals, prctileRange(2));
% % % 
% % %             if globalHigh <= globalLow
% % %                 globalLow  = min(vals);
% % %                 globalHigh = max(vals);
% % %             end
% % % 
% % %             if globalHigh <= globalLow
% % %                 globalLow  = 0;
% % %                 globalHigh = 1;
% % %             end
% % %         end
% % % 
% % %     case {'permotif', 'none'}
% % % 
% % %         globalLow = [];
% % %         globalHigh = [];
% % % 
% % %     otherwise
% % % 
% % %         error('scaleMode must be ''perMotif'', ''global'', or ''none''.');
% % % end
% % % 
% % % %% ------------------------------------------------------------------------
% % % %  Resolve colormap once (used to manually bake RGB below -- see note in
% % % %  the render loop about why we don't rely on continuous AlphaData)
% % % % -------------------------------------------------------------------------
% % % 
% % % cmap = localResolveColormap(colormapName);
% % % 
% % % %% ------------------------------------------------------------------------
% % % %  Main montage loop
% % % % -------------------------------------------------------------------------
% % % 
% % % for figIdx = 1:nFigures
% % % 
% % %     % Indices within selected motif list
% % %     startIdxLocal = (figIdx - 1) * motifsPerFig + 1;
% % %     endIdxLocal   = min(figIdx * motifsPerFig, nMotifSel);
% % %     motifsThisFig = endIdxLocal - startIdxLocal + 1;
% % % 
% % %     % Original motif IDs shown in this figure
% % %     motifIDsThisFig = selectMotifs(startIdxLocal:endIdxLocal);
% % % 
% % %     h = figure('Color', 'w');
% % % 
% % %     % Layout margins (figure-normalized units) reserved for the title,
% % %     % per-column frame headers, and per-row motif labels.
% % %     leftMarginFrac   = 0.06;
% % %     topMarginFrac    = 0.075;
% % %     bottomMarginFrac = 0.01;
% % %     rightMarginFrac  = 0.01;
% % % 
% % %     nCols = nFramesShow;
% % %     nRows = motifsThisFig;
% % % 
% % %     tileW = (1 - leftMarginFrac - rightMarginFrac - (nCols - 1) * colGapFrac) / nCols;
% % %     tileH = (1 - topMarginFrac  - bottomMarginFrac - (nRows - 1) * rowGapFrac) / nRows;
% % % 
% % %     for i = 1:motifsThisFig
% % % 
% % %         motifIdxLocal = startIdxLocal + i - 1;
% % % 
% % %         %% ----------------------------------------------------------------
% % %         %  Reconstruct motif into image stack
% % %         % -----------------------------------------------------------------
% % % 
% % %         if P == 64 * 64
% % % 
% % %             motif = reshape(squeeze(motif_w_kept(:, motifIdxLocal, :)), ...
% % %                 64, 64, []);
% % % 
% % %         else
% % % 
% % %             motif = conditionDffMat( ...
% % %                 squeeze(motif_w_kept(:, motifIdxLocal, :))', ...
% % %                 nanpxs);
% % %         end
% % % 
% % %         %% ----------------------------------------------------------------
% % %         %  Optional Gaussian smoothing (native resolution, before scaling)
% % %         % -----------------------------------------------------------------
% % % 
% % %         if smoothLogic
% % %             motif = applyImgaussfilt(motif, 'sigma', gaussianSigma);
% % %         end
% % % 
% % %         %% ----------------------------------------------------------------
% % %         %  Scaling
% % %         % -----------------------------------------------------------------
% % % 
% % %         switch lower(scaleMode)
% % % 
% % %             case 'permotif'
% % % 
% % %                 motif = localScaleByPercentile(motif, prctileRange);
% % % 
% % %             case 'global'
% % % 
% % %                 motif = localScaleByFixedRange(motif, globalLow, globalHigh);
% % % 
% % %             case 'none'
% % % 
% % %                 % Leave motif as-is.
% % %         end
% % % 
% % %         %% ----------------------------------------------------------------
% % %         %  Render each frame as its own tile (mask-aware upsampling +
% % %         %  transparent background + gray cortex boundary overlay)
% % %         % -----------------------------------------------------------------
% % % 
% % %         for fIdx = 1:nFramesShow
% % % 
% % %             frameData = motif(:, :, fIdx);
% % %             frameData(~isfinite(frameData)) = 0;
% % % 
% % %             hiresData = localUpsampleFrameData(frameData, upsampleFactor, interpMethod);
% % %             hiresMask = hiresMaskShared;
% % % 
% % %             % Bake the soft (anti-aliased) mask directly into an opaque RGB
% % %             % image by blending the colormapped data with a white
% % %             % background, rather than using continuous AlphaData. The
% % %             % 'painters' renderer (used below for -dpdf printing) does not
% % %             % reliably support partial/continuous transparency -- passing
% % %             % it a non-binary AlphaData matrix can crash MATLAB. Producing
% % %             % a fully opaque true-color image sidesteps that entirely while
% % %             % still giving the same smooth, anti-aliased edge.
% % %             rgbTile = localComposeRGBWithWhiteBackground( ...
% % %                 hiresData, hiresMask, displayRange, cmap);
% % % 
% % %             tileX = leftMarginFrac + (fIdx - 1) * (tileW + colGapFrac);
% % %             tileY = 1 - topMarginFrac - i * tileH - (i - 1) * rowGapFrac;
% % % 
% % %             ax = axes('Parent', h, 'Position', [tileX, tileY, tileW, tileH]); %#ok<LAXES>
% % %             image(ax, rgbTile);
% % %             axis(ax, 'image', 'off');
% % %             set(ax, 'Color', 'w');
% % %             hold(ax, 'on');
% % % 
% % %             if cortexBoundaryLogic
% % %                 for bIdx = 1:numel(cortexBoundaryXY)
% % %                     bxy = cortexBoundaryXY{bIdx} * upsampleFactor;
% % %                     % plot() draws point 1 -> N but does not automatically
% % %                     % connect N back to 1 -- explicitly close the loop by
% % %                     % repeating the first point at the end, otherwise a gap
% % %                     % appears at the seam (wherever the raw boundary trace
% % %                     % happened to start).
% % %                     bxyClosed = [bxy; bxy(1, :)];
% % %                     plot(ax, bxyClosed(:,1), bxyClosed(:,2), '-', ...
% % %                         'Color', boundaryColor, ...
% % %                         'LineWidth', boundaryLineWidth);
% % %                 end
% % %             end
% % % 
% % %             hold(ax, 'off');
% % % 
% % %             % Column header (original frame index) on the first row only,
% % %             % and only when meaningful (i.e. a single shared frame axis).
% % %             if i == 1 && strcmp(frameSelectionScope, 'global')
% % %                 title(ax, sprintf('f%d', keptFrameLabels.global(fIdx)), ...
% % %                     'Color', 'k', 'FontSize', 8, 'FontWeight', 'normal');
% % %             end
% % % 
% % %             % Row label (motif ID [+ frame range if perMotif]) on first column
% % %             if fIdx == 1
% % %                 if strcmp(frameSelectionScope, 'permotif')
% % %                     rowLabel = sprintf('M%d (f%d-%d)', ...
% % %                         motifIDsThisFig(i), ...
% % %                         keptFrameLabels.perMotif(motifIdxLocal, 1), ...
% % %                         keptFrameLabels.perMotif(motifIdxLocal, end));
% % %                 else
% % %                     rowLabel = sprintf('M%d', motifIDsThisFig(i));
% % %                 end
% % %                 ylabel(ax, rowLabel, 'Color', 'k', 'FontSize', 8, ...
% % %                     'Rotation', 0, 'HorizontalAlignment', 'right', ...
% % %                     'VerticalAlignment', 'middle', 'Visible', 'on');
% % %                 % axis('off') hides the ylabel too, so force it visible
% % %                 ax.YLabel.Visible = 'on';
% % %             end
% % %         end
% % %     end
% % % 
% % %     if smoothLogic
% % %         smoothLabel = sprintf(' | Gaussian \\sigma = %.2g', gaussianSigma);
% % %     else
% % %         smoothLabel = '';
% % %     end
% % % 
% % %     if strcmp(frameSelectionScope, 'global') && nFramesShow < nFramesTotal
% % %         frameLabel = sprintf(' | frames %d-%d of %d', ...
% % %             keptFrameLabels.global(1), keptFrameLabels.global(end), nFramesTotal);
% % %     elseif strcmp(frameSelectionScope, 'permotif') && nFramesShow < nFramesTotal
% % %         frameLabel = sprintf(' | %d/%d frames (per-motif window)', ...
% % %             nFramesShow, nFramesTotal);
% % %     else
% % %         frameLabel = '';
% % %     end
% % % 
% % %     sgtitle(h, sprintf('Motifs %s%s%s', ...
% % %         compressMotifIDs(motifIDsThisFig), frameLabel, smoothLabel), ...
% % %         'Color', 'k', 'Interpreter', 'tex');
% % % 
% % %     %% --------------------------------------------------------------------
% % %     %  Print montage
% % %     % ---------------------------------------------------------------------
% % % 
% % %     if printLogic
% % % 
% % %         timestampStr  = datestr(now, 'mmddyy_HHMMSS');
% % %         motifLabelStr = compressMotifIDs(motifIDsThisFig);
% % % 
% % %         if smoothLogic
% % %             smoothFileStr = sprintf('_gaussSigma%.2g', gaussianSigma);
% % %             smoothFileStr = strrep(smoothFileStr, '.', 'p');
% % %         else
% % %             smoothFileStr = '';
% % %         end
% % % 
% % %         figSaveName = sprintf('%s_%s%s_%s', ...
% % %             figNamePrefix, ...
% % %             motifLabelStr, ...
% % %             smoothFileStr, ...
% % %             timestampStr);
% % % 
% % %         set(h, 'InvertHardcopy', 'off');  % preserve exact on-screen colors
% % % 
% % %         % Each tile is now a full true-color RGB bitmap (see
% % %         % localComposeRGBWithWhiteBackground), not a small indexed image --
% % %         % 'painters' is a pure vector renderer and is not well-suited to
% % %         % embedding large raster content in a PDF; it has been a source of
% % %         % crashes/instability for exactly this kind of mixed raster+vector
% % %         % (image tiles + boundary lines) figure. 'opengl' handles this
% % %         % mixed content far more robustly.
% % %         print(h, fullfile(saveFigDir, figSaveName), ...
% % %             '-painters', '-bestfit', '-dpdf');
% % %     end
% % % end
% % % 
% % % end
% % % 
% % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % function motifScaled = localScaleByPercentile(motif, prctileRange)
% % % 
% % % vals = motif(:);
% % % vals = vals(isfinite(vals));
% % % vals = vals(vals ~= 0);
% % % 
% % % if isempty(vals)
% % %     motifScaled = zeros(size(motif), 'like', motif);
% % %     return
% % % end
% % % 
% % % lo = prctile(vals, prctileRange(1));
% % % hi = prctile(vals, prctileRange(2));
% % % 
% % % if hi <= lo
% % %     lo = min(vals);
% % %     hi = max(vals);
% % % end
% % % 
% % % if hi <= lo
% % %     motifScaled = zeros(size(motif), 'like', motif);
% % %     return
% % % end
% % % 
% % % motifScaled = (motif - lo) ./ (hi - lo);
% % % motifScaled(motifScaled < 0) = 0;
% % % motifScaled(motifScaled > 1) = 1;
% % % 
% % % end
% % % 
% % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % function motifScaled = localScaleByFixedRange(motif, lo, hi)
% % % 
% % % if hi <= lo
% % %     motifScaled = zeros(size(motif), 'like', motif);
% % %     return
% % % end
% % % 
% % % motifScaled = (motif - lo) ./ (hi - lo);
% % % motifScaled(motifScaled < 0) = 0;
% % % motifScaled(motifScaled > 1) = 1;
% % % 
% % % end
% % % 
% % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % function outStr = compressMotifIDs(ids)
% % % % compressMotifIDs  Convert motif ID vector into compact range string.
% % % %
% % % % Example:
% % % %   [1 2 3 4 5 6 7 9 12 13 14 15] -> '1-7_9_12-15'
% % % 
% % % ids = unique(ids(:))';
% % % 
% % % if isempty(ids)
% % %     outStr = '';
% % %     return;
% % % end
% % % 
% % % rangeParts = {};
% % % rangeStart = ids(1);
% % % prevVal = ids(1);
% % % 
% % % for ii = 2:numel(ids)
% % % 
% % %     if ids(ii) == prevVal + 1
% % % 
% % %         prevVal = ids(ii);
% % % 
% % %     else
% % % 
% % %         rangeParts{end+1} = localRangeToStr(rangeStart, prevVal); %#ok<AGROW>
% % %         rangeStart = ids(ii);
% % %         prevVal = ids(ii);
% % %     end
% % % end
% % % 
% % % rangeParts{end+1} = localRangeToStr(rangeStart, prevVal);
% % % outStr = strjoin(rangeParts, '_');
% % % 
% % % end
% % % 
% % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % function s = localRangeToStr(a, b)
% % % 
% % % if a == b
% % %     s = sprintf('%d', a);
% % % else
% % %     s = sprintf('%d-%d', a, b);
% % % end
% % % 
% % % end
% % % 
% % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % function cmap = localResolveColormap(colormapName)
% % % % localResolveColormap
% % % %
% % % % Resolves a colormap from either:
% % % %   - a string/char name, e.g. 'magma', 'parula', 'hot', 'turbo'
% % % %   - an explicit N x 3 colormap matrix
% % % 
% % % if isnumeric(colormapName)
% % %     cmap = colormapName;
% % %     return
% % % end
% % % 
% % % cmapName = char(colormapName);
% % % 
% % % try
% % %     cmap = feval(cmapName, 256);
% % % catch
% % %     try
% % %         cmap = eval(cmapName); %#ok<EVLDIR>
% % %     catch
% % %         warning('Could not resolve colormap "%s". Falling back to parula.', cmapName);
% % %         cmap = parula(256);
% % %     end
% % % end
% % % 
% % % end
% % % 
% % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % function rgbImage = localComposeRGBWithWhiteBackground(data, alphaMask, displayRange, cmap)
% % % % localComposeRGBWithWhiteBackground
% % % %
% % % % Manually maps `data` through `cmap` (using displayRange as the color
% % % % axis limits, matching what imagesc/clim would normally do) and blends
% % % % the result with a white background using alphaMask (a continuous [0,1]
% % % % field), producing a fully opaque H x W x 3 RGB image.
% % % %
% % % % This exists specifically to AVOID passing a continuous (non-binary)
% % % % AlphaData matrix to imagesc: the 'painters' renderer (used elsewhere in
% % % % this function for -dpdf printing) does not reliably support partial
% % % % transparency, and doing so can crash MATLAB. Baking the blend into plain
% % % % RGB values sidesteps alpha compositing entirely while still producing
% % % % the same smooth, anti-aliased edge.
% % % 
% % % lo = displayRange(1);
% % % hi = displayRange(2);
% % % 
% % % normVal = (data - lo) ./ (hi - lo);
% % % normVal = min(max(normVal, 0), 1);
% % % 
% % % nColors = size(cmap, 1);
% % % colorIdx = round(normVal * (nColors - 1)) + 1;
% % % colorIdx = min(max(colorIdx, 1), nColors);
% % % 
% % % R = reshape(cmap(colorIdx(:), 1), size(data));
% % % G = reshape(cmap(colorIdx(:), 2), size(data));
% % % B = reshape(cmap(colorIdx(:), 3), size(data));
% % % 
% % % a = alphaMask;
% % % 
% % % Rout = R .* a + 1 .* (1 - a);
% % % Gout = G .* a + 1 .* (1 - a);
% % % Bout = B .* a + 1 .* (1 - a);
% % % 
% % % rgbImage = cat(3, Rout, Gout, Bout);
% % % 
% % % end
% % % 
% % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % function [cortexMask, boundaryXY] = localGetCortexMaskAndBoundary( ...
% % %     nanpxs, wantBoundary, numHarmonics, resampleN)
% % % % localGetCortexMaskAndBoundary
% % % %
% % % % Derives a 64x64 logical dorsal-cortex mask directly from nanpxs -- the
% % % % linear indices (or logical mask) of the non-cortex ("NaN") pixels within
% % % % the full 64x64 = 4096 pixel grid -- and its boundary as a cell array of
% % % % [x y] coordinate lists (native 64x64 pixel units, one cell per connected
% % % % boundary component), smoothed via truncated Fourier descriptors.
% % % %
% % % % IMPORTANT: whether a real cortex mask exists is entirely a function of
% % % % nanpxs, NOT of P (the first dimension of motif_w). motif_w can be stored
% % % % either as [nValidPixels x K x L] (reconstructed via conditionDffMat) or
% % % % as an already-reshaped [4096 x K x L] full-grid array -- either way, if
% % % % nanpxs is supplied it still marks the true dorsal-cortex boundary within
% % % % that 64x64 grid and should be used.
% % % %
% % % % If nanpxs is empty, there is no mask information available and the whole
% % % % 64x64 frame is treated as valid (boundaryXY is returned empty).
% % % %
% % % % Boundary smoothness: rather than blurring a rasterized mask (which only
% % % % ever anti-aliases individual pixel-step edges and struggles to remove the
% % % % macroscopic staircase left by the native 64x64 resolution), the actual
% % % % pixel boundary is traced once with bwboundaries, resampled to
% % % % resampleN uniformly arc-length-spaced points (needed for a well-posed
% % % % Fourier truncation), represented as a complex sequence z = x + 1i*y, and
% % % % reconstructed from only its lowest numHarmonics frequency components via
% % % % FFT/IFFT. Truncating high frequencies of a closed, periodic curve
% % % % guarantees a result that is both smooth AND exactly seamlessly closed --
% % % % there is no "sigma in the wrong units" pitfall here, since the smoothing
% % % % is applied directly to the curve's shape, independent of any raster
% % % % resolution or supersampling choice.
% % % 
% % % boundaryXY = {};
% % % 
% % % if isempty(nanpxs)
% % %     cortexMask = true(64, 64);
% % %     return
% % % end
% % % 
% % % nanFlagVec = false(64 * 64, 1);
% % % 
% % % if islogical(nanpxs)
% % %     nanFlagVec(:) = nanpxs(:);
% % % else
% % %     nanFlagVec(nanpxs(:)) = true;
% % % end
% % % 
% % % cortexMask = reshape(~nanFlagVec, 64, 64);
% % % 
% % % if ~wantBoundary
% % %     return
% % % end
% % % 
% % % if exist('bwboundaries', 'file') ~= 2
% % %     warning('bwboundaries (Image Processing Toolbox) not found; skipping cortex boundary overlay.');
% % %     return
% % % end
% % % 
% % % rawBoundaries = bwboundaries(cortexMask, 'noholes');
% % % boundaryXY = cell(size(rawBoundaries));
% % % 
% % % for k = 1:numel(rawBoundaries)
% % %     % bwboundaries returns [row col] = [y x]; convert to [x y]
% % %     rawXY = [rawBoundaries{k}(:,2), rawBoundaries{k}(:,1)];
% % % 
% % %     resampledXY = localResampleClosedCurve(rawXY, resampleN);
% % %     boundaryXY{k} = localSmoothClosedContourFourier(resampledXY, numHarmonics);
% % % end
% % % 
% % % end
% % % 
% % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % function resampledXY = localResampleClosedCurve(xy, nResample)
% % % % localResampleClosedCurve
% % % %
% % % % Resamples a closed 2D polygon (M x 2, [x y]) to nResample points spaced
% % % % uniformly by arc length around the loop. Uniform spacing is required for
% % % % a clean/well-posed Fourier-descriptor truncation afterward.
% % % 
% % % % Ensure explicitly closed (first point repeated at the end) for arc-length
% % % % accumulation, then drop the duplicate after resampling.
% % % if norm(xy(1,:) - xy(end,:)) > 1e-9
% % %     xy = [xy; xy(1,:)];
% % % end
% % % 
% % % % Guard against zero-length (duplicate/repeated) points -- bwboundaries can
% % % % occasionally emit a repeated point where the traced path touches a thin
% % % % or degenerate pixel connection. A zero-length segment here would give
% % % % interp1 a non-strictly-increasing cumulative-distance vector, which can
% % % % produce a small localized artifact in the curve after Fourier smoothing.
% % % segLenRaw = sqrt(sum(diff(xy).^2, 2));
% % % keepPoint = [true; segLenRaw > 1e-9];
% % % xy = xy(keepPoint, :);
% % % 
% % % segLen  = sqrt(sum(diff(xy).^2, 2));
% % % cumDist = [0; cumsum(segLen)];
% % % totalLen = cumDist(end);
% % % 
% % % if totalLen == 0
% % %     resampledXY = repmat(xy(1,:), nResample, 1);
% % %     return
% % % end
% % % 
% % % targetDist = linspace(0, totalLen, nResample + 1);
% % % targetDist(end) = []; % drop duplicate closing point
% % % 
% % % xResampled = interp1(cumDist, xy(:,1), targetDist, 'linear');
% % % yResampled = interp1(cumDist, xy(:,2), targetDist, 'linear');
% % % 
% % % resampledXY = [xResampled(:), yResampled(:)];
% % % 
% % % end
% % % 
% % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % function smoothXY = localSmoothClosedContourFourier(xy, numHarmonics)
% % % % localSmoothClosedContourFourier
% % % %
% % % % Smooths a closed, uniformly arc-length-resampled 2D curve (N x 2, [x y])
% % % % by representing it as a complex sequence z = x + 1i*y, truncating its
% % % % discrete Fourier transform to the lowest numHarmonics positive and
% % % % negative frequencies (plus the DC term), and inverting. This is the
% % % % classic "Fourier descriptor" smoothing approach: it is guaranteed to
% % % % return an exactly closed curve (the representation is inherently
% % % % periodic) that is smooth by construction, with numHarmonics as the sole,
% % % % resolution-independent smoothness knob (fewer harmonics = smoother/
% % % % rounder; more harmonics = closer to the original traced shape).
% % % 
% % % N = size(xy, 1);
% % % z = complex(xy(:,1), xy(:,2));
% % % 
% % % Z = fft(z);
% % % 
% % % numHarmonics = min(numHarmonics, floor((N - 1) / 2));
% % % 
% % % keepMask = false(N, 1);
% % % keepMask(1) = true; % DC term
% % % keepMask(2:(numHarmonics + 1)) = true;               % low positive frequencies
% % % keepMask((N - numHarmonics + 1):N) = true;           % mirrored negative frequencies
% % % 
% % % Z(~keepMask) = 0;
% % % 
% % % zSmooth = ifft(Z);
% % % 
% % % smoothXY = [real(zSmooth), imag(zSmooth)];
% % % 
% % % end
% % % 
% % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % function hiresData = localUpsampleFrameData(frameData, upsampleFactor, interpMethod)
% % % % localUpsampleFrameData
% % % %
% % % % Spatially upsamples a single 64x64 frame's data for high-resolution
% % % % display. The data is interpolated as-is, with no masking applied here --
% % % % masking/fading to white is handled entirely by the shared alpha mask
% % % % (see localRasterizeSmoothMask), computed once from the same smooth
% % % % boundary curve used for the plotted line. Zeroing data outside the mask
% % % % before interpolation was tried and rejected: combined with a separately
% % % % smoothed alpha, it causes visible color fringing (a dark tinge from the
% % % % hard zero shows through wherever alpha is only partially faded).
% % % 
% % % if upsampleFactor == 1
% % %     hiresData = frameData;
% % %     return
% % % end
% % % 
% % % targetSize = size(frameData) * upsampleFactor;
% % % hiresData = imresize(frameData, targetSize, interpMethod);
% % % 
% % % end
% % % 
% % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % function hiresMask = localRasterizeSmoothMask(cortexMask, boundaryXY, upsampleFactor)
% % % % localRasterizeSmoothMask
% % % %
% % % % Builds the high-resolution alpha mask by rasterizing the SAME smooth
% % % % boundary curve(s) used for the plotted gray line (via poly2mask), rather
% % % % than independently bicubic-upsampling the original native-resolution
% % % % binary mask. Using two different smoothing methods for the line (Fourier
% % % % descriptors) and the mask (raster blurring) gives no guarantee they agree
% % % % pixel-for-pixel, which produced a visible fringe/staircase along the
% % % % edge in earlier attempts. Deriving both from the identical curve
% % % % eliminates that mismatch structurally. A light final Gaussian pass
% % % % anti-aliases the (now much finer, supersampled-grid-scale) rasterization
% % % % step, which is a much smaller and less objectionable artifact than the
% % % % native-pixel-grid blockiness this replaces.
% % % %
% % % % Falls back to a plain bicubic-upsampled binary mask if no boundary curve
% % % % is available (e.g. cortexBoundaryLogic was false, or nanpxs was empty).
% % % 
% % % targetSize = size(cortexMask) * upsampleFactor;
% % % 
% % % if isempty(boundaryXY) || exist('poly2mask', 'file') ~= 2
% % %     hiresMask = imresize(double(cortexMask), targetSize, 'bicubic');
% % %     hiresMask = min(max(hiresMask, 0), 1);
% % %     return
% % % end
% % % 
% % % maskAccum = false(targetSize);
% % % 
% % % for k = 1:numel(boundaryXY)
% % %     xy = boundaryXY{k} * upsampleFactor;
% % %     bw = poly2mask(xy(:,1), xy(:,2), targetSize(1), targetSize(2));
% % %     maskAccum = maskAccum | bw;
% % % end
% % % 
% % % hiresMask = double(maskAccum);
% % % 
% % % if exist('imgaussfilt', 'file') == 2
% % %     hiresMask = imgaussfilt(hiresMask, 1.0); % mild anti-aliasing at supersampled-grid scale
% % % end
% % % 
% % % hiresMask = min(max(hiresMask, 0), 1);
% % % 
% % % end
% % % 
% % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % function [motif_w_kept, keptFrameLabels] = localSelectFrameWindow(motif_w_sel, nFramesKeep, scope)
% % % % localSelectFrameWindow
% % % %
% % % % Selects, for each motif (or globally across all motifs), the contiguous
% % % % window of nFramesKeep frames with the highest total energy, out of the
% % % % nFramesTotal frames available. Energy for a frame is the sum of squared
% % % % pixel values across all pixels/motifs being pooled, ignoring NaNs.
% % % %
% % % % Because the chosen window is always contiguous, its complement (the
% % % % dropped frames) is always a prefix and/or suffix of the sequence -- i.e.
% % % % frames are dropped only from the edges, never from the middle, and frame
% % % % order is preserved.
% % % %
% % % % Outputs:
% % % %   motif_w_kept    : [P x nMotifSel x nFramesKeep] tensor with the
% % % %                     dropped frames removed.
% % % %   keptFrameLabels : struct with fields:
% % % %                       .global   - 1 x nFramesKeep vector of original frame
% % % %                                   indices (valid/used when scope=='global')
% % % %                       .perMotif - nMotifSel x nFramesKeep matrix of
% % % %                                   original frame indices per motif (valid/
% % % %                                   used when scope=='perMotif')
% % % 
% % % [P, nMotifSel, nFramesTotal] = size(motif_w_sel);
% % % 
% % % keptFrameLabels = struct('global', [], 'perMotif', []);
% % % 
% % % if nFramesKeep == nFramesTotal
% % %     motif_w_kept = motif_w_sel;
% % %     keptFrameLabels.global   = 1:nFramesTotal;
% % %     keptFrameLabels.perMotif = repmat(1:nFramesTotal, nMotifSel, 1);
% % %     return
% % % end
% % % 
% % % % Energy per motif per frame: [nMotifSel x nFramesTotal]
% % % sq = motif_w_sel .^ 2;
% % % sq(~isfinite(sq)) = 0;
% % % energyMotifFrame = squeeze(sum(sq, 1));       % nMotifSel x nFramesTotal
% % % if nMotifSel == 1
% % %     energyMotifFrame = reshape(energyMotifFrame, 1, nFramesTotal);
% % % end
% % % 
% % % switch scope
% % % 
% % %     case 'global'
% % % 
% % %         aggregateEnergy = sum(energyMotifFrame, 1);         % 1 x nFramesTotal
% % %         winStart = localBestWindowStart(aggregateEnergy, nFramesKeep);
% % %         keepIdx  = winStart:(winStart + nFramesKeep - 1);
% % % 
% % %         motif_w_kept = motif_w_sel(:, :, keepIdx);
% % %         keptFrameLabels.global   = keepIdx;
% % %         keptFrameLabels.perMotif = repmat(keepIdx, nMotifSel, 1);
% % % 
% % %     case 'permotif'
% % % 
% % %         motif_w_kept = zeros(P, nMotifSel, nFramesKeep, 'like', motif_w_sel);
% % %         perMotifIdx  = zeros(nMotifSel, nFramesKeep);
% % % 
% % %         for m = 1:nMotifSel
% % %             winStart = localBestWindowStart(energyMotifFrame(m, :), nFramesKeep);
% % %             keepIdx  = winStart:(winStart + nFramesKeep - 1);
% % % 
% % %             motif_w_kept(:, m, :) = motif_w_sel(:, m, keepIdx);
% % %             perMotifIdx(m, :) = keepIdx;
% % %         end
% % % 
% % %         keptFrameLabels.perMotif = perMotifIdx;
% % % 
% % %     otherwise
% % % 
% % %         error('Unknown frameSelectionScope: %s', scope);
% % % end
% % % 
% % % end
% % 
% % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % function winStart = localBestWindowStart(energyVec, winLen)
% % % % localBestWindowStart
% % % %
% % % % Slides a window of length winLen across energyVec (1 x L) and returns the
% % % % start index of the window with the maximum total energy. Ties are broken
% % % % in favor of the earliest (smallest-index) window.
% % % 
% % % L = numel(energyVec);
% % % nWindows = L - winLen + 1;
% % % 
% % % winSums = zeros(1, nWindows);
% % % for s = 1:nWindows
% % %     winSums(s) = sum(energyVec(s:(s + winLen - 1)));
% % % end
% % % 
% % % [~, winStart] = max(winSums);
% % % 
% % % end
% % 
% % 
% % % function montageMotifsPrintAdvanced(motif_w, nanpxs, varargin)
% % % % MontageMotifsPrintAdvanced  Display and optionally print motif montages.
% % % %
% % % % Usage:
% % % %   MontageMotifsPrintAdvanced(motif_w, nanpxs)
% % % %   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'motifsPerFig', 3)
% % % %   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'printLogic', false)
% % % %   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'selectMotifs', [1 3 6 7 9 12 13 14 15])
% % % %   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'prctileRange', [1 99.5])
% % % %   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'smoothLogic', true, 'gaussianSigma', 1)
% % % %   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'nFrames', 7)
% % % %   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'nFrames', 7, 'frameSelectionScope', 'perMotif')
% % % %   MontageMotifsPrintAdvanced(motif_w, nanpxs, 'upsampleFactor', 6, 'interpMethod', 'bicubic')
% % % %
% % % % Inputs:
% % % %   motif_w      : [P x K x L] array of spatiotemporal motifs
% % % %                  pixels x motifs x frames/lags
% % % %   nanpxs       : NaN pixel information used by conditionDffMat. Also used
% % % %                  here to derive the dorsal-cortex mask/boundary.
% % % %
% % % % Name-value pairs:
% % % %   'motifsPerFig'  : number of motifs per figure, default = 1
% % % %   'printLogic'    : logical scalar, whether to print to PDF, default = true
% % % %   'selectMotifs'  : vector of original motif IDs to display, default = all
% % % %   'prctileRange'  : percentile range for per-motif scaling, default = [1 99.5]
% % % %   'scaleMode'     : 'perMotif', 'global', or 'none', default = 'perMotif'
% % % %   'displayRange'  : display range after scaling, default = [0 1]
% % % %   'colormapName'  : colormap name or N x 3 matrix, default = 'magma'
% % % %   'figNamePrefix' : output figure filename prefix, default = 'Motifs'
% % % %   'smoothLogic'   : apply Gaussian smoothing frame-by-frame, default = false
% % % %   'gaussianSigma' : sigma for Gaussian smoothing, default = 1
% % % %
% % % %   'nFrames'             : number of frames to display per motif, out of the
% % % %                           L available. Default = [] (show all L frames).
% % % %                           When nFrames < L, the dropped frames are chosen
% % % %                           by an energy-based contiguous-window search (see
% % % %                           "Frame-dropping logic" below) -- frames are never
% % % %                           dropped out of the middle of the sequence.
% % % %   'frameSelectionScope' : 'global' (default) or 'perMotif'.
% % % %                           'global'   - a single common window of frames is
% % % %                                        chosen using the pooled (summed)
% % % %                                        energy across all displayed motifs,
% % % %                                        so every row/tile in the montage
% % % %                                        shares the same underlying frame
% % % %                                        indices (recommended: keeps a
% % % %                                        common time axis across motifs).
% % % %                           'perMotif' - each motif independently keeps its
% % % %                                        own best window. Frame indices may
% % % %                                        then differ row-to-row; the kept
% % % %                                        frame range is annotated on each row.
% % % %
% % % %   'cortexBoundaryLogic' : draw the dorsal-cortex boundary in gray,
% % % %                           default = true (only applies when nanpxs actually
% % % %                           crops the image, i.e. P ~= 64*64).
% % % %   'boundaryColor'        : RGB triplet for the boundary line, default = [0.5 0.5 0.5]
% % % %   'boundaryLineWidth'    : line width for the boundary, default = 1.5
% % % %   'boundaryNumHarmonics' : number of low-frequency Fourier harmonics kept
% % % %                           when smoothing the cortex boundary curve,
% % % %                           default = 15. The raw pixel boundary is traced
% % % %                           once, resampled to uniform arc-length spacing,
% % % %                           and reconstructed from only these harmonics --
% % % %                           fewer harmonics = smoother/rounder outline
% % % %                           (small anatomical notches may be smoothed away);
% % % %                           more harmonics = closer to the raw pixel shape.
% % % %   'boundaryResamplePoints' : number of uniformly arc-length-spaced points
% % % %                           the raw boundary is resampled to before Fourier
% % % %                           smoothing, default = 400. Should comfortably
% % % %                           exceed 2x boundaryNumHarmonics.
% % % %
% % % %   'colGapFrac' : gap between adjacent frame columns, as a fraction of the
% % % %                  figure width, default = 0.004 (very tight). Set
% % % %                  independently from 'rowGapFrac' (tiles are laid out with
% % % %                  manually positioned axes rather than tiledlayout, since
% % % %                  tiledlayout's 'TileSpacing' cannot differ by direction).
% % % %   'rowGapFrac' : gap between adjacent motif rows, as a fraction of the
% % % %                  figure height, default = 0.018.
% % % %
% % % %   'upsampleFactor' : spatial upsampling factor applied to each 64x64 frame
% % % %                      before display/printing, default = 4 (i.e. 256x256).
% % % %                      Interpolation is mask-aware (see Notes) to avoid
% % % %                      bleeding intensity across the cortex boundary.
% % % %   'interpMethod'   : interpolation method passed to imresize for the
% % % %                      upsampling step, default = 'bicubic'.
% % % %
% % % % Notes:
% % % %   - Gaussian smoothing and percentile/global scaling are computed at the
% % % %     native 64x64 resolution (as before); spatial upsampling happens last,
% % % %     purely for display/print rendering.
% % % %   - Background is now white. Pixels outside the dorsal-cortex mask are
% % % %     blended toward white (with a soft, anti-aliased edge) and baked
% % % %     directly into an opaque RGB image, rather than relying on continuous
% % % %     AlphaData transparency -- the 'painters' renderer used for -dpdf
% % % %     printing does not reliably support partial alpha and can crash on it.
% % % %   - Per-motif scaling rescales each motif across all lags using the
% % % %     requested percentile range.
% % % %   - Printed filenames preserve original motif IDs.
% % % %
% % % % Frame-dropping logic:
% % % %   Given L available frames and a request to keep nFrames <= L, the
% % % %   function computes, for each frame, an "energy" value (sum of squared
% % % %   pixel values across the cortex, ignoring NaNs). It then slides a window
% % % %   of length nFrames across the 1..L sequence and keeps whichever
% % % %   contiguous window has the greatest total energy. Because the window is
% % % %   contiguous, its complement (the dropped frames) is always a prefix
% % % %   and/or suffix of the sequence -- frames are never dropped from the
% % % %   middle, and frame order is preserved.
% % % 
% % % %% ------------------------------------------------------------------------
% % % %  Save directory
% % % % -------------------------------------------------------------------------
% % % 
% % % saveFigDir = '/Volumes/buschman/Rodent Data/dualImaging_parkj/collectFigure/montage';
% % % 
% % % if ispc
% % %     saveFigDir = 'Z:\Rodent Data\dualImaging_parkj\collectFigure\montage';
% % % end
% % % 
% % % if exist(saveFigDir, 'dir') ~= 7
% % %     mkdir(saveFigDir);
% % % end
% % % 
% % % %% ------------------------------------------------------------------------
% % % %  Parse inputs
% % % % -------------------------------------------------------------------------
% % % 
% % % p = inputParser;
% % % 
% % % p.addParameter('motifsPerFig', 1, ...
% % %     @(x) isnumeric(x) && isscalar(x) && x >= 1);
% % % 
% % % p.addParameter('printLogic', true, ...
% % %     @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
% % % 
% % % p.addParameter('selectMotifs', [], ...
% % %     @(x) isempty(x) || isnumeric(x));
% % % 
% % % p.addParameter('prctileRange', [1 99.5], ...
% % %     @(x) isnumeric(x) && numel(x) == 2 && x(1) < x(2));
% % % 
% % % p.addParameter('scaleMode', 'perMotif', ...
% % %     @(x) ischar(x) || isstring(x));
% % % 
% % % p.addParameter('displayRange', [0 1], ...
% % %     @(x) isnumeric(x) && numel(x) == 2 && x(1) < x(2));
% % % 
% % % p.addParameter('colormapName', 'magma', ...
% % %     @(x) ischar(x) || isstring(x) || isnumeric(x));
% % % 
% % % p.addParameter('figNamePrefix', 'Motifs', ...
% % %     @(x) ischar(x) || isstring(x));
% % % 
% % % p.addParameter('smoothLogic', false, ...
% % %     @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
% % % 
% % % p.addParameter('gaussianSigma', 1, ...
% % %     @(x) isnumeric(x) && isscalar(x) && x > 0);
% % % 
% % % % New: frame-count control + drop logic
% % % p.addParameter('nFrames', [], ...
% % %     @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 1 && mod(x,1) == 0));
% % % 
% % % p.addParameter('frameSelectionScope', 'global', ...
% % %     @(x) ischar(x) || isstring(x));
% % % 
% % % % New: cortex boundary overlay
% % % p.addParameter('cortexBoundaryLogic', true, ...
% % %     @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
% % % 
% % % p.addParameter('boundaryColor', [0.5 0.5 0.5], ...
% % %     @(x) isnumeric(x) && numel(x) == 3);
% % % 
% % % p.addParameter('boundaryLineWidth', 1.5, ...
% % %     @(x) isnumeric(x) && isscalar(x) && x > 0);
% % % 
% % % % Boundary smoothing via Fourier descriptors (see localGetCortexMaskAndBoundary):
% % % % the raw pixel boundary is traced once, resampled to uniform arc-length
% % % % spacing, and reconstructed from only its low-frequency components,
% % % % guaranteeing a smooth, seamlessly closed curve.
% % % % New: sub-pixel boundary smoothness controls (see localGetCortexMaskAndBoundary)
% % % p.addParameter('boundaryNumHarmonics', 15, ...
% % %     @(x) isnumeric(x) && isscalar(x) && x >= 1);
% % % 
% % % p.addParameter('boundaryResamplePoints', 400, ...
% % %     @(x) isnumeric(x) && isscalar(x) && x >= 16);
% % % 
% % % % New: independent horizontal/vertical tile spacing (fraction of figure
% % % % width/height given to the gap between adjacent tiles). Unlike
% % % % tiledlayout's 'TileSpacing', these are set independently per axis.
% % % p.addParameter('colGapFrac', 0.004, ...
% % %     @(x) isnumeric(x) && isscalar(x) && x >= 0 && x < 1);
% % % 
% % % p.addParameter('rowGapFrac', 0.018, ...
% % %     @(x) isnumeric(x) && isscalar(x) && x >= 0 && x < 1);
% % % 
% % % % New: spatial upsampling
% % % p.addParameter('upsampleFactor', 4, ...
% % %     @(x) isnumeric(x) && isscalar(x) && x >= 1);
% % % 
% % % p.addParameter('interpMethod', 'bicubic', ...
% % %     @(x) ischar(x) || isstring(x));
% % % 
% % % p.parse(varargin{:});
% % % 
% % % motifsPerFig    = p.Results.motifsPerFig;
% % % printLogic      = logical(p.Results.printLogic);
% % % selectMotifs    = p.Results.selectMotifs;
% % % prctileRange    = p.Results.prctileRange;
% % % scaleMode       = char(p.Results.scaleMode);
% % % displayRange    = p.Results.displayRange;
% % % colormapName    = p.Results.colormapName;
% % % figNamePrefix   = char(p.Results.figNamePrefix);
% % % smoothLogic     = logical(p.Results.smoothLogic);
% % % gaussianSigma   = p.Results.gaussianSigma;
% % % 
% % % nFramesRequest      = p.Results.nFrames;
% % % frameSelectionScope = lower(char(p.Results.frameSelectionScope));
% % % 
% % % cortexBoundaryLogic = logical(p.Results.cortexBoundaryLogic);
% % % boundaryColor       = p.Results.boundaryColor;
% % % boundaryLineWidth   = p.Results.boundaryLineWidth;
% % % boundaryNumHarmonics   = p.Results.boundaryNumHarmonics;
% % % boundaryResamplePoints = p.Results.boundaryResamplePoints;
% % % colGapFrac          = p.Results.colGapFrac;
% % % rowGapFrac          = p.Results.rowGapFrac;
% % % 
% % % upsampleFactor = p.Results.upsampleFactor;
% % % interpMethod   = char(p.Results.interpMethod);
% % % 
% % % if ~ismember(frameSelectionScope, {'global', 'permotif'})
% % %     error('frameSelectionScope must be ''global'' or ''perMotif''.');
% % % end
% % % 
% % % %% ------------------------------------------------------------------------
% % % %  Validate motif tensor
% % % % -------------------------------------------------------------------------
% % % 
% % % if ~isnumeric(motif_w) || ndims(motif_w) ~= 3
% % %     error('motif_w must be a numeric [P x K x L] array.');
% % % end
% % % 
% % % [P, nMotifTotal, nFramesTotal] = size(motif_w);
% % % 
% % % if isempty(selectMotifs)
% % %     selectMotifs = 1:nMotifTotal;
% % % else
% % %     selectMotifs = unique(selectMotifs(:))';
% % % 
% % %     if any(selectMotifs < 1) || any(selectMotifs > nMotifTotal)
% % %         error('selectMotifs contains indices outside valid motif range 1:%d.', nMotifTotal);
% % %     end
% % % end
% % % 
% % % % Restrict displayed motifs, but preserve original IDs in selectMotifs
% % % motif_w_sel = motif_w(:, selectMotifs, :);
% % % 
% % % nMotifSel = numel(selectMotifs);
% % % 
% % % %% ------------------------------------------------------------------------
% % % %  Cortex mask + boundary (used for transparency and the gray outline)
% % % % -------------------------------------------------------------------------
% % % 
% % % [cortexMask, cortexBoundaryXY] = localGetCortexMaskAndBoundary( ...
% % %     nanpxs, cortexBoundaryLogic, boundaryNumHarmonics, boundaryResamplePoints);
% % % 
% % % if isempty(cortexBoundaryXY)
% % %     cortexBoundaryLogic = false;
% % % end
% % % 
% % % % Precompute the high-resolution alpha mask ONCE, by rasterizing the exact
% % % % same smooth boundary curve used for the plotted line (via poly2mask)
% % % % rather than separately bicubic-upsampling the original blocky binary
% % % % mask. Deriving the mask and the line from two different smoothing
% % % % methods (Fourier descriptors for one, raster blurring for the other) has
% % % % no guarantee of agreeing pixel-for-pixel, which is what kept showing up
% % % % as a residual fringe/staircase along the edge. Using the identical curve
% % % % for both eliminates that mismatch entirely. This is also computed once
% % % % here rather than per motif/frame, since it doesn't depend on the data.
% % % hiresMaskShared = localRasterizeSmoothMask( ...
% % %     cortexMask, cortexBoundaryXY, upsampleFactor);
% % % 
% % % %% ------------------------------------------------------------------------
% % % %  Frame selection (drop only from the edges, keep highest-energy window)
% % % % -------------------------------------------------------------------------
% % % 
% % % if isempty(nFramesRequest)
% % %     nFramesKeep = nFramesTotal;
% % % else
% % %     nFramesKeep = nFramesRequest;
% % % end
% % % 
% % % if nFramesKeep > nFramesTotal
% % %     error('nFrames (%d) cannot exceed the number of available frames (%d).', ...
% % %         nFramesKeep, nFramesTotal);
% % % end
% % % 
% % % [motif_w_kept, keptFrameLabels] = localSelectFrameWindow( ...
% % %     motif_w_sel, nFramesKeep, frameSelectionScope);
% % % 
% % % nMotifSel_check = size(motif_w_kept, 2); %#ok<NASGU>
% % % nFramesShow     = size(motif_w_kept, 3);
% % % 
% % % nFigures = ceil(nMotifSel / motifsPerFig);
% % % 
% % % %% ------------------------------------------------------------------------
% % % %  Optional global scaling range (computed on the frames actually shown)
% % % % -------------------------------------------------------------------------
% % % 
% % % switch lower(scaleMode)
% % % 
% % %     case 'global'
% % % 
% % %         vals = motif_w_kept(:);
% % %         vals = vals(isfinite(vals));
% % %         vals = vals(vals ~= 0);
% % % 
% % %         if isempty(vals)
% % %             globalLow = 0;
% % %             globalHigh = 1;
% % %         else
% % %             globalLow  = prctile(vals, prctileRange(1));
% % %             globalHigh = prctile(vals, prctileRange(2));
% % % 
% % %             if globalHigh <= globalLow
% % %                 globalLow  = min(vals);
% % %                 globalHigh = max(vals);
% % %             end
% % % 
% % %             if globalHigh <= globalLow
% % %                 globalLow  = 0;
% % %                 globalHigh = 1;
% % %             end
% % %         end
% % % 
% % %     case {'permotif', 'none'}
% % % 
% % %         globalLow = [];
% % %         globalHigh = [];
% % % 
% % %     otherwise
% % % 
% % %         error('scaleMode must be ''perMotif'', ''global'', or ''none''.');
% % % end
% % % 
% % % %% ------------------------------------------------------------------------
% % % %  Resolve colormap once (used to manually bake RGB below -- see note in
% % % %  the render loop about why we don't rely on continuous AlphaData)
% % % % -------------------------------------------------------------------------
% % % 
% % % cmap = localResolveColormap(colormapName);
% % % 
% % % %% ------------------------------------------------------------------------
% % % %  Main montage loop
% % % % -------------------------------------------------------------------------
% % % 
% % % for figIdx = 1:nFigures
% % % 
% % %     % Indices within selected motif list
% % %     startIdxLocal = (figIdx - 1) * motifsPerFig + 1;
% % %     endIdxLocal   = min(figIdx * motifsPerFig, nMotifSel);
% % %     motifsThisFig = endIdxLocal - startIdxLocal + 1;
% % % 
% % %     % Original motif IDs shown in this figure
% % %     motifIDsThisFig = selectMotifs(startIdxLocal:endIdxLocal);
% % % 
% % %     h = figure('Color', 'w');
% % % 
% % %     % Layout margins (figure-normalized units) reserved for the title,
% % %     % per-column frame headers, and per-row motif labels.
% % %     leftMarginFrac   = 0.06;
% % %     topMarginFrac    = 0.075;
% % %     bottomMarginFrac = 0.01;
% % %     rightMarginFrac  = 0.01;
% % % 
% % %     nCols = nFramesShow;
% % %     nRows = motifsThisFig;
% % % 
% % %     tileW = (1 - leftMarginFrac - rightMarginFrac - (nCols - 1) * colGapFrac) / nCols;
% % %     tileH = (1 - topMarginFrac  - bottomMarginFrac - (nRows - 1) * rowGapFrac) / nRows;
% % % 
% % %     for i = 1:motifsThisFig
% % % 
% % %         motifIdxLocal = startIdxLocal + i - 1;
% % % 
% % %         %% ----------------------------------------------------------------
% % %         %  Reconstruct motif into image stack
% % %         % -----------------------------------------------------------------
% % % 
% % %         if P == 64 * 64
% % % 
% % %             motif = reshape(squeeze(motif_w_kept(:, motifIdxLocal, :)), ...
% % %                 64, 64, []);
% % % 
% % %         else
% % % 
% % %             motif = conditionDffMat( ...
% % %                 squeeze(motif_w_kept(:, motifIdxLocal, :))', ...
% % %                 nanpxs);
% % %         end
% % % 
% % %         %% ----------------------------------------------------------------
% % %         %  Optional Gaussian smoothing (native resolution, before scaling)
% % %         % -----------------------------------------------------------------
% % % 
% % %         if smoothLogic
% % %             motif = applyImgaussfilt(motif, 'sigma', gaussianSigma);
% % %         end
% % % 
% % %         %% ----------------------------------------------------------------
% % %         %  Scaling
% % %         % -----------------------------------------------------------------
% % % 
% % %         switch lower(scaleMode)
% % % 
% % %             case 'permotif'
% % % 
% % %                 motif = localScaleByPercentile(motif, prctileRange);
% % % 
% % %             case 'global'
% % % 
% % %                 motif = localScaleByFixedRange(motif, globalLow, globalHigh);
% % % 
% % %             case 'none'
% % % 
% % %                 % Leave motif as-is.
% % %         end
% % % 
% % %         %% ----------------------------------------------------------------
% % %         %  Render each frame as its own tile (mask-aware upsampling +
% % %         %  transparent background + gray cortex boundary overlay)
% % %         % -----------------------------------------------------------------
% % % 
% % %         for fIdx = 1:nFramesShow
% % % 
% % %             frameData = motif(:, :, fIdx);
% % %             frameData(~isfinite(frameData)) = 0;
% % % 
% % %             hiresData = localUpsampleFrameData(frameData, upsampleFactor, interpMethod);
% % %             hiresMask = hiresMaskShared;
% % % 
% % %             % Bake the soft (anti-aliased) mask directly into an opaque RGB
% % %             % image by blending the colormapped data with a white
% % %             % background, rather than using continuous AlphaData. The
% % %             % 'painters' renderer (used below for -dpdf printing) does not
% % %             % reliably support partial/continuous transparency -- passing
% % %             % it a non-binary AlphaData matrix can crash MATLAB. Producing
% % %             % a fully opaque true-color image sidesteps that entirely while
% % %             % still giving the same smooth, anti-aliased edge.
% % %             rgbTile = localComposeRGBWithWhiteBackground( ...
% % %                 hiresData, hiresMask, displayRange, cmap);
% % % 
% % %             tileX = leftMarginFrac + (fIdx - 1) * (tileW + colGapFrac);
% % %             tileY = 1 - topMarginFrac - i * tileH - (i - 1) * rowGapFrac;
% % % 
% % %             ax = axes('Parent', h, 'Position', [tileX, tileY, tileW, tileH]); %#ok<LAXES>
% % %             image(ax, rgbTile);
% % %             axis(ax, 'image', 'off');
% % %             set(ax, 'Color', 'w');
% % %             hold(ax, 'on');
% % % 
% % %             if cortexBoundaryLogic
% % %                 for bIdx = 1:numel(cortexBoundaryXY)
% % %                     bxy = cortexBoundaryXY{bIdx} * upsampleFactor;
% % %                     plot(ax, bxy(:,1), bxy(:,2), '-', ...
% % %                         'Color', boundaryColor, ...
% % %                         'LineWidth', boundaryLineWidth);
% % %                 end
% % %             end
% % % 
% % %             hold(ax, 'off');
% % % 
% % %             % Column header (original frame index) on the first row only,
% % %             % and only when meaningful (i.e. a single shared frame axis).
% % %             if i == 1 && strcmp(frameSelectionScope, 'global')
% % %                 title(ax, sprintf('f%d', keptFrameLabels.global(fIdx)), ...
% % %                     'Color', 'k', 'FontSize', 8, 'FontWeight', 'normal');
% % %             end
% % % 
% % %             % Row label (motif ID [+ frame range if perMotif]) on first column
% % %             if fIdx == 1
% % %                 if strcmp(frameSelectionScope, 'permotif')
% % %                     rowLabel = sprintf('M%d (f%d-%d)', ...
% % %                         motifIDsThisFig(i), ...
% % %                         keptFrameLabels.perMotif(motifIdxLocal, 1), ...
% % %                         keptFrameLabels.perMotif(motifIdxLocal, end));
% % %                 else
% % %                     rowLabel = sprintf('M%d', motifIDsThisFig(i));
% % %                 end
% % %                 ylabel(ax, rowLabel, 'Color', 'k', 'FontSize', 8, ...
% % %                     'Rotation', 0, 'HorizontalAlignment', 'right', ...
% % %                     'VerticalAlignment', 'middle', 'Visible', 'on');
% % %                 % axis('off') hides the ylabel too, so force it visible
% % %                 ax.YLabel.Visible = 'on';
% % %             end
% % %         end
% % %     end
% % % 
% % %     if smoothLogic
% % %         smoothLabel = sprintf(' | Gaussian \\sigma = %.2g', gaussianSigma);
% % %     else
% % %         smoothLabel = '';
% % %     end
% % % 
% % %     if strcmp(frameSelectionScope, 'global') && nFramesShow < nFramesTotal
% % %         frameLabel = sprintf(' | frames %d-%d of %d', ...
% % %             keptFrameLabels.global(1), keptFrameLabels.global(end), nFramesTotal);
% % %     elseif strcmp(frameSelectionScope, 'permotif') && nFramesShow < nFramesTotal
% % %         frameLabel = sprintf(' | %d/%d frames (per-motif window)', ...
% % %             nFramesShow, nFramesTotal);
% % %     else
% % %         frameLabel = '';
% % %     end
% % % 
% % %     sgtitle(h, sprintf('Motifs %s%s%s', ...
% % %         compressMotifIDs(motifIDsThisFig), frameLabel, smoothLabel), ...
% % %         'Color', 'k', 'Interpreter', 'tex');
% % % 
% % %     %% --------------------------------------------------------------------
% % %     %  Print montage
% % %     % ---------------------------------------------------------------------
% % % 
% % %     if printLogic
% % % 
% % %         timestampStr  = datestr(now, 'mmddyy_HHMMSS');
% % %         motifLabelStr = compressMotifIDs(motifIDsThisFig);
% % % 
% % %         if smoothLogic
% % %             smoothFileStr = sprintf('_gaussSigma%.2g', gaussianSigma);
% % %             smoothFileStr = strrep(smoothFileStr, '.', 'p');
% % %         else
% % %             smoothFileStr = '';
% % %         end
% % % 
% % %         figSaveName = sprintf('%s_%s%s_%s', ...
% % %             figNamePrefix, ...
% % %             motifLabelStr, ...
% % %             smoothFileStr, ...
% % %             timestampStr);
% % % 
% % %         set(h, 'InvertHardcopy', 'off');  % preserve exact on-screen colors
% % % 
% % %         print(h, fullfile(saveFigDir, figSaveName), ...
% % %             '-painters', '-bestfit', '-dpdf');
% % %     end
% % % end
% % % 
% % % end
% % % 
% % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % function motifScaled = localScaleByPercentile(motif, prctileRange)
% % % 
% % % vals = motif(:);
% % % vals = vals(isfinite(vals));
% % % vals = vals(vals ~= 0);
% % % 
% % % if isempty(vals)
% % %     motifScaled = zeros(size(motif), 'like', motif);
% % %     return
% % % end
% % % 
% % % lo = prctile(vals, prctileRange(1));
% % % hi = prctile(vals, prctileRange(2));
% % % 
% % % if hi <= lo
% % %     lo = min(vals);
% % %     hi = max(vals);
% % % end
% % % 
% % % if hi <= lo
% % %     motifScaled = zeros(size(motif), 'like', motif);
% % %     return
% % % end
% % % 
% % % motifScaled = (motif - lo) ./ (hi - lo);
% % % motifScaled(motifScaled < 0) = 0;
% % % motifScaled(motifScaled > 1) = 1;
% % % 
% % % end
% % % 
% % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % function motifScaled = localScaleByFixedRange(motif, lo, hi)
% % % 
% % % if hi <= lo
% % %     motifScaled = zeros(size(motif), 'like', motif);
% % %     return
% % % end
% % % 
% % % motifScaled = (motif - lo) ./ (hi - lo);
% % % motifScaled(motifScaled < 0) = 0;
% % % motifScaled(motifScaled > 1) = 1;
% % % 
% % % end
% % % 
% % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % function outStr = compressMotifIDs(ids)
% % % % compressMotifIDs  Convert motif ID vector into compact range string.
% % % %
% % % % Example:
% % % %   [1 2 3 4 5 6 7 9 12 13 14 15] -> '1-7_9_12-15'
% % % 
% % % ids = unique(ids(:))';
% % % 
% % % if isempty(ids)
% % %     outStr = '';
% % %     return;
% % % end
% % % 
% % % rangeParts = {};
% % % rangeStart = ids(1);
% % % prevVal = ids(1);
% % % 
% % % for ii = 2:numel(ids)
% % % 
% % %     if ids(ii) == prevVal + 1
% % % 
% % %         prevVal = ids(ii);
% % % 
% % %     else
% % % 
% % %         rangeParts{end+1} = localRangeToStr(rangeStart, prevVal); %#ok<AGROW>
% % %         rangeStart = ids(ii);
% % %         prevVal = ids(ii);
% % %     end
% % % end
% % % 
% % % rangeParts{end+1} = localRangeToStr(rangeStart, prevVal);
% % % outStr = strjoin(rangeParts, '_');
% % % 
% % % end
% % % 
% % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % function s = localRangeToStr(a, b)
% % % 
% % % if a == b
% % %     s = sprintf('%d', a);
% % % else
% % %     s = sprintf('%d-%d', a, b);
% % % end
% % % 
% % % end
% % % 
% % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % function cmap = localResolveColormap(colormapName)
% % % % localResolveColormap
% % % %
% % % % Resolves a colormap from either:
% % % %   - a string/char name, e.g. 'magma', 'parula', 'hot', 'turbo'
% % % %   - an explicit N x 3 colormap matrix
% % % 
% % % if isnumeric(colormapName)
% % %     cmap = colormapName;
% % %     return
% % % end
% % % 
% % % cmapName = char(colormapName);
% % % 
% % % try
% % %     cmap = feval(cmapName, 256);
% % % catch
% % %     try
% % %         cmap = eval(cmapName); %#ok<EVLDIR>
% % %     catch
% % %         warning('Could not resolve colormap "%s". Falling back to parula.', cmapName);
% % %         cmap = parula(256);
% % %     end
% % % end
% % % 
% % % end
% % % 
% % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % function rgbImage = localComposeRGBWithWhiteBackground(data, alphaMask, displayRange, cmap)
% % % % localComposeRGBWithWhiteBackground
% % % %
% % % % Manually maps `data` through `cmap` (using displayRange as the color
% % % % axis limits, matching what imagesc/clim would normally do) and blends
% % % % the result with a white background using alphaMask (a continuous [0,1]
% % % % field), producing a fully opaque H x W x 3 RGB image.
% % % %
% % % % This exists specifically to AVOID passing a continuous (non-binary)
% % % % AlphaData matrix to imagesc: the 'painters' renderer (used elsewhere in
% % % % this function for -dpdf printing) does not reliably support partial
% % % % transparency, and doing so can crash MATLAB. Baking the blend into plain
% % % % RGB values sidesteps alpha compositing entirely while still producing
% % % % the same smooth, anti-aliased edge.
% % % 
% % % lo = displayRange(1);
% % % hi = displayRange(2);
% % % 
% % % normVal = (data - lo) ./ (hi - lo);
% % % normVal = min(max(normVal, 0), 1);
% % % 
% % % nColors = size(cmap, 1);
% % % colorIdx = round(normVal * (nColors - 1)) + 1;
% % % colorIdx = min(max(colorIdx, 1), nColors);
% % % 
% % % R = reshape(cmap(colorIdx(:), 1), size(data));
% % % G = reshape(cmap(colorIdx(:), 2), size(data));
% % % B = reshape(cmap(colorIdx(:), 3), size(data));
% % % 
% % % a = alphaMask;
% % % 
% % % Rout = R .* a + 1 .* (1 - a);
% % % Gout = G .* a + 1 .* (1 - a);
% % % Bout = B .* a + 1 .* (1 - a);
% % % 
% % % rgbImage = cat(3, Rout, Gout, Bout);
% % % 
% % % end
% % % 
% % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % function [cortexMask, boundaryXY] = localGetCortexMaskAndBoundary( ...
% % %     nanpxs, wantBoundary, numHarmonics, resampleN)
% % % % localGetCortexMaskAndBoundary
% % % %
% % % % Derives a 64x64 logical dorsal-cortex mask directly from nanpxs -- the
% % % % linear indices (or logical mask) of the non-cortex ("NaN") pixels within
% % % % the full 64x64 = 4096 pixel grid -- and its boundary as a cell array of
% % % % [x y] coordinate lists (native 64x64 pixel units, one cell per connected
% % % % boundary component), smoothed via truncated Fourier descriptors.
% % % %
% % % % IMPORTANT: whether a real cortex mask exists is entirely a function of
% % % % nanpxs, NOT of P (the first dimension of motif_w). motif_w can be stored
% % % % either as [nValidPixels x K x L] (reconstructed via conditionDffMat) or
% % % % as an already-reshaped [4096 x K x L] full-grid array -- either way, if
% % % % nanpxs is supplied it still marks the true dorsal-cortex boundary within
% % % % that 64x64 grid and should be used.
% % % %
% % % % If nanpxs is empty, there is no mask information available and the whole
% % % % 64x64 frame is treated as valid (boundaryXY is returned empty).
% % % %
% % % % Boundary smoothness: rather than blurring a rasterized mask (which only
% % % % ever anti-aliases individual pixel-step edges and struggles to remove the
% % % % macroscopic staircase left by the native 64x64 resolution), the actual
% % % % pixel boundary is traced once with bwboundaries, resampled to
% % % % resampleN uniformly arc-length-spaced points (needed for a well-posed
% % % % Fourier truncation), represented as a complex sequence z = x + 1i*y, and
% % % % reconstructed from only its lowest numHarmonics frequency components via
% % % % FFT/IFFT. Truncating high frequencies of a closed, periodic curve
% % % % guarantees a result that is both smooth AND exactly seamlessly closed --
% % % % there is no "sigma in the wrong units" pitfall here, since the smoothing
% % % % is applied directly to the curve's shape, independent of any raster
% % % % resolution or supersampling choice.
% % % 
% % % boundaryXY = {};
% % % 
% % % if isempty(nanpxs)
% % %     cortexMask = true(64, 64);
% % %     return
% % % end
% % % 
% % % nanFlagVec = false(64 * 64, 1);
% % % 
% % % if islogical(nanpxs)
% % %     nanFlagVec(:) = nanpxs(:);
% % % else
% % %     nanFlagVec(nanpxs(:)) = true;
% % % end
% % % 
% % % cortexMask = reshape(~nanFlagVec, 64, 64);
% % % 
% % % if ~wantBoundary
% % %     return
% % % end
% % % 
% % % if exist('bwboundaries', 'file') ~= 2
% % %     warning('bwboundaries (Image Processing Toolbox) not found; skipping cortex boundary overlay.');
% % %     return
% % % end
% % % 
% % % rawBoundaries = bwboundaries(cortexMask, 'noholes');
% % % boundaryXY = cell(size(rawBoundaries));
% % % 
% % % for k = 1:numel(rawBoundaries)
% % %     % bwboundaries returns [row col] = [y x]; convert to [x y]
% % %     rawXY = [rawBoundaries{k}(:,2), rawBoundaries{k}(:,1)];
% % % 
% % %     resampledXY = localResampleClosedCurve(rawXY, resampleN);
% % %     boundaryXY{k} = localSmoothClosedContourFourier(resampledXY, numHarmonics);
% % % end
% % % 
% % % end
% % % 
% % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % function resampledXY = localResampleClosedCurve(xy, nResample)
% % % % localResampleClosedCurve
% % % %
% % % % Resamples a closed 2D polygon (M x 2, [x y]) to nResample points spaced
% % % % uniformly by arc length around the loop. Uniform spacing is required for
% % % % a clean/well-posed Fourier-descriptor truncation afterward.
% % % 
% % % % Ensure explicitly closed (first point repeated at the end) for arc-length
% % % % accumulation, then drop the duplicate after resampling.
% % % if norm(xy(1,:) - xy(end,:)) > 1e-9
% % %     xy = [xy; xy(1,:)];
% % % end
% % % 
% % % % Guard against zero-length (duplicate/repeated) points -- bwboundaries can
% % % % occasionally emit a repeated point where the traced path touches a thin
% % % % or degenerate pixel connection. A zero-length segment here would give
% % % % interp1 a non-strictly-increasing cumulative-distance vector, which can
% % % % produce a small localized artifact in the curve after Fourier smoothing.
% % % segLenRaw = sqrt(sum(diff(xy).^2, 2));
% % % keepPoint = [true; segLenRaw > 1e-9];
% % % xy = xy(keepPoint, :);
% % % 
% % % segLen  = sqrt(sum(diff(xy).^2, 2));
% % % cumDist = [0; cumsum(segLen)];
% % % totalLen = cumDist(end);
% % % 
% % % if totalLen == 0
% % %     resampledXY = repmat(xy(1,:), nResample, 1);
% % %     return
% % % end
% % % 
% % % targetDist = linspace(0, totalLen, nResample + 1);
% % % targetDist(end) = []; % drop duplicate closing point
% % % 
% % % xResampled = interp1(cumDist, xy(:,1), targetDist, 'linear');
% % % yResampled = interp1(cumDist, xy(:,2), targetDist, 'linear');
% % % 
% % % resampledXY = [xResampled(:), yResampled(:)];
% % % 
% % % end
% % % 
% % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % function smoothXY = localSmoothClosedContourFourier(xy, numHarmonics)
% % % % localSmoothClosedContourFourier
% % % %
% % % % Smooths a closed, uniformly arc-length-resampled 2D curve (N x 2, [x y])
% % % % by representing it as a complex sequence z = x + 1i*y, truncating its
% % % % discrete Fourier transform to the lowest numHarmonics positive and
% % % % negative frequencies (plus the DC term), and inverting. This is the
% % % % classic "Fourier descriptor" smoothing approach: it is guaranteed to
% % % % return an exactly closed curve (the representation is inherently
% % % % periodic) that is smooth by construction, with numHarmonics as the sole,
% % % % resolution-independent smoothness knob (fewer harmonics = smoother/
% % % % rounder; more harmonics = closer to the original traced shape).
% % % 
% % % N = size(xy, 1);
% % % z = complex(xy(:,1), xy(:,2));
% % % 
% % % Z = fft(z);
% % % 
% % % numHarmonics = min(numHarmonics, floor((N - 1) / 2));
% % % 
% % % keepMask = false(N, 1);
% % % keepMask(1) = true; % DC term
% % % keepMask(2:(numHarmonics + 1)) = true;               % low positive frequencies
% % % keepMask((N - numHarmonics + 1):N) = true;           % mirrored negative frequencies
% % % 
% % % Z(~keepMask) = 0;
% % % 
% % % zSmooth = ifft(Z);
% % % 
% % % smoothXY = [real(zSmooth), imag(zSmooth)];
% % % 
% % % end
% % % 
% % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % function hiresData = localUpsampleFrameData(frameData, upsampleFactor, interpMethod)
% % % % localUpsampleFrameData
% % % %
% % % % Spatially upsamples a single 64x64 frame's data for high-resolution
% % % % display. The data is interpolated as-is, with no masking applied here --
% % % % masking/fading to white is handled entirely by the shared alpha mask
% % % % (see localRasterizeSmoothMask), computed once from the same smooth
% % % % boundary curve used for the plotted line. Zeroing data outside the mask
% % % % before interpolation was tried and rejected: combined with a separately
% % % % smoothed alpha, it causes visible color fringing (a dark tinge from the
% % % % hard zero shows through wherever alpha is only partially faded).
% % % 
% % % if upsampleFactor == 1
% % %     hiresData = frameData;
% % %     return
% % % end
% % % 
% % % targetSize = size(frameData) * upsampleFactor;
% % % hiresData = imresize(frameData, targetSize, interpMethod);
% % % 
% % % end
% % % 
% % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % function hiresMask = localRasterizeSmoothMask(cortexMask, boundaryXY, upsampleFactor)
% % % % localRasterizeSmoothMask
% % % %
% % % % Builds the high-resolution alpha mask by rasterizing the SAME smooth
% % % % boundary curve(s) used for the plotted gray line (via poly2mask), rather
% % % % than independently bicubic-upsampling the original native-resolution
% % % % binary mask. Using two different smoothing methods for the line (Fourier
% % % % descriptors) and the mask (raster blurring) gives no guarantee they agree
% % % % pixel-for-pixel, which produced a visible fringe/staircase along the
% % % % edge in earlier attempts. Deriving both from the identical curve
% % % % eliminates that mismatch structurally. A light final Gaussian pass
% % % % anti-aliases the (now much finer, supersampled-grid-scale) rasterization
% % % % step, which is a much smaller and less objectionable artifact than the
% % % % native-pixel-grid blockiness this replaces.
% % % %
% % % % Falls back to a plain bicubic-upsampled binary mask if no boundary curve
% % % % is available (e.g. cortexBoundaryLogic was false, or nanpxs was empty).
% % % 
% % % targetSize = size(cortexMask) * upsampleFactor;
% % % 
% % % if isempty(boundaryXY) || exist('poly2mask', 'file') ~= 2
% % %     hiresMask = imresize(double(cortexMask), targetSize, 'bicubic');
% % %     hiresMask = min(max(hiresMask, 0), 1);
% % %     return
% % % end
% % % 
% % % maskAccum = false(targetSize);
% % % 
% % % for k = 1:numel(boundaryXY)
% % %     xy = boundaryXY{k} * upsampleFactor;
% % %     bw = poly2mask(xy(:,1), xy(:,2), targetSize(1), targetSize(2));
% % %     maskAccum = maskAccum | bw;
% % % end
% % % 
% % % hiresMask = double(maskAccum);
% % % 
% % % if exist('imgaussfilt', 'file') == 2
% % %     hiresMask = imgaussfilt(hiresMask, 1.0); % mild anti-aliasing at supersampled-grid scale
% % % end
% % % 
% % % hiresMask = min(max(hiresMask, 0), 1);
% % % 
% % % end
% % % 
% % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % function [motif_w_kept, keptFrameLabels] = localSelectFrameWindow(motif_w_sel, nFramesKeep, scope)
% % % % localSelectFrameWindow
% % % %
% % % % Selects, for each motif (or globally across all motifs), the contiguous
% % % % window of nFramesKeep frames with the highest total energy, out of the
% % % % nFramesTotal frames available. Energy for a frame is the sum of squared
% % % % pixel values across all pixels/motifs being pooled, ignoring NaNs.
% % % %
% % % % Because the chosen window is always contiguous, its complement (the
% % % % dropped frames) is always a prefix and/or suffix of the sequence -- i.e.
% % % % frames are dropped only from the edges, never from the middle, and frame
% % % % order is preserved.
% % % %
% % % % Outputs:
% % % %   motif_w_kept    : [P x nMotifSel x nFramesKeep] tensor with the
% % % %                     dropped frames removed.
% % % %   keptFrameLabels : struct with fields:
% % % %                       .global   - 1 x nFramesKeep vector of original frame
% % % %                                   indices (valid/used when scope=='global')
% % % %                       .perMotif - nMotifSel x nFramesKeep matrix of
% % % %                                   original frame indices per motif (valid/
% % % %                                   used when scope=='perMotif')
% % % 
% % % [P, nMotifSel, nFramesTotal] = size(motif_w_sel);
% % % 
% % % keptFrameLabels = struct('global', [], 'perMotif', []);
% % % 
% % % if nFramesKeep == nFramesTotal
% % %     motif_w_kept = motif_w_sel;
% % %     keptFrameLabels.global   = 1:nFramesTotal;
% % %     keptFrameLabels.perMotif = repmat(1:nFramesTotal, nMotifSel, 1);
% % %     return
% % % end
% % % 
% % % % Energy per motif per frame: [nMotifSel x nFramesTotal]
% % % sq = motif_w_sel .^ 2;
% % % sq(~isfinite(sq)) = 0;
% % % energyMotifFrame = squeeze(sum(sq, 1));       % nMotifSel x nFramesTotal
% % % if nMotifSel == 1
% % %     energyMotifFrame = reshape(energyMotifFrame, 1, nFramesTotal);
% % % end
% % % 
% % % switch scope
% % % 
% % %     case 'global'
% % % 
% % %         aggregateEnergy = sum(energyMotifFrame, 1);         % 1 x nFramesTotal
% % %         winStart = localBestWindowStart(aggregateEnergy, nFramesKeep);
% % %         keepIdx  = winStart:(winStart + nFramesKeep - 1);
% % % 
% % %         motif_w_kept = motif_w_sel(:, :, keepIdx);
% % %         keptFrameLabels.global   = keepIdx;
% % %         keptFrameLabels.perMotif = repmat(keepIdx, nMotifSel, 1);
% % % 
% % %     case 'permotif'
% % % 
% % %         motif_w_kept = zeros(P, nMotifSel, nFramesKeep, 'like', motif_w_sel);
% % %         perMotifIdx  = zeros(nMotifSel, nFramesKeep);
% % % 
% % %         for m = 1:nMotifSel
% % %             winStart = localBestWindowStart(energyMotifFrame(m, :), nFramesKeep);
% % %             keepIdx  = winStart:(winStart + nFramesKeep - 1);
% % % 
% % %             motif_w_kept(:, m, :) = motif_w_sel(:, m, keepIdx);
% % %             perMotifIdx(m, :) = keepIdx;
% % %         end
% % % 
% % %         keptFrameLabels.perMotif = perMotifIdx;
% % % 
% % %     otherwise
% % % 
% % %         error('Unknown frameSelectionScope: %s', scope);
% % % end
% % % 
% % % end
% % % 
% % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % function winStart = localBestWindowStart(energyVec, winLen)
% % % % localBestWindowStart
% % % %
% % % % Slides a window of length winLen across energyVec (1 x L) and returns the
% % % % start index of the window with the maximum total energy. Ties are broken
% % % % in favor of the earliest (smallest-index) window.
% % % 
% % % L = numel(energyVec);
% % % nWindows = L - winLen + 1;
% % % 
% % % winSums = zeros(1, nWindows);
% % % for s = 1:nWindows
% % %     winSums(s) = sum(energyVec(s:(s + winLen - 1)));
% % % end
% % % 
% % % [~, winStart] = max(winSums);
% % % 
% % % end
