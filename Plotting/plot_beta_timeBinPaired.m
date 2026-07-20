function h = plot_beta_timeBinPaired(beta, X_names, motifId, varargin)
%PLOT_BETA_TIMEBINPAIRED  Paired Go/NoGo beta weights by time bin (basis index).
%
% SYNOPSIS
%   h = plot_beta_timeBinPaired(beta, X_names, motifId, ...)
%
% DESCRIPTION
%   Visualizes the raised-cosine basis weights (beta) for one motif,
%   organized by time bin (basis index) so Go and NoGo predictors sharing
%   the same time bin can be compared directly. Two rendering styles are
%   supported via 'plotStyle':
%
%   'bar' (default):
%       Bars for Go/NoGo sharing a time bin are placed right next to each
%       other; different time bins get a small gap between them; onset
%       (filled bars) and offset (open/edge-only bars) are grouped into
%       two separate blocks with a larger gap and a dashed divider.
%       Layout: [Go1 NoGo1] gap [Go2 NoGo2] gap ... (onset block)
%                    block gap
%               [Go1 NoGo1] gap [Go2 NoGo2] gap ... (offset block)
%
%   'curve':
%       The same 9 Go and 9 NoGo coefficients per block are instead shown
%       as a smoothed curve connecting the points (shape-preserving
%       interpolation, default 'pchip' -- no overshoot beyond the data),
%       with the ACTUAL coefficients overlaid as markers. The smoothing is
%       purely a visual aid to ease comparison of the Go vs NoGo profile
%       across bins -- it does NOT reconstruct a physically-timed kernel;
%       the x-axis remains basis (time-bin) index, not seconds. Onset =
%       solid line + filled markers; offset = dashed line + open markers.
%
% INPUTS
%   beta       : [P x K] ridge coefficients (P predictors, K motifs)
%   X_names    : 1xP cellstr of predictor names (expects '..._rc<N>'
%                suffix for basis index, e.g. 'toneOnGo_rc3', as produced
%                by convolve_events_with_basis)
%   motifId    : scalar motif index to plot (column of beta)
%
% NAME-VALUE ARGS
%   'plotStyle'       : 'bar' or 'curve'. Default: 'bar'.
%   'predictorI'      : logical(1xP) mask to pre-select candidate columns
%                       (e.g. contains(X_names,'tone')). Default: all.
%   'includeOffset'   : logical, plot the toneOff block too. Default: true.
%   'pairSpacing'     : (bar style only) center-to-center distance between
%                       Go/NoGo bars within the same time bin. Default: 1.
%   'groupGap'        : (bar style only) edge-to-edge gap between adjacent
%                       time-bin pairs. Default: 0.5.
%   'blockGap'        : edge-to-edge gap between onset and offset blocks
%                       (both styles). Default: 1.5.
%   'barWidth'        : (bar style only) width of each individual bar.
%                       Default: 0.9.
%   'interpMethod'    : (curve style only) interpolation method passed to
%                       interp1 for the smoothed display curve. Default:
%                       'pchip' (shape-preserving, avoids overshoot).
%   'nInterp'         : (curve style only) number of points used to draw
%                       the smoothed curve. Default: 200.
%   'title'           : axes title. Default: auto.
%   'figureScaleFactor', 'figureWidthFactor', 'visible' : figure sizing /
%                       visibility, as in plot_beta_with_labels.m.
%   'figSaveDir', 'header', 'figSaveKeyword' : save-to-PDF options, same
%                       convention as plot_beta_with_labels.m /
%                       plot_beta_groupMeanBars.m.
%   'yLim'            : [min max] override for the y-axis. If provided,
%                       this fixed range is used for the axis limits AND
%                       for positioning the "Tone Onset"/"Tone Offset"
%                       labels and divider line -- use this instead of
%                       overriding ylim() after the fact, since the
%                       label/divider positions are baked in at plot
%                       time. Default: [] (auto-scaled to the data, as
%                       before).
%   'saveMotifId'     : override for the motif number shown in the Y-AXIS
%                       LABEL, DEFAULT TITLE, and SAVED FILENAME (data
%                       indexing/plotting itself still uses the real
%                       motifId argument). Use this when beta is
%                       synthetic/single-column data for some other real
%                       motif (e.g. a group-mean profile), so the figure
%                       is labeled and tagged with the motif it actually
%                       represents instead of its column index. Default:
%                       [] (use motifId, as before).
%
% OUTPUT
%   h : struct with fig/ax handles, plotted-object handles, computed
%       x-positions and matched basis indices for onset/offset blocks,
%       and opt (parsed args).
%
% EXAMPLE
%   predictorI = cell2mat(cellfun(@(a) contains(a,'tone'), glmRez.X_names, 'uni', 0));
%   h = plot_beta_timeBinPaired(glmRez.beta, glmRez.X_names, 12, ...
%           'predictorI', predictorI, ...
%           'plotStyle', 'curve', ...
%           'header', 'm1092_100924', ...
%           'figSaveKeyword', 'tonePredictorsPairedCurve', ...
%           'figSaveDir', figSaveDir);
%
% See also: plot_beta_with_labels, plot_beta_groupMeanBars

% -------- parse args --------
p = inputParser;
p.addParameter('plotStyle', 'bar', @(s) any(strcmpi(s, {'bar','curve'})));
p.addParameter('predictorI', [], @(x) isempty(x) || (islogical(x) && isvector(x)));
p.addParameter('includeOffset', true, @(x) islogical(x) && isscalar(x));
p.addParameter('pairSpacing', 1,   @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('groupGap',    0.5, @(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('blockGap',    1.5, @(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('barWidth',    0.9, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('interpMethod', 'pchip', @(s) any(strcmpi(s, {'linear','pchip','spline','makima'})));
p.addParameter('nInterp', 200, @(x) isnumeric(x) && isscalar(x) && x>=2);
p.addParameter('title', '', @(s) ischar(s) || isstring(s));
p.addParameter('figureScaleFactor', 2, @(x) isnumeric(x) && x>0);
p.addParameter('figureWidthFactor', 2.1, @(x) isnumeric(x) && x>0);
p.addParameter('visible', 'on', @(s) any(strcmpi(s, {'on','off'})));
p.addParameter('figSaveDir', {}, @(x) isempty(x) || ischar(x) || isstring(x) || iscell(x));
p.addParameter('header', '', @(s) ischar(s) || isstring(s));
p.addParameter('figSaveKeyword', '', @(s) ischar(s) || isstring(s));
p.addParameter('yLim', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
p.addParameter('saveMotifId', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
p.parse(varargin{:});
opt = p.Results;

% Single source of truth for "which motif number to show/save this as" --
% used by BOTH the y-axis label/default title AND the saved filename.
% Falls back to motifId (the real column-index argument) unless overridden.
motifForDisplay = motifId;
if ~isempty(opt.saveMotifId)
    motifForDisplay = opt.saveMotifId;
end

% normalize figSaveDir
figSaveDir = opt.figSaveDir;
if iscell(figSaveDir)
    if isempty(figSaveDir), figSaveDir = ""; else, figSaveDir = string(figSaveDir{1}); end
else
    figSaveDir = string(figSaveDir);
end
header         = string(opt.header);
figSaveKeyword = string(opt.figSaveKeyword);

if strcmpi(opt.plotStyle,'bar') && opt.pairSpacing <= opt.barWidth
    warning('plot_beta_timeBinPaired:tightPairSpacing', ...
        'pairSpacing (%.2f) <= barWidth (%.2f): Go/NoGo bars within a time bin may overlap.', ...
        opt.pairSpacing, opt.barWidth);
end

% -------- sanity checks --------
[P, K] = size(beta);
assert(iscellstr(X_names) && numel(X_names)==P, 'X_names must be 1xP cellstr.');
assert(isscalar(motifId) && motifId>=1 && motifId<=K, 'motifId out of range.');

if ~isempty(opt.predictorI)
    sel = opt.predictorI(:).';
    assert(numel(sel)==P, 'predictorI must be logical(1xP).');
    assert(any(sel), 'predictorI selects zero predictors.');
else
    sel = true(1,P);
end

bAll     = beta(:, motifId);
bSel     = bAll(sel);
namesSel = X_names(sel);

% -------- classify predictors: type + basis (time-bin) index --------
[typeKey, idxNum] = classify_tone_predictor_(namesSel);

isOnG  = typeKey=="tOnG";   isOnNg  = typeKey=="tOnNg";
isOffG = typeKey=="tOffG";  isOffNg = typeKey=="tOffNg";

idxOnG  = idxNum(isOnG);   valOnG  = bSel(isOnG);
idxOnNg = idxNum(isOnNg);  valOnNg = bSel(isOnNg);
idxOffG = idxNum(isOffG);  valOffG = bSel(isOffG);
idxOffNg= idxNum(isOffNg); valOffNg= bSel(isOffNg);

assert(~isempty(idxOnG) && ~isempty(idxOnNg), ...
    'No toneOnGo/toneOnNoGo predictors found -- check predictorI/X_names.');

% drop any unmatched (NaN index) entries defensively
keepOnG  = ~isnan(idxOnG);   idxOnG  = idxOnG(keepOnG);   valOnG  = valOnG(keepOnG);
keepOnNg = ~isnan(idxOnNg);  idxOnNg = idxOnNg(keepOnNg); valOnNg = valOnNg(keepOnNg);

% sort each by basis index so Go/NoGo align by time bin
[idxOnG,  ordOnG]  = sort(idxOnG);   valOnG  = valOnG(ordOnG);
[idxOnNg, ordOnNg] = sort(idxOnNg);  valOnNg = valOnNg(ordOnNg);
assert(isequal(idxOnG(:), idxOnNg(:)), ...
    'toneOnGo and toneOnNoGo basis indices do not match -- cannot pair by time bin.');
nBinsOn = numel(idxOnG);

doOffset = opt.includeOffset && ~isempty(idxOffG) && ~isempty(idxOffNg);
if doOffset
    keepOffG  = ~isnan(idxOffG);  idxOffG  = idxOffG(keepOffG);  valOffG  = valOffG(keepOffG);
    keepOffNg = ~isnan(idxOffNg); idxOffNg = idxOffNg(keepOffNg);valOffNg = valOffNg(keepOffNg);
    [idxOffG,  ordOffG]  = sort(idxOffG);   valOffG  = valOffG(ordOffG);
    [idxOffNg, ordOffNg] = sort(idxOffNg);  valOffNg = valOffNg(ordOffNg);
    assert(isequal(idxOffG(:), idxOffNg(:)), ...
        'toneOffGo and toneOffNoGo basis indices do not match -- cannot pair by time bin.');
    nBinsOff = numel(idxOffG);
end

% -------- shared style handles --------
blueish = [0.30 0.50 0.90];
redish  = [0.90 0.30 0.30];
edgeOnlyLineWidth = 2;

% -------- figure/axes setup --------
h = struct();
h.fig = figure('Color','w', 'Visible', opt.visible);
set(h.fig, 'Units', 'normalized');
pos = get(h.fig, 'Position');
pos(3:4) = pos(3:4) * opt.figureScaleFactor;
set(h.fig, 'Position', pos);
h.ax = axes('Parent', h.fig); hold(h.ax, 'on');

switch lower(opt.plotStyle)

    %% ===================== BAR STYLE =====================
    case 'bar'
        unit     = opt.pairSpacing;
        gapPair  = opt.groupGap;
        gapBlock = opt.blockGap;
        bw       = opt.barWidth;

        xOnG  = nan(nBinsOn,1);
        xOnNg = nan(nBinsOn,1);
        cursor = 1;
        for b = 1:nBinsOn
            xOnG(b)  = cursor;
            xOnNg(b) = cursor + unit;
            cursor   = xOnNg(b) + bw + gapPair;
        end
        lastOnsetX = xOnNg(end);

        if doOffset
            xOffG  = nan(nBinsOff,1);
            xOffNg = nan(nBinsOff,1);
            cursor = lastOnsetX + bw + gapBlock;
            for b = 1:nBinsOff
                xOffG(b)  = cursor;
                xOffNg(b) = cursor + unit;
                cursor    = xOffNg(b) + bw + gapPair;
            end
        end

        h.barsOnG  = gobjects(nBinsOn,1);
        h.barsOnNg = gobjects(nBinsOn,1);
        for b = 1:nBinsOn
            h.barsOnG(b)  = bar(h.ax, xOnG(b),  valOnG(b),  bw, 'FaceColor', blueish, 'EdgeColor', 'none');
            h.barsOnNg(b) = bar(h.ax, xOnNg(b), valOnNg(b), bw, 'FaceColor', redish,  'EdgeColor', 'none');
        end
        if doOffset
            h.barsOffG  = gobjects(nBinsOff,1);
            h.barsOffNg = gobjects(nBinsOff,1);
            for b = 1:nBinsOff
                h.barsOffG(b)  = bar(h.ax, xOffG(b),  valOffG(b),  bw, 'FaceColor', 'none', 'EdgeColor', blueish, 'LineWidth', edgeOnlyLineWidth);
                h.barsOffNg(b) = bar(h.ax, xOffNg(b), valOffNg(b), bw, 'FaceColor', 'none', 'EdgeColor', redish,  'LineWidth', edgeOnlyLineWidth);
            end
        end

        if doOffset
            xAll = [xOnG; xOnNg; xOffG; xOffNg];
        else
            xAll = [xOnG; xOnNg];
        end
        xPad = unit;
        xlim(h.ax, [min(xAll)-xPad, max(xAll)+xPad]);   % locks XLimMode='manual'
        if ~isempty(opt.yLim)
            yl = opt.yLim;
        else
            yl = ylim(h.ax);
        end
        ylim(h.ax, yl);                                  % locks YLimMode='manual'

        plot(h.ax, xlim(h.ax), [0 0], 'k-', 'LineWidth', 0.8);

        xtickPosOn = (xOnG + xOnNg) / 2;
        xtickLblOn = arrayfun(@(v) sprintf('%d', v), idxOnG, 'uni', 0);
        if nBinsOn > 4
            for b = 1:numel(xtickLblOn)
                if mod(idxOnG(b),2)==0, xtickLblOn{b} = ''; end
            end
        end
        if doOffset
            xtickPosOff = (xOffG + xOffNg) / 2;
            xtickLblOff = arrayfun(@(v) sprintf('%d', v), idxOffG, 'uni', 0);
            if nBinsOff > 4
                for b = 1:numel(xtickLblOff)
                    if mod(idxOffG(b),2)==0, xtickLblOff{b} = ''; end
                end
            end
            xtickPos = [xtickPosOn; xtickPosOff];
            xtickLbl = [xtickLblOn; xtickLblOff];
        else
            xtickPos = xtickPosOn;
            xtickLbl = xtickLblOn;
        end

        xticks(h.ax, xtickPos);
        xticklabels(h.ax, xtickLbl);
        set(h.ax, 'TickLabelInterpreter', 'none');

        yTop = yl(2) - 0.05*(yl(2)-yl(1));
        text(h.ax, mean(xtickPosOn), yTop, 'Tone Onset', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        if doOffset
            xDivider = (lastOnsetX + xOffG(1)) / 2;
            plot(h.ax, [xDivider xDivider], yl, 'k--', 'LineWidth', 0.5);
            text(h.ax, mean(xtickPosOff), yTop, 'Tone Offset', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        end

        % legend via fixed-xlim off-screen proxies (keeps real bars untouched)
        xProxy = min(xAll) - 10*xPad;
        pGoOn   = bar(h.ax, xProxy, 0, 'FaceColor', blueish, 'EdgeColor', 'none');
        pNogoOn = bar(h.ax, xProxy, 0, 'FaceColor', redish,  'EdgeColor', 'none');
        legH   = [pGoOn, pNogoOn];
        legLbl = {'Go (Onset)', 'NoGo (Onset)'};
        if doOffset
            pGoOff   = bar(h.ax, xProxy, 0, 'FaceColor', 'none', 'EdgeColor', blueish, 'LineWidth', edgeOnlyLineWidth);
            pNogoOff = bar(h.ax, xProxy, 0, 'FaceColor', 'none', 'EdgeColor', redish,  'LineWidth', edgeOnlyLineWidth);
            legH   = [legH, pGoOff, pNogoOff];
            legLbl = [legLbl, {'Go (Offset)', 'NoGo (Offset)'}];
        end
        legend(h.ax, legH, legLbl, 'Location', 'bestoutside');

        h.xOnG = xOnG; h.xOnNg = xOnNg; h.idxOnG = idxOnG;
        if doOffset
            h.xOffG = xOffG; h.xOffNg = xOffNg; h.idxOffG = idxOffG;
        end

    %% ===================== CURVE STYLE =====================
    case 'curve'
        markerSize = 6;
        lineWidthCurve = 2;
        im = opt.interpMethod;

        % onset block: natural bin-index x positions (Go/NoGo share x)
        xOn = idxOnG(:);
        if numel(xOn) >= 2
            xOnDense    = linspace(xOn(1), xOn(end), opt.nInterp);
            curveGoOn   = interp1(xOn, valOnG,  xOnDense, im);
            curveNogoOn = interp1(xOn, valOnNg, xOnDense, im);
        else
            xOnDense = xOn; curveGoOn = valOnG; curveNogoOn = valOnNg;
        end

        if doOffset
            % offset block: shift x so it sits after onset block + blockGap,
            % preserving unit-1 spacing between consecutive bins
            xOff = idxOffG(:) - idxOffG(1) + xOn(end) + opt.blockGap + 1;
            if numel(xOff) >= 2
                xOffDense     = linspace(xOff(1), xOff(end), opt.nInterp);
                curveGoOff    = interp1(xOff, valOffG,  xOffDense, im);
                curveNogoOff  = interp1(xOff, valOffNg, xOffDense, im);
            else
                xOffDense = xOff; curveGoOff = valOffG; curveNogoOff = valOffNg;
            end
        end

        % onset: solid line, filled markers (actual coefficients)
        h.curveGoOn   = plot(h.ax, xOnDense, curveGoOn,   '-', 'Color', blueish, 'LineWidth', lineWidthCurve);
        plot(h.ax, xOn, valOnG,  'o', 'MarkerFaceColor', blueish, 'MarkerEdgeColor', blueish, 'MarkerSize', markerSize, 'LineStyle', 'none');
        h.curveNogoOn = plot(h.ax, xOnDense, curveNogoOn, '-', 'Color', redish,  'LineWidth', lineWidthCurve);
        plot(h.ax, xOn, valOnNg, 'o', 'MarkerFaceColor', redish,  'MarkerEdgeColor', redish,  'MarkerSize', markerSize, 'LineStyle', 'none');

        if doOffset
            % offset: dashed line, open markers (actual coefficients)
            h.curveGoOff   = plot(h.ax, xOffDense, curveGoOff,   '--', 'Color', blueish, 'LineWidth', lineWidthCurve);
            plot(h.ax, xOff, valOffG,  'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', blueish, 'MarkerSize', markerSize, 'LineWidth', edgeOnlyLineWidth, 'LineStyle', 'none');
            h.curveNogoOff = plot(h.ax, xOffDense, curveNogoOff, '--', 'Color', redish,  'LineWidth', lineWidthCurve);
            plot(h.ax, xOff, valOffNg, 'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', redish,  'MarkerSize', markerSize, 'LineWidth', edgeOnlyLineWidth, 'LineStyle', 'none');
        end

        if doOffset
            xAll = [xOn; xOff];
        else
            xAll = xOn;
        end
        xPad = 0.5;
        xlim(h.ax, [min(xAll)-xPad, max(xAll)+xPad]);
        if ~isempty(opt.yLim)
            yl = opt.yLim;
        else
            yl = ylim(h.ax);
        end
        ylim(h.ax, yl);

        plot(h.ax, xlim(h.ax), [0 0], 'k-', 'LineWidth', 0.8);

        xtickLblOn = arrayfun(@(v) sprintf('%d', v), idxOnG, 'uni', 0);
        if nBinsOn > 4
            for b = 1:numel(xtickLblOn)
                if mod(idxOnG(b),2)==0, xtickLblOn{b} = ''; end
            end
        end
        if doOffset
            xtickLblOff = arrayfun(@(v) sprintf('%d', v), idxOffG, 'uni', 0);
            if nBinsOff > 4
                for b = 1:numel(xtickLblOff)
                    if mod(idxOffG(b),2)==0, xtickLblOff{b} = ''; end
                end
            end
            xtickPos = [xOn; xOff];
            xtickLbl = [xtickLblOn; xtickLblOff];
        else
            xtickPos = xOn;
            xtickLbl = xtickLblOn;
        end

        xticks(h.ax, xtickPos);
        xticklabels(h.ax, xtickLbl);
        set(h.ax, 'TickLabelInterpreter', 'none');

        yTop = yl(2) - 0.05*(yl(2)-yl(1));
        text(h.ax, mean(xOn), yTop, 'Tone Onset', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        if doOffset
            xDivider = (xOn(end) + xOff(1)) / 2;
            plot(h.ax, [xDivider xDivider], yl, 'k--', 'LineWidth', 0.5);
            text(h.ax, mean(xOff), yTop, 'Tone Offset', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        end

        legH   = [h.curveGoOn, h.curveNogoOn];
        legLbl = {'Go (Onset)', 'NoGo (Onset)'};
        if doOffset
            legH   = [legH, h.curveGoOff, h.curveNogoOff];
            legLbl = [legLbl, {'Go (Offset)', 'NoGo (Offset)'}];
        end
        legend(h.ax, legH, legLbl, 'Location', 'bestoutside');

        h.xOn = xOn; h.idxOnG = idxOnG;
        if doOffset
            h.xOff = xOff; h.idxOffG = idxOffG;
        end
end

xlabel(h.ax, 'Time bin (basis index)');
ylabel(h.ax, sprintf('\\beta (motif %d)', motifForDisplay));

if strlength(string(opt.title)) > 0
    title(h.ax, opt.title, 'Interpreter', 'none');
else
    styleStr = ternary_(strcmpi(opt.plotStyle,'curve'), 'smoothed \beta profile', 'paired \beta');
    title(h.ax, sprintf('Motif %d: %s by time bin', motifForDisplay, styleStr), 'Interpreter', 'tex');
end
box(h.ax, 'off');
grid(h.ax, 'on');
set(h.ax, 'TickDir', 'out');

r = pbaspect;
pbaspect([opt.figureWidthFactor*r(1) r(2) r(3)]);

% -------- save figure (optional) --------
if strlength(figSaveDir) > 0
    if ~isfolder(figSaveDir)
        mkdir(figSaveDir);
    end
    dateStr = char(datetime("today","Format","MMddyy"));

    parts = strings(0,1);
    if strlength(header) > 0,         parts(end+1,1) = header; end
    if strcmpi(opt.plotStyle, 'curve')
        parts(end+1,1) = "glmBetaCurve";
    else
        parts(end+1,1) = "glmBetaPaired";
    end
    if strlength(figSaveKeyword) > 0, parts(end+1,1) = figSaveKeyword; end
    parts(end+1,1) = "motif" + string(motifForDisplay);
    parts(end+1,1) = dateStr;

    figSaveName = strjoin(parts, "_");
    print(h.fig, fullfile(figSaveDir, figSaveName), '-dpdf', '-painters', '-bestfit');
end

h.opt = opt;

end % function


% ===== helper: classify predictor by tone type + basis (time-bin) index =====
function [typeKey, idxNum] = classify_tone_predictor_(nameC)
n = numel(nameC);
typeKey = repmat("other", n, 1);
idxNum  = nan(n,1);

for i = 1:n
    nm = string(nameC{i});

    if contains(nm, "toneOnGo", "IgnoreCase", true)
        typeKey(i) = "tOnG";
    elseif contains(nm, "toneOffGo", "IgnoreCase", true)
        typeKey(i) = "tOffG";
    elseif contains(nm, "toneOnNoGo", "IgnoreCase", true)
        typeKey(i) = "tOnNg";
    elseif contains(nm, "toneOffNoGo", "IgnoreCase", true)
        typeKey(i) = "tOffNg";
    end

    tok = regexp(nm, '_rc(\d+)$', 'tokens', 'once');
    if ~isempty(tok)
        idxNum(i) = str2double(tok{1});
    end
end
end

% ===== helper: tiny inline ternary for the title string =====
function out = ternary_(cond, a, b)
if cond, out = a; else, out = b; end
end