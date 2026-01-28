function h = plotMDSConvergenceAndAmongDistance(Y, sessInfo, ell, cmap, mIdC, varargin)
% plotMDSConvergenceAndAmongDistance
%   Plot per-mouse convergence to a reference ellipsoid center (Mahalanobis or Euclidean),
%   with optional "among-mice" dispersion at each (aligned) session index.
%
%   IMPORTANT (color + ordering is now STRICTLY driven by mIdC):
%     - mIdC is REQUIRED and defines the global mouse ordering.
%     - cmap is interpreted in the SAME row order as mIdC.
%     - selectMice only subsets AFTER mIdC ordering (stable).
%
% INPUTS
%   Y        : [S x dim] MDS coordinates
%   sessInfo : table with S rows (same order as Y). Needs variables:
%              - mouseId (or header for inference) and sessWithin (or inferred)
%   ell      : struct with fields .mu (1xdim) and .Sigma (dimxdim)
%   cmap     : [] OR [numel(mIdC) x 3] RGB colormap (row i corresponds to mIdC{i})
%   mIdC     : REQUIRED cellstr/string of mouse IDs, e.g. {'m1237','m1092',...}
%
% NAME–VALUE PAIRS
%   'dims'            : dims of Y to use (default: 1:size(Y,2))
%   'selectMice'      : mouse IDs to include (default: {}) (subset of mIdC)
%   'alignment'       : 'left' (default) or 'right'
%
%   'plotAmong'       : true/false (default: false)
%   'amongMetric'     : 'euclidean' (default) or 'mahal'
%   'amongMahSource'  : 'global' (default) or 'ellipsoid' (used if amongMetric='mahal')
%   'amongMahReg'     : regularization for SigmaAmong (default: 1e-10)
%   'amongNormalize'  : 'points' (default) or 'pairs'
%   'amongMode'       : 'overlay' (default) or 'only'
%
%   'refMode'         : 'mahal' (default) or 'euclidean'
%
%   'highlightInside' : true/false (default: true)
%   'confLevel'       : chi2 confidence level (default: 0.95)
%   'insideEdgeColor' : default 'k'
%   'insideLineWidth' : default 1.5
%
%   'ax'              : axis handle to plot into (default: [])
%   'titleStr'        : title string for convergence panel (default: '')
%   'lineWidth'       : per-mouse line width (default: 1.8)
%   'markerSize'      : marker size (default: 5)
%   'gridOn'          : true/false (default: true)
%   'figWidthScale'   : scale factor for figure width (default: 1.5)
%
% NEW (saving)
%   'figSaveDir'      : directory to save PDF (default: '')
%   'reprint'         : overwrite if exists (default: false)
%   'figureNameBase'  : filename prefix (default: 'plotXcorrPosLagMDSConvergence_')
%   'trialType'       : 'both' | 'go' | 'nogo' (default: 'both')
%
%   Filename:
%     <figSaveDir>/<figureNameBase><trialType>_<MMDDYY>.pdf
%
% OUTPUT
%   h.fig
%   h.ax.convergence, h.ax.among
%   h.lines.mouse, h.lines.among
%   h.data (computed vectors, options, membership flags)
%   h.save (save info)

% -------------------- parse options --------------------
p = inputParser;

p.addParameter('dims', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
p.addParameter('selectMice', {}, @(x) iscell(x) || isstring(x));
p.addParameter('alignment', 'left', @(s) ischar(s) || isstring(s));

p.addParameter('plotAmong', false, @(x) islogical(x) && isscalar(x));
p.addParameter('amongMetric', 'euclidean', @(s) ischar(s) || isstring(s));
p.addParameter('amongMahSource', 'global', @(s) ischar(s) || isstring(s));
p.addParameter('amongMahReg', 1e-10, @(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('amongNormalize', 'points', @(s) ischar(s) || isstring(s));
p.addParameter('amongMode', 'overlay', @(s) ischar(s) || isstring(s));

p.addParameter('refMode', 'mahal', @(s) ischar(s) || isstring(s));

% highlight options
p.addParameter('highlightInside', true, @(x) islogical(x) && isscalar(x));
p.addParameter('confLevel', 0.95, @(x) isnumeric(x) && isscalar(x) && x>0 && x<1);
p.addParameter('insideEdgeColor', 'k');
p.addParameter('insideLineWidth', 1.5, @(x) isnumeric(x) && isscalar(x) && x>=0);

p.addParameter('ax', [], @(x) isempty(x) || isgraphics(x,'axes'));
p.addParameter('titleStr', '', @(s) ischar(s) || isstring(s));
p.addParameter('lineWidth', 1.8, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('markerSize', 5, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('gridOn', true, @(x) islogical(x) && isscalar(x));
p.addParameter('figWidthScale', 1.5, @(x) isnumeric(x) && isscalar(x) && x>0);

% -------- NEW: saving --------
p.addParameter('figSaveDir', '', @(s) ischar(s) || isstring(s));
p.addParameter('reprint', false, @(x) islogical(x) && isscalar(x));
p.addParameter('figureNameBase', 'plotXcorrPosLagMDSConvergence_', @(s) ischar(s) || isstring(s));
p.addParameter('trialType', 'both', @(s) any(strcmpi(string(s), ["both","go","nogo"])));

p.parse(varargin{:});
opt = p.Results;

alignment      = lower(string(opt.alignment));
amongMetric    = lower(string(opt.amongMetric));
amongMahSource = lower(string(opt.amongMahSource));
amongNormalize = lower(string(opt.amongNormalize));
amongMode      = lower(string(opt.amongMode));
refMode        = lower(string(opt.refMode));
trialType      = lower(string(opt.trialType));

if ~ismember(alignment, ["left","right"])
    error('alignment must be ''left'' or ''right''.');
end
if ~ismember(amongMetric, ["euclidean","mahal"])
    error('amongMetric must be ''euclidean'' or ''mahal''.');
end
if ~ismember(amongMahSource, ["global","ellipsoid"])
    error('amongMahSource must be ''global'' or ''ellipsoid''.');
end
if ~ismember(amongNormalize, ["points","pairs"])
    error('amongNormalize must be ''points'' or ''pairs''.');
end
if ~ismember(amongMode, ["overlay","only"])
    error('amongMode must be ''overlay'' or ''only''.');
end
if ~ismember(refMode, ["mahal","euclidean"])
    error('refMode must be ''mahal'' or ''euclidean''.');
end
if ~ismember(trialType, ["both","go","nogo"])
    error('trialType must be ''both'', ''go'', or ''nogo''.');
end

wantSideBySide = opt.plotAmong && amongMode == "overlay";

% -------------------- validate + normalize mIdC --------------------
if nargin < 5 || isempty(mIdC)
    error('mIdC is REQUIRED (cell array / string array of mouse IDs).');
end
mIdC = string(mIdC(:));
if any(strlength(mIdC)==0)
    error('mIdC contains empty mouse IDs.');
end

% -------------------- dims selection --------------------
if isempty(opt.dims)
    dims = 1:size(Y,2);
else
    dims = opt.dims(:)';
    if any(dims > size(Y,2))
        error('dims exceeds available columns in Y.');
    end
end
Yd = Y(:, dims);
d = numel(dims);

% -------------------- infer mouseId + sessWithin --------------------
mouseId     = local_inferMouseId(sessInfo);
sessWithin  = local_inferSessWithin(sessInfo, mouseId);

% -------------------- strict mouse universe based on mIdC --------------------
presentInSess = ismember(mIdC, unique(mouseId, 'stable'));
mouseU_all = mIdC(presentInSess); % mIdC order (filtered)

if isempty(mouseU_all)
    error('None of the mIdC mice were found in sessInfo.mouseId.');
end

% -------------------- selectMice filtering (preserve mIdC order) --------------------
selectMice = string(opt.selectMice);
if ~isempty(selectMice)
    keepMouse = ismember(mouseU_all, selectMice);
    if ~any(keepMouse)
        warning('selectMice did not match any of the available mice (mIdC ∩ sessInfo). Nothing plotted.');
        h = struct();
        return;
    end
    mouseU = mouseU_all(keepMouse);
else
    mouseU = mouseU_all;
end

% -------------------- strict colormap mapping by mIdC --------------------
if isempty(cmap)
    cmap_all = lines(numel(mIdC));
else
    if size(cmap,2) ~= 3
        error('cmap must be N x 3.');
    end
    if size(cmap,1) ~= numel(mIdC)
        error('cmap must have exactly numel(mIdC) rows (strict mapping).');
    end
    cmap_all = cmap;
end

cmap_avail = cmap_all(presentInSess, :);
cmap_plot  = cmap_avail(ismember(mouseU_all, mouseU), :);

% -------------------- reference (ell) sliced to dims --------------------
if isempty(ell) || ~isfield(ell,'mu') || ~isfield(ell,'Sigma')
    error('ell must contain fields mu and Sigma.');
end
mu_full  = ell.mu(:)';
Sig_full = ell.Sigma;

if numel(mu_full) < max(dims) || any(size(Sig_full) < max(dims))
    error('ell.mu / ell.Sigma do not have enough dimensions for requested dims.');
end

mu    = mu_full(dims);
Sigma = Sig_full(dims, dims);
Sigma = (Sigma + Sigma')/2 + 1e-10*eye(d);
R_ell = local_cholStable(Sigma);

% -------------------- build aligned x per session --------------------
maxFinal = 0;
perMouseMax = zeros(numel(mouseU),1);
for im = 1:numel(mouseU)
    idx = find(mouseId == mouseU(im));
    perMouseMax(im) = max(sessWithin(idx));
    maxFinal = max(maxFinal, perMouseMax(im));
end

xAligned = sessWithin;
if alignment == "right"
    for im = 1:numel(mouseU)
        idx = find(mouseId == mouseU(im));
        shift = maxFinal - perMouseMax(im);
        xAligned(idx) = sessWithin(idx) + shift;
    end
end

% -------------------- compute per-session distance to ref --------------------
S = size(Yd,1);
distToRef = nan(S,1);

switch refMode
    case "euclidean"
        diff = Yd - mu;
        distToRef = sqrt(sum(diff.^2, 2));
    case "mahal"
        X = Yd - mu;
        Z = X / R_ell;
        distToRef = sqrt(sum(Z.^2, 2));
end

% -------------------- inside-ellipsoid membership --------------------
insideRef = false(S,1);
if opt.highlightInside
    X = Yd - mu;
    Z = X / R_ell;
    md2 = sum(Z.^2, 2);
    thr = chi2inv(opt.confLevel, d);
    insideRef = (md2 <= thr);
end

% -------------------- Mahalanobis setup for among-distance (optional) --------------------
W_among = [];
if opt.plotAmong && amongMetric == "mahal"
    switch amongMahSource
        case "global"
            SigmaAmong = cov(Yd, 'omitrows');
        case "ellipsoid"
            SigmaAmong = Sigma;
    end
    SigmaAmong = (SigmaAmong + SigmaAmong')/2 + opt.amongMahReg*eye(size(SigmaAmong,1));
    W_among = local_cholStable(SigmaAmong);
end

% -------------------- figure/axes setup (+ width scaling) --------------------
userAx = opt.ax;

if wantSideBySide
    h.fig = figure('Color','w');
else
    if isempty(userAx)
        h.fig = figure('Color','w');
    else
        h.fig = ancestor(userAx, 'figure');
    end
end

pos = h.fig.Position;
pos(3) = pos(3) * opt.figWidthScale;
h.fig.Position = pos;

h.ax = struct();

if wantSideBySide
    tl = tiledlayout(h.fig, 1, 2, 'TileSpacing','compact', 'Padding','compact');
    h.ax.convergence = nexttile(tl, 1);
    h.ax.among       = nexttile(tl, 2);
else
    if amongMode == "only"
        if isempty(userAx)
            h.ax.among = axes('Parent', h.fig);
        else
            h.ax.among = userAx;
        end
        h.ax.convergence = [];
    else
        if isempty(userAx)
            h.ax.convergence = axes('Parent', h.fig);
        else
            h.ax.convergence = userAx;
        end
        h.ax.among = [];
    end
end

% -------------------- plotting --------------------
h.lines = struct();
h.lines.mouse = gobjects(numel(mouseU),1);
h.lines.among = gobjects(1);

% --- (A) Convergence plot ---
if ~isempty(h.ax.convergence)
    axC = h.ax.convergence;
    hold(axC,'on');

    for im = 1:numel(mouseU)
        mid = mouseU(im);
        idx = find(mouseId == mid);

        [~, ord] = sortrows([xAligned(idx), sessWithin(idx)], [1 2]);
        idx = idx(ord);

        x = xAligned(idx);
        y = distToRef(idx);

        hl = plot(axC, x, y, '-', 'LineWidth', opt.lineWidth, 'Color', cmap_plot(im,:));
        h.lines.mouse(im) = hl;

        hp = plot(axC, x, y, 'o', ...
            'MarkerSize', opt.markerSize, ...
            'Color', cmap_plot(im,:), ...
            'MarkerFaceColor', cmap_plot(im,:), ...
            'MarkerEdgeColor', 'none');
        hp.HandleVisibility = 'off';

        if opt.highlightInside && any(insideRef(idx))
            idxIn = idx(insideRef(idx));
            hIn = plot(axC, xAligned(idxIn), distToRef(idxIn), 'o', ...
                'MarkerSize', opt.markerSize, ...
                'MarkerFaceColor', cmap_plot(im,:), ...
                'MarkerEdgeColor', opt.insideEdgeColor, ...
                'LineWidth', opt.insideLineWidth);
            hIn.HandleVisibility = 'off';
        end
    end

    xlabel(axC, sprintf('sessWithin (%s aligned)', alignment));
    ylabel(axC, sprintf('Distance to reference (%s)', refMode));

    if strlength(string(opt.titleStr)) > 0
        title(axC, opt.titleStr, 'Interpreter','none');
    end
    if opt.gridOn, grid(axC,'on'); end
    box(axC,'off');

    legH  = h.lines.mouse(isgraphics(h.lines.mouse));
    legLb = cellstr(mouseU);
    legend(axC, legH, legLb, 'Location','eastoutside');

    hold(axC,'off');
end

% --- (B) Among-mice dispersion plot ---
amongX = [];
amongY = [];

if opt.plotAmong
    if wantSideBySide
        axA = h.ax.among;
    else
        axA = h.ax.among;
        if isempty(axA) && amongMode == "overlay"
            axA = h.ax.convergence;
        end
    end

    if ~isempty(axA)
        hold(axA,'on');

        xGrid = unique(xAligned, 'sorted');
        amongY = nan(numel(xGrid),1);

        for xi = 1:numel(xGrid)
            x0 = xGrid(xi);

            pts = [];
            for im = 1:numel(mouseU)
                idx = find(mouseId == mouseU(im) & xAligned == x0);
                if isempty(idx), continue; end
                [~, kk] = max(sessWithin(idx));
                idx = idx(kk);
                pts = [pts; Yd(idx,:)]; %#ok<AGROW>
            end

            nPts = size(pts,1);
            if nPts < 2
                amongY(xi) = nan;
                continue;
            end

            switch amongMetric
                case "euclidean"
                    dists = pdist(pts, 'euclidean');
                case "mahal"
                    ptsW  = pts / W_among;
                    dists = pdist(ptsW, 'euclidean');
            end

            m = mean(dists);
            if amongNormalize == "points"
                m = m / nPts;
            end
            amongY(xi) = m;
        end

        amongX = xGrid;

        h.lines.among = plot(axA, amongX, amongY, '-', 'LineWidth', 2.5, 'Color', [0 0 0]);

        xlabel(axA, sprintf('sessWithin (%s aligned)', alignment));
        if amongMetric == "mahal"
            ylabel(axA, sprintf('Among dispersion (%s:%s/%s)', amongMetric, amongMahSource, amongNormalize));
        else
            ylabel(axA, sprintf('Among dispersion (%s/%s)', amongMetric, amongNormalize));
        end

        if opt.gridOn, grid(axA,'on'); end
        box(axA,'off');

        hold(axA,'off');
    end
end

% -------------------- package outputs --------------------
h.data = struct();
h.data.mIdC          = mIdC;
h.data.presentInSess = presentInSess;
h.data.mouseU_all    = mouseU_all;
h.data.mouseU        = mouseU;
h.data.mouseId       = mouseId;
h.data.sessWithin    = sessWithin;
h.data.xAligned      = xAligned;
h.data.dims          = dims;
h.data.distToRef     = distToRef;
h.data.insideRef     = insideRef;
h.data.amongX        = amongX;
h.data.amongY        = amongY;
h.data.cmap_all      = cmap_all;
h.data.cmap_plot     = cmap_plot;
h.data.params        = opt;
h.data.params.alignment      = alignment;
h.data.params.amongMetric    = amongMetric;
h.data.params.amongMahSource = amongMahSource;
h.data.params.amongNormalize = amongNormalize;
h.data.params.amongMode      = amongMode;
h.data.params.refMode        = refMode;
h.data.params.trialType      = trialType;

% -------------------- NEW: saving --------------------
h.save = struct('didSave', false, 'file', '', 'figSaveDir', '', 'base', '', 'trialType', '', 'dateStr', '');

figSaveDir = strtrim(string(opt.figSaveDir));
if strlength(figSaveDir) > 0
    figSaveDir = string(figSaveDir); %#ok<NASGU>
end

if strlength(strtrim(string(opt.figSaveDir))) > 0
    % user provided figSaveDir -> eligible for saving
    % but only save if they explicitly requested by passing figSaveDir? (no printFig flag requested)
    % We'll save whenever figSaveDir is non-empty (as requested) and figure exists.
    outDir = char(string(opt.figSaveDir));
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    base = char(string(opt.figureNameBase));
    dstr = datestr(now, 'mmddyy');

    % requested full filename: plotMDSConvergenceAndAmongDistance_(trialType)_(date).pdf
    % (use user's base unless you want exactly fixed string; you asked baseNameBase FYI)
    fnBase = sprintf('%s%s_%s', base, trialType, dstr);
    outFile = fullfile(outDir, [fnBase '.pdf']);

    if exist(outFile, 'file') && ~opt.reprint
        % skip
    else
        set(h.fig, 'InvertHardcopy', 'off');
        print(h.fig, outFile, '-dpdf', '-painters', '-bestfit');
        h.save.didSave  = true;
        h.save.file     = outFile;
        h.save.figSaveDir = outDir;
        h.save.base     = fnBase;
        h.save.trialType = char(trialType);
        h.save.dateStr  = dstr;
    end
end

end

% ========================= local helpers =========================
function mouseId = local_inferMouseId(sessInfo)
S = height(sessInfo);

if ismember('mouseId', sessInfo.Properties.VariableNames)
    mouseId = string(sessInfo.mouseId);
    if all(strlength(mouseId) > 0)
        return;
    end
end

if ismember('header', sessInfo.Properties.VariableNames)
    hdr = string(sessInfo.header);
    mouseId = strings(S,1);
    for ii = 1:S
        tok = regexp(hdr(ii), '(m\d{3,5})', 'tokens', 'once');
        if isempty(tok)
            mouseId(ii) = "mouse_" + string(ii);
        else
            mouseId(ii) = string(tok{1});
        end
    end
    return;
end

if ismember('mouseIdx', sessInfo.Properties.VariableNames)
    mouseId = "mouse_" + string(sessInfo.mouseIdx);
    return;
end

mouseId = "mouse_" + string((1:S)');
end

function sessWithin = local_inferSessWithin(sessInfo, mouseId)
S = height(sessInfo);

if ismember('sessWithin', sessInfo.Properties.VariableNames) && all(isfinite(sessInfo.sessWithin))
    sessWithin = sessInfo.sessWithin;
    return;
end

sessWithin = nan(S,1);
mouseU = unique(mouseId, 'stable');
for im = 1:numel(mouseU)
    idx = find(mouseId == mouseU(im));
    if ismember('header', sessInfo.Properties.VariableNames)
        [~,ord] = sort(string(sessInfo.header(idx)));
        idx = idx(ord);
    end
    sessWithin(idx) = (1:numel(idx))';
end
end

function R = local_cholStable(Sigma)
Sigma = (Sigma + Sigma')/2;
reg = 0;
for it = 1:10
    try
        R = chol(Sigma + reg*eye(size(Sigma,1)));
        return;
    catch
        if reg == 0
            reg = 1e-12;
        else
            reg = reg * 10;
        end
    end
end
error('Cholesky failed even after regularization. Sigma may be badly conditioned.');
end
