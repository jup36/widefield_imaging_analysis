function out = plotSessionConvergenceToFastLearnerRef(Y, sessInfo, fastLearnerIdC, varargin)
%PLOTSESSIONCONVERGENCETOFASTLEARNERREF
%  Compute + visualize per-session convergence (distance) to a reference point
%  defined as the mean of the final N sessions of a set of fast learners.
%
%  out = plotSessionConvergenceToFastLearnerRef(Y, sessInfo, fastLearnerIdC, 'Name', value, ...)
%
% INPUTS
%  Y              : [S x dim] MDS coordinates (dim>=1). Rows correspond to sessInfo rows.
%  sessInfo       : table with S rows. Needs vars: mouseId, sessWithin (or header).
%  fastLearnerIdC : cellstr or string array of mouse IDs, e.g., {'m1044','m1045',...}
%
% NAME–VALUE PAIRS
%  'nFinalSess'      : number of final sessions per fast learner to define ref cloud (default 3)
%  'dims'            : dims used for distance, e.g. [1 2] or [1 2 3]. Default [] -> 1:min(3,size(Y,2))
%  'metric'          : 'euclidean' (default) or 'mahalanobis'
%  'Sigma'           : covariance for Mahalanobis. Default [] -> computed from ref points in chosen dims
%  'mu'              : reference mean. Default [] -> computed from ref points in chosen dims
%  'useEllipsoidThr' : true/false (default true). If true and metric='mahalanobis',
%                      mark points inside chi2 conf region.
%  'confLevel'       : chi2 conf level for in-cloud marking (default 0.95)
%  'cmap'            : [nMice x 3] colormap. Default [] -> lines(nMice)
%  'plotFastLearners': true/false (default true) highlight fast learners with thicker lines
%  'normalize'       : 'none' (default) | 'subtractFirst' | 'zscoreWithinMouse'
%  'xField'          : 'sessWithin' (default) or any numeric sessInfo field for x-axis
%  'titleStr'        : figure title (default '')
%
% OUTPUT (struct)
%  out.ref.idxPts    : indices used for ref
%  out.ref.mu        : 1 x d reference mean (in chosen dims)
%  out.ref.Sigma     : d x d covariance (in chosen dims)
%  out.ref.confLevel : conf level used (if applicable)
%  out.distTbl       : table with one row per session: mouseId, sessWithin, dist, insideRef
%  out.hFig, out.ax  : figure handles
%
% NOTE
%  If metric='euclidean', "insideRef" is always false unless you provide your own rule downstream.

% -------------------- parse inputs --------------------
p = inputParser;
p.addParameter('nFinalSess', 3, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('dims', [], @(x) isempty(x) || (isnumeric(x)&&isvector(x)&&all(x>=1)));
p.addParameter('metric', 'euclidean', @(s)ischar(s)||isstring(s));
p.addParameter('Sigma', [], @(x) isempty(x) || (isnumeric(x)&&ismatrix(x)));
p.addParameter('mu', [], @(x) isempty(x) || (isnumeric(x)&&isvector(x)));
p.addParameter('useEllipsoidThr', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('confLevel', 0.95, @(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);
p.addParameter('cmap', [], @(x) isempty(x) || (isnumeric(x)&&size(x,2)==3));
p.addParameter('plotFastLearners', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('normalize', 'none', @(s)ischar(s)||isstring(s));
p.addParameter('xField', 'sessWithin', @(s)ischar(s)||isstring(s));
p.addParameter('titleStr', '', @(s)ischar(s)||isstring(s));
p.parse(varargin{:});
opt = p.Results;

metric = lower(string(opt.metric));
normalizeMode = lower(string(opt.normalize));
xField = string(opt.xField);

% -------------------- infer required sessInfo fields --------------------
S = size(Y,1);
assert(height(sessInfo)==S, 'sessInfo must have same number of rows as Y.');

if ~ismember('mouseId', sessInfo.Properties.VariableNames)
    error('sessInfo must contain variable "mouseId".');
end
mouseId = string(sessInfo.mouseId);

if ismember(xField, sessInfo.Properties.VariableNames)
    xVals = sessInfo.(xField);
    if ~isnumeric(xVals)
        error('sessInfo.%s must be numeric for x-axis.', xField);
    end
else
    error('sessInfo does not contain xField "%s".', xField);
end

% If sessWithin missing but xField is sessWithin, try to infer
if xField=="sessWithin" && ~all(isfinite(xVals))
    % attempt fallback inference by sorting headers
    if ~ismember('header', sessInfo.Properties.VariableNames)
        error('sessWithin has NaNs and sessInfo has no "header" for inference.');
    end
    xVals = local_inferSessWithinFromHeader(sessInfo, mouseId);
end

% -------------------- choose dims --------------------
if isempty(opt.dims)
    dims = 1:min(3, size(Y,2));
else
    dims = opt.dims(:)';
    if any(dims > size(Y,2))
        error('dims exceeds available columns in Y.');
    end
end
Yd = Y(:, dims);

d = numel(dims);

% -------------------- select reference points (final N sessions of fast learners) --------------------
fastLearnerIdC = string(fastLearnerIdC(:));
idxRef = [];

for k = 1:numel(fastLearnerIdC)
    mid = fastLearnerIdC(k);
    idxM = find(mouseId == mid);
    if isempty(idxM), continue; end

    [~, ord] = sort(xVals(idxM), 'ascend');
    idxM = idxM(ord);

    nTake = min(opt.nFinalSess, numel(idxM));
    idxRef = [idxRef; idxM(end-nTake+1:end)]; %#ok<AGROW>
end
idxRef = unique(idxRef, 'stable');

if isempty(idxRef)
    error('No reference sessions found. Check fastLearnerIdC and sessInfo.mouseId.');
end

% -------------------- compute reference mu / Sigma (in chosen dims) --------------------
if isempty(opt.mu)
    mu = mean(Yd(idxRef,:), 1, 'omitnan');
else
    mu = opt.mu(:)';
    if numel(mu) ~= d
        error('Provided mu must have length %d to match dims.', d);
    end
end

if metric=="mahalanobis"
    if isempty(opt.Sigma)
        Sigma = cov(Yd(idxRef,:), 'omitrows');
    else
        Sigma = opt.Sigma;
    end
    Sigma = (Sigma + Sigma')/2;
    Sigma = Sigma + 1e-8*eye(d);
else
    Sigma = [];
end

% -------------------- compute distance per session --------------------
dist = nan(S,1);
insideRef = false(S,1);

switch metric
    case "euclidean"
        dif = Yd - mu;
        dist = sqrt(sum(dif.^2, 2, 'omitnan'));

    case "mahalanobis"
        % md^2 = (x-mu) * inv(Sigma) * (x-mu)'
        R = chol(Sigma);
        Z = (Yd - mu) / R;
        md2 = sum(Z.^2, 2, 'omitnan');
        dist = sqrt(md2);

        if opt.useEllipsoidThr
            thr = chi2inv(opt.confLevel, d);
            insideRef = (md2 <= thr);
        end

    otherwise
        error('Unknown metric "%s". Use "euclidean" or "mahalanobis".', metric);
end

% -------------------- normalization (optional) --------------------
mouseU = unique(mouseId, 'stable');
distPlot = dist;

switch normalizeMode
    case "none"
        % nothing

    case "subtractfirst"
        for im = 1:numel(mouseU)
            idx = find(mouseId == mouseU(im));
            [~, ord] = sort(xVals(idx), 'ascend');
            idx = idx(ord);
            if isempty(idx), continue; end
            distPlot(idx) = dist(idx) - dist(idx(1));
        end

    case "zscorewithinmouse"
        for im = 1:numel(mouseU)
            idx = find(mouseId == mouseU(im));
            v = dist(idx);
            muM = mean(v, 'omitnan');
            sdM = std(v, 'omitnan');
            if sdM==0 || ~isfinite(sdM), sdM = 1; end
            distPlot(idx) = (v - muM) ./ sdM;
        end

    otherwise
        error('Unknown normalize mode "%s".', normalizeMode);
end

% -------------------- colormap --------------------
if isempty(opt.cmap)
    cmap = lines(numel(mouseU));
else
    cmap = opt.cmap;
    if size(cmap,1) < numel(mouseU)
        error('Provided cmap has %d rows but need >= %d.', size(cmap,1), numel(mouseU));
    end
    cmap = cmap(1:numel(mouseU),:);
end

isFast = ismember(mouseU, fastLearnerIdC);

% -------------------- plot --------------------
hFig = figure('Color','w');
ax = axes('Parent', hFig); hold(ax,'on');
hLeg = gobjects(numel(mouseU),1);

for im = 1:numel(mouseU)
    idx = find(mouseId == mouseU(im));
    [xSorted, ord] = sort(xVals(idx), 'ascend');
    idx = idx(ord);
    ySorted = distPlot(idx);

    lw = 1.5;
    if opt.plotFastLearners && isFast(im)
        lw = 2.8;
    end

    % LINE (this will be used for legend)
    hLine = plot(ax, xSorted, ySorted, '-', 'Color', cmap(im,:), 'LineWidth', lw);
    hLeg(im) = hLine;

    % POINTS (hide from legend)
    hPts = plot(ax, xSorted, ySorted, 'o', ...
        'Color', cmap(im,:), 'MarkerFaceColor', cmap(im,:), 'MarkerSize', 5);
    hPts.HandleVisibility = 'off';

    % INSIDE points (also hide from legend)
    if any(insideRef(idx))
        idxIn = idx(insideRef(idx));
        hIn = plot(ax, xVals(idxIn), distPlot(idxIn), 'o', ...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', cmap(im,:), ...
            'LineWidth', 1.5, 'MarkerSize', 6);
        hIn.HandleVisibility = 'off';
    end
end

% --- legend (lines only) ---
legend(ax, hLeg, cellstr(mouseU), 'Location','eastoutside');

switch normalizeMode
    case "subtractfirst"
        ylab = sprintf('Distance to fast-learner ref (%s; dims=%s) [minus first]', metric, mat2str(dims));
    case "zscorewithinmouse"
        ylab = sprintf('Distance to fast-learner ref (%s; dims=%s) [z within mouse]', metric, mat2str(dims));
    otherwise
        ylab = sprintf('Distance to fast-learner ref (%s; dims=%s)', metric, mat2str(dims));
end
ylabel(ax, ylab);

if strlength(opt.titleStr) > 0
    title(ax, opt.titleStr, 'Interpreter','none');
end

grid(ax,'on');
box(ax,'off');

% legend (mouse IDs)
legStr = cellstr(mouseU);
legend(ax, legStr, 'Location','eastoutside');

% -------------------- pack outputs --------------------
distTbl = table((1:S)', mouseId, xVals, dist, distPlot, insideRef, ...
    'VariableNames', {'sessIdx','mouseId', char(xField), 'distRaw', 'distPlot', 'insideRef'});

out = struct();
out.ref = struct();
out.ref.idxPts    = idxRef;
out.ref.mu        = mu;
out.ref.Sigma     = Sigma;
out.ref.confLevel = opt.confLevel;
out.ref.dims      = dims;
out.ref.metric    = metric;

out.distTbl = distTbl;
out.hFig = hFig;
out.ax = ax;

end

% -------------------------------------------------------------------------
function sessWithin = local_inferSessWithinFromHeader(sessInfo, mouseId)
% fallback: infer sessWithin by sorting headers per mouse
S = height(sessInfo);
sessWithin = nan(S,1);

if ~ismember('header', sessInfo.Properties.VariableNames)
    error('Cannot infer sessWithin: sessInfo lacks "header".');
end

hdr = string(sessInfo.header);
mouseU = unique(mouseId, 'stable');

for im = 1:numel(mouseU)
    idx = find(mouseId == mouseU(im));
    [~, ord] = sort(hdr(idx), 'ascend');
    idx = idx(ord);
    sessWithin(idx) = (1:numel(idx))';
end
end
