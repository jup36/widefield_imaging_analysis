function h = plotMDSWithBehaviorFaces(Y, sessInfo, cmap, ell, behTbl, varargin)
% plotMDSWithBehaviorFaces
%   Plot MDS embedding (plot3mds style) + ellipsoid overlay, but encode a behavioral
%   metric in scatter FACE color (finite metric only). NaN metric sessions are OPEN circles.
%   Optionally scale marker size by the metric.
%
% REQUIRED
%   Y        : [S x D] MDS coordinates
%   sessInfo : table with S rows; must contain 'header' (or something parsable)
%   cmap     : [Nmice x 3] mouse colormap (same ordering used by plot3mds)
%   ell      : struct with fields mu, Sigma (optional for non-3D)
%   behTbl   : table with variable 'header' + metric columns (e.g., biasC)
%
% NAME–VALUE (this function)
%   'dims'      : [] or vector of dims to plot (default [])
%   'view'      : [] or [az el]
%   'metric'    : metric name or alias (default 'biasC'; alias 'bias' -> 'biasC')
%   'cLim'      : [] or [lo hi] for metric color scaling (default symmetric auto)
%   'colormap'  : colormap for metric coloring (default 'parula')
%   'flipCmap'  : true/false (default true). If true, flips metric colormap so higher metric -> "low end" colors.
%   'showColorbar' : true/false (default true)
%   'cbTickStep'    : colorbar tick spacing (default 0.5)
%
%   'baseSize'  : baseline marker size (default 28)
%   'sizeScale' : scalar controlling magnitude of size modulation (default 120)
%   'sizeBy'    : 'abs' (default) | 'raw' | 'none'
%   'openIfNaN' : true/false open circles for NaN metric (default true)
%   'openLineWidth' : default 1.5
%
% Ellipsoid opts
%   'ellAlpha','ellEdgeAlpha','ellN','ellLineWidth','ellColor'
%
% Saving (same as plotMDSWithEllipsoid pattern)
%   'printFig','figSaveDir','trialTag'
%
% NOTE
%   We DO NOT forward unknown NV pairs to plot3mds.

% -------------------- parse --------------------
p = inputParser;
p.KeepUnmatched = false;   % IMPORTANT: don't forward unknown stuff to plot3mds

% subspace
p.addParameter('dims', [], @(v) isempty(v) || (isnumeric(v)&&isvector(v)));

% behavior metric
p.addParameter('metric', 'biasC', @(s)ischar(s)||isstring(s));
p.addParameter('cLim', [], @(x) isempty(x) || (isnumeric(x)&&numel(x)==2));
p.addParameter('colormap', 'parula');
p.addParameter('flipCmap', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('showColorbar', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('cbTickStep', 0.5, @(x)isnumeric(x)&&isscalar(x)&&x>0);

% size rules
p.addParameter('baseSize', 28, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('sizeScale', 120, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('sizeBy', 'abs', @(s) any(strcmpi(string(s), ["abs","raw","none"])));
p.addParameter('openIfNaN', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('openLineWidth', 1.5, @(x)isnumeric(x)&&isscalar(x)&&x>0);

% view
p.addParameter('view', [], @(v) isempty(v) || (isnumeric(v)&&numel(v)==2));

% ellipsoid
p.addParameter('ellAlpha', 0.12, @(v)isnumeric(v)&&isscalar(v)&&v>=0&&v<=1);
p.addParameter('ellEdgeAlpha', 0.05, @(v)isnumeric(v)&&isscalar(v)&&v>=0&&v<=1);
p.addParameter('ellN', 40, @(v)isnumeric(v)&&isscalar(v)&&v>=8);
p.addParameter('ellLineWidth', 0.5, @(v)isnumeric(v)&&isscalar(v)&&v>=0);
p.addParameter('ellColor', [0 0 0], @(v)isnumeric(v)&&numel(v)==3);

% saving
p.addParameter('figSaveDir', compatiblepath("Z:\Rodent Data\dualImaging_parkj\collectFigure\motifTDR\xcorr_mds"), @(s)ischar(s)||isstring(s));
p.addParameter('printFig', false, @(x)islogical(x)&&isscalar(x));
p.addParameter('trialTag', "", @(s)ischar(s)||isstring(s));

p.parse(varargin{:});
opt = p.Results;

% Preserve the ANIMAL colormap (do NOT overwrite this later)
cmapAnimal = cmap;

% -------------------- subspace --------------------
dims = opt.dims;
if isempty(dims)
    Ysub = Y;
else
    dims = opt.dims(:)'; %#ok<NASGU>
    Ysub = Y(:, dims);
end
dim = size(Ysub,2);

% -------------------- align behTbl -> sessInfo.header --------------------
assert(istable(behTbl) && ismember('header', behTbl.Properties.VariableNames), ...
    'behTbl must be a table with variable ''header''.');

hdrSess = string(sessInfo.header(:));
[tf, loc] = ismember(hdrSess, string(behTbl.header));

% ----- resolve metric name (aliases) -----
metricReq = lower(strtrim(string(opt.metric)));
varNamesL = lower(string(behTbl.Properties.VariableNames));

metricResolved = metricReq;
if any(metricReq == ["bias","criterion","c"])
    metricResolved = "biasc";
end

% d' aliases (optional)
if any(metricReq == ["d'","dprime","dprm"])
    cand = ["dprime","dprm","d_prm","d"];
    hit = cand(ismember(cand, varNamesL));
    if ~isempty(hit), metricResolved = hit(1); end
end

if ~ismember(metricResolved, varNamesL)
    error('Requested metric "%s" not found in behTbl. Available: %s', ...
        char(opt.metric), strjoin(string(behTbl.Properties.VariableNames), ", "));
end

metricName = string(behTbl.Properties.VariableNames(varNamesL == metricResolved));
metricName = metricName(1);

metricVal = nan(numel(hdrSess),1);
metricVal(tf) = behTbl{loc(tf), metricName};

% -------------------- base plot (plot3mds) --------------------
% IMPORTANT: we ONLY send parameters plot3mds understands.
h.base = plot3mds(Ysub, sessInfo, cmapAnimal, ...
    'markerEdgeColor','none');   % keep your prior style
ax = h.base.ax;
hold(ax,'on');

if ~isempty(opt.view)
    view(ax, opt.view(:).');
end

% -------------------- ellipsoid overlay (only in 3D) --------------------
h.ell = struct('surf', [], 'muMarker', []);
if dim == 3 && ~isempty(ell) && isfield(ell,'mu') && isfield(ell,'Sigma')
    muFull  = ell.mu(:)';  SigFull = ell.Sigma;

    if ~isempty(opt.dims) && numel(opt.dims)==3 && numel(muFull)>=max(opt.dims) && size(SigFull,1)>=max(opt.dims)
        mu = muFull(opt.dims);
        Sigma = SigFull(opt.dims,opt.dims);
    else
        mu = muFull(1:3);
        Sigma = SigFull(1:3,1:3);
    end

    confLevel = 0.95;
    if isfield(ell,'confLevel') && ~isempty(ell.confLevel), confLevel = ell.confLevel; end
    thr = chi2inv(confLevel, 3);
    scale = sqrt(thr);

    Sigma = (Sigma+Sigma')/2 + 1e-8*eye(3);
    [V,L] = eig(Sigma);
    radii = scale * sqrt(max(diag(L),0));

    [xs,ys,zs] = sphere(opt.ellN);
    U = [xs(:) ys(:) zs(:)]';
    E  = V * (diag(radii) * U);

    Ex = reshape(E(1,:) + mu(1), size(xs));
    Ey = reshape(E(2,:) + mu(2), size(ys));
    Ez = reshape(E(3,:) + mu(3), size(zs));

    h.ell.surf = surf(ax, Ex, Ey, Ez, ...
        'FaceAlpha', opt.ellAlpha, ...
        'EdgeAlpha', opt.ellEdgeAlpha, ...
        'LineWidth', opt.ellLineWidth, ...
        'FaceColor', opt.ellColor, ...
        'EdgeColor', opt.ellColor);

    h.ell.muMarker = scatter3(ax, mu(1), mu(2), mu(3), 110, ...
        'Marker','p', 'MarkerFaceColor', opt.ellColor, 'MarkerEdgeColor', opt.ellColor);

    try
        set(get(get(h.ell.surf,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        set(get(get(h.ell.muMarker,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        uistack(h.ell.surf,'bottom');
    catch
    end
end

% -------------------- metric colormap + limits --------------------
% Build metric colormap WITHOUT overwriting animal cmap
cmapMetric = feval(opt.colormap, 256);
if opt.flipCmap
    cmapMetric = flipud(cmapMetric);
end
colormap(ax, cmapMetric);

if isempty(opt.cLim)
    mx = max(abs(metricVal), [], 'omitnan');
    if ~isfinite(mx) || mx==0, mx = 1; end
    clim = [-mx mx];
else
    clim = opt.cLim;
end
caxis(ax, clim);

% -------------------- size scaling --------------------
modeSize = lower(string(opt.sizeBy));
switch modeSize
    case "none"
        sz = opt.baseSize * ones(size(metricVal));
    case "raw"
        m0 = metricVal;
        mx = max(abs(m0), [], 'omitnan'); if ~isfinite(mx) || mx==0, mx=1; end
        sz = opt.baseSize + opt.sizeScale * (m0./mx);
    otherwise % "abs"
        m0 = abs(metricVal);
        mx = max(m0, [], 'omitnan'); if ~isfinite(mx) || mx==0, mx=1; end
        sz = opt.baseSize + opt.sizeScale * (m0./mx);
end
sz(~isfinite(sz)) = opt.baseSize;

% -------------------- recolor per-session scatters --------------------
mouseId    = local_inferMouseId(sessInfo);
sessWithin = local_inferSessWithin(sessInfo, mouseId);
mouseU = unique(mouseId, 'stable');

for im = 1:numel(mouseU)
    idx = find(mouseId == mouseU(im));
    [~, ord] = sort(sessWithin(idx), 'ascend');
    idx = idx(ord);

    sc = h.base.scatters{im};
    if isempty(sc) || ~all(isgraphics(sc)), continue; end

    for s = 1:numel(idx)
        ii = idx(s);
        m  = metricVal(ii);

        set(sc(s), 'SizeData', sz(ii));

        if ~isfinite(m) && opt.openIfNaN
            % OPEN circle: edge color = ANIMAL color
            set(sc(s), ...
                'MarkerFaceColor','none', ...
                'MarkerEdgeColor', cmapAnimal(im,:), ...
                'LineWidth', opt.openLineWidth);
        else
            % FILLED circle: face by metric, edge by ANIMAL color
            set(sc(s), ...
                'CData', m, ...
                'MarkerFaceColor','flat', ...
                'MarkerEdgeColor', cmapAnimal(im,:), ...
                'LineWidth', 0.75);
        end
    end
end

% -------------------- colorbar formatting --------------------
if opt.showColorbar
    h.cb = colorbar(ax);
    h.cb.Label.String = char(metricName);

    % Set ticks every cbTickStep (default 0.5)
    lo = clim(1); hi = clim(2);
    step = opt.cbTickStep;

    % ensure we include 0 if inside range (nice for bias)
    t0 = ceil(lo/step)*step : step : floor(hi/step)*step;
    if isempty(t0)
        t0 = [lo hi];
    end
    h.cb.Ticks = t0;
else
    h.cb = gobjects(0);
end

title(ax, sprintf('MDS colored by %s (NaN = open circles)', char(metricName)), 'Interpreter','none');
hold(ax,'off');

% -------------------- save --------------------
h.save = localMaybeSaveFig(opt, ax, metricName);

end

% ========================= local helpers =========================
function mouseId = local_inferMouseId(sessInfo)
S = height(sessInfo);

if ismember('mouseId', sessInfo.Properties.VariableNames)
    mouseId = string(sessInfo.mouseId);
    if all(strlength(mouseId)>0), return; end
end

if ismember('header', sessInfo.Properties.VariableNames)
    hdr = string(sessInfo.header);
    mouseId = strings(S,1);
    for ii = 1:S
        tok = regexp(hdr(ii), '(m\d{3,5})', 'tokens', 'once');
        if isempty(tok), mouseId(ii) = "mouse_" + string(ii);
        else,           mouseId(ii) = string(tok{1});
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

function saveStruct = localMaybeSaveFig(opt, ax, metricName)
saveStruct = struct('didSave', false, 'saveDir', '', 'saveFile', '', 'tag', '');

if ~opt.printFig, return; end

outDir = char(string(opt.figSaveDir));
if ~exist(outDir,'dir'), mkdir(outDir); end

tag = strtrim(string(opt.trialTag));
if strlength(tag)==0, tag = "unknownTrials"; end

[az, el] = view(ax);
dstr = datestr(now,'mmddyy');
fnBase = sprintf('xcorr_mds_behaviorFaces_%s_%s_%s_view%d_%d', char(tag), char(metricName), dstr, round(az), round(el));
outFile = fullfile(outDir, [fnBase '.pdf']);

fig = ancestor(ax,'figure');
set(fig,'InvertHardcopy','off');
print(fig, outFile, '-dpdf','-painters','-bestfit');

saveStruct.didSave  = true;
saveStruct.saveDir  = outDir;
saveStruct.saveFile = outFile;
saveStruct.tag      = char(tag);
end
