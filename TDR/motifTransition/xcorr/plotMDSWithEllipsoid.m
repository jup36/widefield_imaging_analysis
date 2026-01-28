function h = plotMDSWithEllipsoid(Y, sessInfo, cmap, ell, varargin)
%PLOTMDSWITHELLIPSOID  Plot MDS embedding using plot3mds styling + ellipsoid overlay.
%
%   h = plotMDSWithEllipsoid(Y, sessInfo, cmap, ell, 'Name', value, ...)
%
% Subspace plotting
%   'dims' : [] (default) or integer vector selecting columns of Y to plot.
%            Examples:
%              []      -> use all columns of Y (original behavior)
%              [1 2]   -> plot in MDS dimensions 1 and 2
%              [2]     -> plot in 1D using MDS dimension 2 (x-axis only)
%              [1 3 4] -> plot in 3D using dims 1,3,4 (if available)
%
% Figure saving
%   'figSaveDir' : base directory for saving figures (default below)
%   'printFig'   : true/false, whether to save figure (default false)
%   'trialTag'   : "" (default) OR one of {"bothTrials","goTrials","nogoTrials"}.
%                  If empty, tag is inferred from input variable name of Y via inputname(1).
%                  (inputname is fragile; passing trialTag is recommended for reliability.)
%
% View control (important for reproducible saving + filename suffix)
%   'view' : [] (default) or [az el]. If provided, sets axes view before saving.
%
% Saving rule:
%   saveBase = "xcorr_mds_acrossSessionTrj_"
%   suffix   = trialTag (or inferred from inputname(1) if trialTag is empty)
%   final    = saveBase + suffix + "_" + datestr(now,'mmddyy') + "_viewAZ_EL.pdf"
%
% NOTES
%   • Ellipsoid overlay/highlighting is only applied when numel(dims)==3
%     and ell has compatible mu/Sigma (3D in the selected subspace).
%   • Any unmatched Name/Value pairs are forwarded to plot3mds.

% -------------------- parse options --------------------
p = inputParser;
p.KeepUnmatched = true; % forward anything unknown to plot3mds

% subspace dims
p.addParameter('dims', [], @(v) isempty(v) || (isnumeric(v) && isvector(v) && all(v==round(v)) && all(v>=1)));

% view control (NEW)
p.addParameter('view', [], @(v) isempty(v) || (isnumeric(v) && numel(v)==2));

% ellipsoid overlay options
p.addParameter('ellAlpha', 0.12, @(v)isnumeric(v)&&isscalar(v)&&v>=0&&v<=1);
p.addParameter('ellEdgeAlpha', 0.05, @(v)isnumeric(v)&&isscalar(v)&&v>=0&&v<=1);
p.addParameter('ellN', 40, @(v)isnumeric(v)&&isscalar(v)&&v>=8);
p.addParameter('ellLineWidth', 0.5, @(v)isnumeric(v)&&isscalar(v)&&v>=0);
p.addParameter('ellColor', [0 0 0], @(v)isnumeric(v)&&numel(v)==3);
p.addParameter('highlightInside', true, @(v)islogical(v)||ismember(v,[0 1]));
p.addParameter('insideEdgeColor', 'k'); % can be 'k' or [0 0 0], etc.

% figure saving
p.addParameter('figSaveDir', compatiblepath("Z:\Rodent Data\dualImaging_parkj\collectFigure\motifTDR\xcorr_mds"), ...
    @(s) ischar(s) || isstring(s));
p.addParameter('printFig', false, @(x)islogical(x)&&isscalar(x));

% NEW: robust tag for saving (fixes inputname fragility)
p.addParameter('trialTag', "", @(s)ischar(s) || isstring(s));

p.parse(varargin{:});
opt = p.Results;

% -------------------- select subspace --------------------
dims = opt.dims;

if isempty(dims)
    Ysub = Y;
else
    dims = dims(:)'; % row
    if any(dims > size(Y,2))
        error('Requested dims exceed available columns of Y (size(Y,2)=%d).', size(Y,2));
    end
    Ysub = Y(:, dims);
end

% We'll need dim after subspace selection
[~, dim] = size(Ysub);

% -------------------- base plot: EXACT styling from plot3mds --------------------
nv = local_namedargs2cell(p.Unmatched);

% enforce no marker edge unless user supplied it
if ~isfield(p.Unmatched, 'markerEdgeColor')
    nv = [nv, {'markerEdgeColor','none'}]; %#ok<AGROW>
end

% Call plot3mds, but with Ysub instead of Y
h.base = plot3mds(Ysub, sessInfo, cmap, nv{:});
ax = h.base.ax;
hold(ax, 'on');

% Apply requested view early (so it also affects filename suffix on save)
if ~isempty(opt.view)
    try
        view(ax, opt.view(:).'); % [az el]
        drawnow;
    catch
    end
end

h.ell = struct('surf', [], 'muMarker', []);
h.meta = struct();
h.meta.dims = dims;
h.meta.Ysub = Ysub;

% -------------------- 1D/2D: no ellipsoid overlay; return --------------------
if dim ~= 3
    % Optional: label axes by chosen dims (helps interpretation)
    if ~isempty(dims)
        try
            if dim == 1
                xlabel(ax, sprintf('MDS%d', dims(1)));
            elseif dim == 2
                xlabel(ax, sprintf('MDS%d', dims(1)));
                ylabel(ax, sprintf('MDS%d', dims(2)));
            end
        catch
        end
    end
    hold(ax, 'off');

    % save (also for 1D/2D)
    h.save = localMaybeSaveFig(opt, ax);
    return;
end

% -------------------- sanity checks for ellipsoid --------------------
if isempty(ell) || ~isfield(ell,'mu') || ~isfield(ell,'Sigma')
    warning('plotMDSWithEllipsoid:MissingEll', ...
        'ell must contain fields mu and Sigma for 3D ellipsoid overlay. Skipping ellipsoid.');
    hold(ax, 'off');

    % save even if ellipsoid skipped
    h.save = localMaybeSaveFig(opt, ax);
    return;
end

% Handle ell.mu / ell.Sigma dimensioning:
muFull  = ell.mu(:)';
SigFull = ell.Sigma;

if ~isempty(dims) && numel(dims) == 3
    if numel(muFull) >= max(dims) && size(SigFull,1) >= max(dims)
        mu    = muFull(dims);
        Sigma = SigFull(dims, dims);
    else
        % fall back: assume provided ell already matches plotted space
        mu    = muFull(1:3);
        Sigma = SigFull(1:3,1:3);
    end
else
    mu    = muFull(1:3);
    Sigma = SigFull(1:3,1:3);
end

% confidence scaling
confLevel = 0.95;
if isfield(ell,'confLevel') && ~isempty(ell.confLevel)
    confLevel = ell.confLevel;
end
thr   = chi2inv(confLevel, 3);
scale = sqrt(thr);

% regularize covariance defensively
Sigma = (Sigma + Sigma')/2;
Sigma = Sigma + 1e-8*eye(3);

% eigendecomposition -> axes directions + radii
[V, L] = eig(Sigma);
radii  = scale * sqrt(max(diag(L), 0));

% sphere mesh
[xs, ys, zs] = sphere(opt.ellN);
U = [xs(:) ys(:) zs(:)]';

% transform sphere -> ellipsoid
E  = V * (diag(radii) * U);
Ex = reshape(E(1,:) + mu(1), size(xs));
Ey = reshape(E(2,:) + mu(2), size(ys));
Ez = reshape(E(3,:) + mu(3), size(zs));

% draw ellipsoid surface
h.ell.surf = surf(ax, Ex, Ey, Ez, ...
    'FaceAlpha', opt.ellAlpha, ...
    'EdgeAlpha', opt.ellEdgeAlpha, ...
    'LineWidth', opt.ellLineWidth, ...
    'FaceColor', opt.ellColor, ...
    'EdgeColor', opt.ellColor);

% mu marker
h.ell.muMarker = scatter3(ax, mu(1), mu(2), mu(3), 110, ...
    'Marker', 'p', ...
    'MarkerFaceColor', opt.ellColor, ...
    'MarkerEdgeColor', opt.ellColor, ...
    'MarkerFaceAlpha', 1.0, ...
    'MarkerEdgeAlpha', 1.0);

% keep ellipsoid out of legend
try
    set(get(get(h.ell.surf,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');
    set(get(get(h.ell.muMarker,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');
catch
end

% push ellipsoid behind points
try
    uistack(h.ell.surf, 'bottom');
catch
end

% -------------------- highlight points INSIDE ellipsoid --------------------
if opt.highlightInside
    inside = local_pointsInsideEllipsoid(Ysub, mu, Sigma, thr);

    mouseId    = local_inferMouseId(sessInfo);
    sessWithin = local_inferSessWithin(sessInfo, mouseId);

    mouseU = unique(mouseId, 'stable');

    for im = 1:numel(mouseU)
        idx = find(mouseId == mouseU(im));
        [~, ord] = sort(sessWithin(idx), 'ascend');
        idx = idx(ord);

        sc = h.base.scatters{im}; % handles ordered by sessWithin
        if isempty(sc) || ~all(isgraphics(sc)), continue; end

        for s = 1:numel(idx)
            if inside(idx(s))
                set(sc(s), ...
                    'MarkerEdgeColor', opt.insideEdgeColor, ...
                    'LineWidth', 2);
                if isprop(sc(s), 'MarkerEdgeAlpha')
                    set(sc(s), 'MarkerEdgeAlpha', 1.0);
                end
            else
                set(sc(s), ...
                    'MarkerEdgeColor', 'none', ...
                    'LineWidth', 0.5);
            end
        end
    end
end

% label axes by dims if provided
if ~isempty(dims)
    xlabel(ax, sprintf('MDS%d', dims(1)));
    ylabel(ax, sprintf('MDS%d', dims(2)));
    zlabel(ax, sprintf('MDS%d', dims(3)));
end

hold(ax, 'off');

% save
h.save = localMaybeSaveFig(opt, ax);

end

% ========================= local helpers =========================
function nv = local_namedargs2cell(s)
fn = fieldnames(s);
nv = cell(1, 2*numel(fn));
for i = 1:numel(fn)
    nv{2*i-1} = fn{i};
    nv{2*i}   = s.(fn{i});
end
end

function inside = local_pointsInsideEllipsoid(Y, mu, Sigma, thr)
X = Y - mu;           % Sx3
R = chol(Sigma);      % Sigma = R'*R
Z = X / R;            % Sx3
md2 = sum(Z.^2, 2);
inside = (md2 <= thr);
end

function mouseId = local_inferMouseId(sessInfo)
S = height(sessInfo);

if ismember('mouseId', sessInfo.Properties.VariableNames) && ~all(cellfun(@isempty, cellstr(string(sessInfo.mouseId))))
    mouseId = string(sessInfo.mouseId);
    return;
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

% -------- saving logic (robust) --------
function saveStruct = localMaybeSaveFig(opt, ax)
saveStruct = struct('didSave', false, 'saveDir', '', 'saveFile', '', 'tag', '');

if ~opt.printFig
    return;
end

figSaveDir = char(string(opt.figSaveDir));
if ~exist(figSaveDir, 'dir')
    mkdir(figSaveDir);
end

% Robust tag: prefer opt.trialTag, else fall back to inputname-based inference
tag = string(opt.trialTag);
if strlength(strtrim(tag)) == 0
    yName = inputname(1); % may be empty/fragile
    tag = localMapYNameToTag(yName);
else
    tag = strtrim(tag);
end

dstr = datestr(now, 'mmddyy');
base = 'xcorr_mds_acrossSessionTrj_';

% view suffix from actual axes view at save time
[az, el] = view(ax);
viewStr = sprintf('_view%d_%d', round(az), round(el));

fnBase  = sprintf('%s%s_%s%s', base, char(tag), dstr, viewStr);
outFile = fullfile(figSaveDir, [fnBase '.pdf']);

fig = ancestor(ax, 'figure');
set(fig, 'InvertHardcopy', 'off');
print(fig, outFile, '-dpdf', '-painters', '-bestfit');

saveStruct.didSave  = true;
saveStruct.saveDir  = figSaveDir;
saveStruct.saveFile = outFile;
saveStruct.tag      = char(tag);
end

function tag = localMapYNameToTag(yName)
% Map variable name of Y to save tag (best-effort fallback)
yName = string(yName);

switch yName
    case "Y_full"
        tag = "bothTrials";
    case "Y_go"
        tag = "goTrials";
    case "Y_nogo"
        tag = "nogoTrials";
    otherwise
        tag = "unknownTrials";
end
end
