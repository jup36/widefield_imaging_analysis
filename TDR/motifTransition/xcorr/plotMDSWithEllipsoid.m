function h = plotMDSWithEllipsoid(Y, sessInfo, cmap, ell, varargin)
%PLOTMDSWITHELLIPSOID  Plot MDS embedding using plot3mds styling + ellipsoid overlay.
%
%   h = plotMDSWithEllipsoid(Y, sessInfo, cmap, ell, 'Name', value, ...)
%
% SYNOPSIS
%   This is a thin wrapper around plot3mds that preserves all of its styling:
%     - per-mouse trajectory colors (cmap)
%     - per-session alpha ramps (early -> late)
%     - per-session dot-size ramps (early -> late)
%   Then it overlays a Gaussian ellipsoid (from ell.mu / ell.Sigma) and
%   highlights ONLY those session points that fall inside the ellipsoid by
%   turning ON a black marker edge. By default, no marker edges are shown.
%
% INPUTS
%   Y        : [S x dim] embedding coordinates (dim = 2 or 3).
%   sessInfo : table with S rows, same order as Y. Recommended variables:
%                - mouseId (string/cellstr) OR header (for inference)
%                - sessWithin (numeric) (for within-mouse ordering)
%   cmap     : [nMice x 3] RGB colormap, one row per mouse (as in plot3mds).
%   ell      : struct with (at minimum):
%                - mu    : 1x3 mean in MDS space
%                - Sigma : 3x3 covariance in MDS space
%              Optional:
%                - confLevel : scalar in (0,1), default 0.95
%
% NAME–VALUE OPTIONS
%   All name–value args supported by plot3mds are forwarded as-is.
%
%   Ellipsoid overlay options:
%     'ellAlpha'      : surface transparency (default 0.12)
%     'ellEdgeAlpha'  : edge transparency   (default 0.05)
%     'ellN'          : sphere mesh resolution (default 40)
%     'ellLineWidth'  : ellipsoid edge line width (default 0.5)
%     'ellColor'      : ellipsoid face/edge RGB (default [0 0 0])
%     'highlightInside' : true/false, outline points inside ellipsoid (default true)
%     'insideEdgeColor' : marker edge color for inside points (default 'k')
%
% OUTPUT
%   h.base : handles from plot3mds (fig, ax, scatters, etc.)
%   h.ell  : struct with fields .surf and .muMarker
%
% NOTES
%   • This function assumes plot3mds defaults to markerEdgeColor='none'
%     (or you pass 'markerEdgeColor','none' into plot3mds).
%   • Ellipsoid + inside-point highlighting are implemented for dim==3.
%     For dim~=3, the base plot is drawn and the function returns.

% -------------------- parse options --------------------
p = inputParser;
p.KeepUnmatched = true; % forward anything unknown to plot3mds

% ellipsoid overlay options
p.addParameter('ellAlpha', 0.12, @(v)isnumeric(v)&&isscalar(v)&&v>=0&&v<=1);
p.addParameter('ellEdgeAlpha', 0.05, @(v)isnumeric(v)&&isscalar(v)&&v>=0&&v<=1);
p.addParameter('ellN', 40, @(v)isnumeric(v)&&isscalar(v)&&v>=8);
p.addParameter('ellLineWidth', 0.5, @(v)isnumeric(v)&&isscalar(v)&&v>=0);
p.addParameter('ellColor', [0 0 0], @(v)isnumeric(v)&&numel(v)==3);
p.addParameter('highlightInside', true, @(v)islogical(v)||ismember(v,[0 1]));
p.addParameter('insideEdgeColor', 'k'); % can be 'k' or [0 0 0], etc.

p.parse(varargin{:});
opt = p.Results;

% -------------------- base plot: EXACT styling from plot3mds --------------------
% Forward unmatched name-value pairs to plot3mds
nv = local_namedargs2cell(p.Unmatched);

% IMPORTANT: default should be no marker edge; enforce here if user didn't supply it
if ~isfield(p.Unmatched, 'markerEdgeColor')
    nv = [nv, {'markerEdgeColor','none'}]; %#ok<AGROW>
end

h.base = plot3mds(Y, sessInfo, cmap, nv{:});
ax = h.base.ax;
hold(ax, 'on');

h.ell = struct('surf', [], 'muMarker', []);

% -------------------- if not 3D, return after base plot --------------------
[~, dim] = size(Y);
if dim ~= 3
    return;
end

% -------------------- sanity checks for ellipsoid --------------------
if isempty(ell) || ~isfield(ell,'mu') || ~isfield(ell,'Sigma')
    warning('plotMDSWithEllipsoid:MissingEll', ...
        'ell must contain fields mu and Sigma for 3D ellipsoid overlay. Skipping ellipsoid.');
    return;
end

mu = ell.mu(:)';           % 1x3
Sigma = ell.Sigma;         % 3x3

% confidence scaling (as-is; user said keep this)
confLevel = 0.95;
if isfield(ell,'confLevel') && ~isempty(ell.confLevel)
    confLevel = ell.confLevel;
end
thr = chi2inv(confLevel, 3);        % threshold for squared Mahalanobis dist
scale = sqrt(thr);                  % radius scaling factor

% regularize covariance defensively
Sigma = (Sigma + Sigma')/2;
Sigma = Sigma + 1e-8*eye(3);

% eigendecomposition -> axes directions + radii
[V, L] = eig(Sigma);
radii = scale * sqrt(max(diag(L), 0));

% sphere mesh
[xs, ys, zs] = sphere(opt.ellN);
U = [xs(:) ys(:) zs(:)]';

% transform sphere -> ellipsoid
E = V * (diag(radii) * U);
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

% mu marker (small, visible)
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

% push ellipsoid behind points (so it doesn't visually "mask" dots)
try
    uistack(h.ell.surf, 'bottom');
catch
end

% -------------------- highlight points INSIDE ellipsoid --------------------
% Only turn on a marker edge for points whose Mahalanobis distance is within the ellipsoid.
if opt.highlightInside
    inside = local_pointsInsideEllipsoid(Y, mu, Sigma, thr);

    % Map global session rows -> scatter handles in h.base.scatters{im}
    mouseId    = local_inferMouseId(sessInfo);
    sessWithin = local_inferSessWithin(sessInfo, mouseId);

    mouseU = unique(mouseId, 'stable');

    for im = 1:numel(mouseU)
        idx = find(mouseId == mouseU(im));
        [~, ord] = sort(sessWithin(idx), 'ascend');
        idx = idx(ord);

        sc = h.base.scatters{im}; % handles, ordered by sessWithin in plot3mds
        if isempty(sc) || ~all(isgraphics(sc)), continue; end

        for s = 1:numel(idx)
            if inside(idx(s))
                set(sc(s), ...
                    'MarkerEdgeColor', opt.insideEdgeColor, ...
                    'LineWidth', 2);   % <-- thicker black outline

                if isprop(sc(s), 'MarkerEdgeAlpha')
                    set(sc(s), 'MarkerEdgeAlpha', 1.0);
                end
            else
                set(sc(s), ...
                    'MarkerEdgeColor', 'none', ...
                    'LineWidth', 0.5); % reset to a minimal default
            end

        end
    end
end

hold(ax, 'off');
end

% ========================= local helpers =========================
function nv = local_namedargs2cell(s)
% Convert struct of name/value pairs to 1x(2N) cell array
fn = fieldnames(s);
nv = cell(1, 2*numel(fn));
for i = 1:numel(fn)
    nv{2*i-1} = fn{i};
    nv{2*i}   = s.(fn{i});
end
end

function inside = local_pointsInsideEllipsoid(Y, mu, Sigma, thr)
% Return logical vector inside (Sx1): Mahalanobis^2 <= thr
X = Y - mu;           % Sx3
R = chol(Sigma);      % Sigma = R'*R
Z = X / R;            % Sx3
md2 = sum(Z.^2, 2);   % Sx1
inside = (md2 <= thr);
end

function mouseId = local_inferMouseId(sessInfo)
% Return Sx1 string mouse IDs, matching plot3mds inference.
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
% Return Sx1 numeric session indices within each mouse (1..nSess per mouse).
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
