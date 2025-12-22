function h = plot3mds(Y, sessInfo, cmap, varargin)
%PLOT3MDS  Plot session embeddings (2D/3D) with per-mouse trajectories.
%
% h = plot3mds(Y, sessInfo, cmap, 'Name', value, ...)
%
% INPUTS
%   Y        : [S x dim] embedding coordinates (dim=2 or 3).
%   sessInfo : table with at least:
%              - mouseIdx OR mouseId OR header (used to infer mouse identity)
%              - sessWithin (optional; otherwise inferred from header order)
%              - header (recommended; used for labels / fallback inference)
%   cmap     : [nMice x 3] RGB colormap provided by user, one row per mouse.
%
% NAME–VALUE OPTIONS
%   'ax'           : existing axes handle (default: new fig/axes)
%   'view'         : 1x2 view angles for 3D (default [35 20])
%   'lineWidth'    : trajectory line width (default 1.5)
%   'markerEdgeColor' : marker edge color (default 'k')
%   'minAlpha'     : alpha for earliest sessions (default 0.25)
%   'maxAlpha'     : alpha for latest sessions   (default 1.00)
%   'minDotSize'   : dot size for earliest session (default 36)
%   'dotScaleMax'  : scale factor for latest dot vs earliest (default 2.0)
%   'title'        : plot title (default 'Session state space')
%   'showLegend'   : true/false (default true)
%
% OUTPUT (handles)
%   h.fig, h.ax, h.lines, h.scatters, h.mouseIds

% -------------------- parse inputs --------------------
p = inputParser;
p.addParameter('ax', [], @(x) isempty(x) || ishghandle(x));
p.addParameter('view', [35 20], @(v) isnumeric(v) && numel(v)==2);
p.addParameter('lineWidth', 1.5, @(v) isnumeric(v) && isscalar(v));
p.addParameter('markerEdgeColor', 'k');
p.addParameter('minAlpha', 0.25, @(v) isnumeric(v) && isscalar(v) && v>=0 && v<=1);
p.addParameter('maxAlpha', 1.00, @(v) isnumeric(v) && isscalar(v) && v>=0 && v<=1);
p.addParameter('minDotSize', 36, @(v) isnumeric(v) && isscalar(v) && v>0);
p.addParameter('dotScaleMax', 2.0, @(v) isnumeric(v) && isscalar(v) && v>=1);
p.addParameter('title', 'Session state space (directional motif xcorr structure)', @(s) ischar(s) || isstring(s));
p.addParameter('showLegend', true, @(v) islogical(v) || ismember(v,[0 1]));
p.parse(varargin{:});
opt = p.Results;

[S, dim] = size(Y);
assert(dim==2 || dim==3, 'Y must be Sx2 or Sx3.');
assert(height(sessInfo)==S, 'sessInfo must have S rows (same as Y).');

% -------------------- infer mouseId per session --------------------
mouseId = inferMouseId(sessInfo);           % Sx1 string
mouseIdU = unique(mouseId, 'stable');       % preserve order of appearance
nMice = numel(mouseIdU);

if size(cmap,1) < nMice || size(cmap,2) ~= 3
    error('cmap must be [nMice x 3] with nMice >= number of unique mice (%d).', nMice);
end

% -------------------- infer session order within mouse --------------------
sessWithin = inferSessWithin(sessInfo, mouseId);  % Sx1 numeric, 1..nSess(mouse)

% -------------------- setup figure/axes --------------------
if isempty(opt.ax)
    h.fig = figure('Color','w');
    h.ax  = axes('Parent', h.fig);
else
    h.ax = opt.ax;
    h.fig = ancestor(h.ax, 'figure');
end
hold(h.ax, 'on');

h.lines = gobjects(nMice,1);
h.scatters = cell(nMice,1);
h.mouseIds = mouseIdU;

% -------------------- plot per mouse --------------------
for im = 1:nMice
    thisId = mouseIdU(im);
    idx = (mouseId == thisId);

    % order within mouse by sessWithin
    [~, ord] = sort(sessWithin(idx), 'ascend');
    Ym = Y(idx, :);
    Ym = Ym(ord, :);
    sw = sessWithin(idx);
    sw = sw(ord);

    nSess = size(Ym,1);
    if nSess == 0, continue; end

    % alpha and dot size ramps (earliest -> latest)
    if nSess == 1
        alphaVec = opt.maxAlpha;
        sizeVec  = opt.minDotSize * opt.dotScaleMax;
    else
        alphaVec = linspace(opt.minAlpha, opt.maxAlpha, nSess);
        sizeVec  = opt.minDotSize * linspace(1, opt.dotScaleMax, nSess);
    end

    baseCol = cmap(im,:);

    % --- trajectory as segmented line so alpha can vary per segment ---
    for s = 1:(nSess-1)
        colA = [baseCol alphaVec(s+1)];  % use later alpha for segment visibility
        if dim == 3
            plot3(h.ax, Ym(s:s+1,1), Ym(s:s+1,2), Ym(s:s+1,3), '-', ...
                'LineWidth', opt.lineWidth, 'Color', colA);
        else
            plot(h.ax, Ym(s:s+1,1), Ym(s:s+1,2), '-', ...
                'LineWidth', opt.lineWidth, 'Color', colA);
        end
    end

    % --- markers (session points) with per-point alpha and size ---
    sc = gobjects(nSess,1);
    for s = 1:nSess
        colA = [baseCol alphaVec(s)];
        if dim == 3
            sc(s) = scatter3(h.ax, Ym(s,1), Ym(s,2), Ym(s,3), sizeVec(s), ...
                'MarkerFaceColor', colA(1:3), ...
                'MarkerFaceAlpha', colA(4), ...
                'MarkerEdgeColor', opt.markerEdgeColor, ...
                'MarkerEdgeAlpha', colA(4));
        else
            sc(s) = scatter(h.ax, Ym(s,1), Ym(s,2), sizeVec(s), ...
                'MarkerFaceColor', colA(1:3), ...
                'MarkerFaceAlpha', colA(4), ...
                'MarkerEdgeColor', opt.markerEdgeColor);
        end
    end
    h.scatters{im} = sc;

    % dummy handle for legend (solid color)
    if dim == 3
        h.lines(im) = plot3(h.ax, nan, nan, nan, '-o', ...
            'Color', baseCol, 'MarkerFaceColor', baseCol, ...
            'MarkerEdgeColor', opt.markerEdgeColor, ...
            'LineWidth', opt.lineWidth, 'DisplayName', char(thisId));
    else
        h.lines(im) = plot(h.ax, nan, nan, '-o', ...
            'Color', baseCol, 'MarkerFaceColor', baseCol, ...
            'MarkerEdgeColor', opt.markerEdgeColor, ...
            'LineWidth', opt.lineWidth, 'DisplayName', char(thisId));
    end
end

% -------------------- cosmetics --------------------
xlabel(h.ax, 'MDS1'); ylabel(h.ax, 'MDS2');
if dim == 3, zlabel(h.ax, 'MDS3'); end
title(h.ax, opt.title, 'Interpreter','none');
grid(h.ax, 'on'); axis(h.ax, 'equal');
set(h.ax, 'TickDir','out', 'LineWidth', 1);

if dim == 3
    view(h.ax, opt.view(1), opt.view(2));
end

if opt.showLegend
    legend(h.ax, h.lines, 'Location', 'bestoutside', 'Interpreter','none');
end

end % plot3mds


% ========================= helpers =========================

function mouseId = inferMouseId(sessInfo)
% Return Sx1 string mouse IDs.
S = height(sessInfo);

if ismember('mouseId', sessInfo.Properties.VariableNames) && ~all(cellfun(@isempty, cellstr(string(sessInfo.mouseId))))
    mouseId = string(sessInfo.mouseId);
    return;
end

% Fallback: try parse from header
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

% Fallback: try parse from mouseIdx if present
if ismember('mouseIdx', sessInfo.Properties.VariableNames)
    mouseId = "mouse_" + string(sessInfo.mouseIdx);
    return;
end

% Ultimate fallback
mouseId = "mouse_" + string((1:S)');
end


function sessWithin = inferSessWithin(sessInfo, mouseId)
% Return Sx1 numeric session indices within each mouse (1..nSess per mouse).
S = height(sessInfo);

if ismember('sessWithin', sessInfo.Properties.VariableNames) && all(isfinite(sessInfo.sessWithin))
    sessWithin = sessInfo.sessWithin;
    return;
end

% If no sessWithin, infer by order of appearance within each mouse.
sessWithin = nan(S,1);
mouseU = unique(mouseId, 'stable');

for im = 1:numel(mouseU)
    idx = find(mouseId == mouseU(im));
    % if headers contain dates, try sorting by header (often chronological)
    if ismember('header', sessInfo.Properties.VariableNames)
        [~,ord] = sort(string(sessInfo.header(idx)));
        idx = idx(ord);
    end
    sessWithin(idx) = (1:numel(idx))';
end
end
