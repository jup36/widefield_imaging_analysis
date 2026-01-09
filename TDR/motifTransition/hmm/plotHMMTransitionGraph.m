function hOut = plotHMMTransitionGraph(A, varargin)
% plotHMMTransitionGraph
% MATLAB-native digraph visualization of HMM transition matrix A.
%
% NEW (readability upgrades):
%  - Edge labels: 3 decimals + font scaling
%  - Node size scaling
%  - Nonlinear line-width mapping to emphasize small probabilities
%
% Name-Value options (new + old)
%   'Layout'            : 'circle' (default)
%   'NodeLabels'        : {} -> {'S1','S2',...}
%   'ShowEdgeLabels'    : true/false (default true)
%   'EdgeLabelDecimals' : decimals for edge labels (default 3)
%   'EdgeLabelFontScale': multiply edge-label font size (default 2)
%   'LineWidthRange'    : [min max] (default [0.5 6])
%   'LineWidthMap'      : 'linear' | 'power' | 'log' (default 'power')
%   'GammaExp'          : exponent for 'power' mapping (<1 boosts small) (default 0.4)
%   'ArrowSize'         : default 12
%   'NodeSizeScale'     : multiply node MarkerSize (default 3)
%   'PerState'          : true/false (default false)
%   'HighlightColor'    : RGB (default [0 0 0])
%   'FadeColor'         : RGB (default [0 0 0])
%   'TitlePrefix'       : string (default '')
%
% Output
%   hOut : struct with .G .edgeTable .hFig .hPlot

% -------------------- parse inputs --------------------
p = inputParser;
p.addRequired('A', @(x) isnumeric(x) && ismatrix(x) && size(x,1)==size(x,2));

p.addParameter('Layout', 'circle', @(x) ischar(x) || isstring(x));
p.addParameter('NodeLabels', {}, @(x) iscellstr(x) || isstring(x));
p.addParameter('ShowEdgeLabels', true, @(x) islogical(x) && isscalar(x));
p.addParameter('EdgeLabelDecimals', 3, @(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('EdgeLabelFontScale', 2, @(x) isnumeric(x) && isscalar(x) && x>0);

p.addParameter('LineWidthRange', [0.5 6], @(x) isnumeric(x) && numel(x)==2);
p.addParameter('LineWidthMap', 'power', @(x) any(strcmpi(x, {'linear','power','log'})));
p.addParameter('GammaExp', 0.4, @(x) isnumeric(x) && isscalar(x) && x>0);

p.addParameter('ArrowSize', 12, @(x) isnumeric(x) && isscalar(x));
p.addParameter('NodeSizeScale', 3, @(x) isnumeric(x) && isscalar(x) && x>0);

p.addParameter('PerState', false, @(x) islogical(x) && isscalar(x));
p.addParameter('HighlightColor', [0 0 0], @(x) isnumeric(x) && numel(x)==3);
p.addParameter('FadeColor', [0 0 0], @(x) isnumeric(x) && numel(x)==3);
p.addParameter('TitlePrefix', '', @(x) ischar(x) || isstring(x));
p.parse(A, varargin{:});
opt = p.Results;

A = opt.A;
S = size(A,1);

% default node labels
if isempty(opt.NodeLabels)
    nodeLabels = arrayfun(@(i) sprintf('S%d', i), 1:S, 'UniformOutput', false);
else
    nodeLabels = cellstr(opt.NodeLabels);
    assert(numel(nodeLabels)==S, 'NodeLabels must have length S.');
end

% -------------------- build digraph --------------------
G = digraph(A);
w = G.Edges.Weight;

% -------------------- map weights -> line widths --------------------
lwMin = opt.LineWidthRange(1);
lwMax = opt.LineWidthRange(2);

if isempty(w) || all(w==0)
    lw = lwMin * ones(numedges(G),1);
else
    wMax = max(w);
    wn = w / wMax;

    switch lower(opt.LineWidthMap)
        case 'linear'
            f = wn;

        case 'power'
            % gammaExp < 1 expands low end
            f = wn .^ opt.GammaExp;

        case 'log'
            % log mapping; protects zero
            eps0 = 1e-12;
            f = log10(wn + eps0) - log10(eps0);
            f = f ./ max(f); % normalize to [0,1]

        otherwise
            f = wn;
    end

    lw = lwMin + (lwMax - lwMin) .* f;
end

% -------------------- edge labels (format) --------------------
if opt.ShowEdgeLabels
    fmt = sprintf('%%.%df', opt.EdgeLabelDecimals); % e.g. '%.3f'
    edgeLab = arrayfun(@(x) sprintf(fmt, x), w, 'UniformOutput', false);
else
    edgeLab = {};
end

% -------------------- output struct init --------------------
hOut = struct();
hOut.G = G;
hOut.edgeTable = G.Edges;

% -------------------- plot (all edges) --------------------
if ~opt.PerState
    hFig = figure('Color','w');
    h = plot(G, ...
        'Layout', char(opt.Layout), ...
        'NodeLabel', nodeLabels, ...
        'ArrowSize', opt.ArrowSize);

    % --- node styling ---
    baseNodeSize = h.MarkerSize;
    h.MarkerSize = baseNodeSize * opt.NodeSizeScale;

    % --- edge styling ---
    h.LineWidth = lw;
    h.EdgeColor = repmat(opt.HighlightColor, numedges(G), 1);

    if opt.ShowEdgeLabels
        h.EdgeLabel = edgeLab;

        % Increase edge label font size (MATLAB stores as text objects)
        ax = ancestor(h, 'axes');
        txt = findall(ax, 'Type', 'text');
        % Filter to edge labels by matching any label string (safe for small graphs)
        edgeSet = string(edgeLab);
        for ii = 1:numel(txt)
            if any(strcmp(string(txt(ii).String), edgeSet))
                txt(ii).FontSize = txt(ii).FontSize * opt.EdgeLabelFontScale;
                txt(ii).FontWeight = 'bold';
            end
        end
    end

    ttl = 'HMM transition graph (all edges)';
    if strlength(string(opt.TitlePrefix)) > 0
        ttl = sprintf('%s | %s', string(opt.TitlePrefix), ttl);
    end
    title(ttl, 'Interpreter','none');

    hOut.hFig = hFig;
    hOut.hPlot = h;
    return;
end

% -------------------- per-state figures (highlight outgoing) --------------------
hFig = gobjects(1,S);
hPlt = gobjects(1,S);

ends = G.Edges.EndNodes; % M×2
src  = ends(:,1);

for s = 1:S
    hFig(s) = figure('Color','w');
    hPlt(s) = plot(G, ...
        'Layout', char(opt.Layout), ...
        'NodeLabel', nodeLabels, ...
        'ArrowSize', opt.ArrowSize);

    baseNodeSize = hPlt(s).MarkerSize;
    hPlt(s).MarkerSize = baseNodeSize * opt.NodeSizeScale;

    isOut = (src == s);

    lw2 = lw;
    lw2(~isOut) = max(lwMin, lw2(~isOut) * 0.6);
    hPlt(s).LineWidth = lw2;

    edgeCol = repmat(opt.FadeColor, numedges(G), 1);
    edgeCol(isOut,:) = repmat(opt.HighlightColor, sum(isOut), 1);
    hPlt(s).EdgeColor = edgeCol;

    if opt.ShowEdgeLabels
        hPlt(s).EdgeLabel = edgeLab;

        ax = ancestor(hPlt(s), 'axes');
        txt = findall(ax, 'Type', 'text');
        edgeSet = string(edgeLab);
        for ii = 1:numel(txt)
            if any(strcmp(string(txt(ii).String), edgeSet))
                txt(ii).FontSize = txt(ii).FontSize * opt.EdgeLabelFontScale;
                txt(ii).FontWeight = 'bold';
            end
        end
    end

    ttl = sprintf('HMM transitions: outgoing from %s', nodeLabels{s});
    if strlength(string(opt.TitlePrefix)) > 0
        ttl = sprintf('%s | %s', string(opt.TitlePrefix), ttl);
    end
    title(ttl, 'Interpreter','none');
end

hOut.hFig = hFig;
hOut.hPlot = hPlt;
end
