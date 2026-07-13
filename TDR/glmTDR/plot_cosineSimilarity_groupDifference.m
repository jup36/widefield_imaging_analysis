function h = plot_cosineSimilarity_groupDifference(csMean_byGroup, leafOrder, varargin)
%PLOT_COSINESIMILARITY_GROUPDIFFERENCE  Heatmap of (groupA - groupB) mean cosine similarity.
%
% SYNOPSIS
%   h = plot_cosineSimilarity_groupDifference(csMean_byGroup, leafOrder, ...)
%
% DESCRIPTION
%   Plots groupA's mean cosine similarity matrix minus groupB's, using a
%   diverging (blue-white-red) colormap centered at zero, reordered by
%   the shared leaf order (e.g. from the grand-mean clustering) so it's
%   directly comparable to the other heatmaps in this pipeline. This is
%   the direct visual complement to computeMotifClusterContrast_byGroup.m
%   -- it shows WHERE across the whole matrix groupA and groupB differ,
%   not just within a single pre-chosen cluster pair.
%
% INPUTS
%   csMean_byGroup : scalar struct with at least two fields, each holding
%                    a [K x K] mean cosine similarity matrix in NATURAL
%                    (unreordered) motif order (e.g. the csMean_byGroup
%                    output of analyzeMotifBetaCosineSimilarity_acrossAnimals.m).
%   leafOrder      : [K x 1] or [1 x K] permutation vector giving the
%                    shared display order (e.g. h.leafOrder from
%                    analyzeMotifBetaCosineSimilarity_acrossAnimals.m).
%
% NAME-VALUE ARGS
%   'groupA', 'groupB' : field names in csMean_byGroup to subtract as
%                        (groupA - groupB). Default: the first two
%                        fieldnames of csMean_byGroup, in that order.
%   'clim'              : [min max] color limits. Default: [] (symmetric,
%                        auto-scaled to the max absolute difference).
%   'figureScaleFactor' : figure size multiplier. Default: 1.2.
%   'visible'           : 'on'/'off'. Default: 'on'.
%   'figSaveDir', 'figSaveKeyword' : save-to-PDF options, same convention
%                        as other functions in this pipeline.
%
% OUTPUT
%   h : struct with fig/ax handles, the computed (reordered) difference
%       matrix, and opt.
%
% EXAMPLE
%   h = plot_cosineSimilarity_groupDifference(csMean_byGroup, h_cs.leafOrder, ...
%           'groupA', 'fast', 'groupB', 'slow', ...
%           'figSaveDir', figSaveDir, 'figSaveKeyword', 'clusterContrast');
%
% See also: analyzeMotifBetaCosineSimilarity_acrossAnimals, computeMotifClusterContrast_byGroup

grpNamesAll = fieldnames(csMean_byGroup);
assert(numel(grpNamesAll) >= 2, 'csMean_byGroup must have at least 2 group fields.');

p = inputParser;
p.addParameter('groupA', grpNamesAll{1}, @(s) ischar(s) || isstring(s));
p.addParameter('groupB', grpNamesAll{2}, @(s) ischar(s) || isstring(s));
p.addParameter('clim', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
p.addParameter('figureScaleFactor', 1.2, @(x) isnumeric(x) && x>0);
p.addParameter('visible', 'on', @(s) any(strcmpi(s, {'on','off'})));
p.addParameter('figSaveDir', {}, @(x) isempty(x) || ischar(x) || isstring(x) || iscell(x));
p.addParameter('figSaveKeyword', '', @(s) ischar(s) || isstring(s));
p.parse(varargin{:});
opt = p.Results;
opt.groupA = char(opt.groupA);
opt.groupB = char(opt.groupB);

assert(isfield(csMean_byGroup, opt.groupA), 'csMean_byGroup has no field "%s".', opt.groupA);
assert(isfield(csMean_byGroup, opt.groupB), 'csMean_byGroup has no field "%s".', opt.groupB);

A = csMean_byGroup.(opt.groupA);
B = csMean_byGroup.(opt.groupB);
assert(~isempty(A) && ~isempty(B), 'Both "%s" and "%s" must have non-empty matrices.', opt.groupA, opt.groupB);
assert(isequal(size(A), size(B)), '"%s" and "%s" matrices must be the same size.', opt.groupA, opt.groupB);

K = size(A,1);
leafOrder = leafOrder(:)';
assert(numel(leafOrder)==K, 'leafOrder must have K=%d elements.', K);

% normalize figSaveDir
figSaveDir = opt.figSaveDir;
if iscell(figSaveDir)
    if isempty(figSaveDir), figSaveDir = ""; else, figSaveDir = string(figSaveDir{1}); end
else
    figSaveDir = string(figSaveDir);
end
figSaveKeyword = string(opt.figSaveKeyword);

% -------- compute difference, reordered --------
Diff = A(leafOrder, leafOrder) - B(leafOrder, leafOrder);

if isempty(opt.clim)
    m = max(abs(Diff(:)));
    if m == 0, m = 1e-6; end   % guard against a fully-zero difference
    climVals = [-m m];
else
    climVals = opt.clim;
end

% -------- plot --------
h = struct();
h.opt = opt;
h.Diff = Diff;
h.leafOrder = leafOrder;

h.fig = figure('Color','w', 'Visible', opt.visible);
set(h.fig, 'Units', 'normalized');
pos = get(h.fig, 'Position');
pos(3:4) = pos(3:4) * opt.figureScaleFactor;
set(h.fig, 'Position', pos);

h.ax = axes('Parent', h.fig);
imagesc(h.ax, Diff, climVals);
axis(h.ax, 'square');
colormap(h.ax, diverging_redblue_(256));
colorbar(h.ax);
set(h.ax, 'XTick', 1:K, 'XTickLabel', string(leafOrder), ...
          'YTick', 1:K, 'YTickLabel', string(leafOrder));
xtickangle(h.ax, 90);
xlabel(h.ax, 'Motif (reordered)');
ylabel(h.ax, 'Motif (reordered)');
title(h.ax, sprintf('\\beta cosine similarity difference (%s - %s)', opt.groupA, opt.groupB), 'Interpreter', 'tex');

% -------- save (optional) --------
if strlength(figSaveDir) > 0
    if ~isfolder(figSaveDir)
        mkdir(figSaveDir);
    end
    dateStr = char(datetime("today","Format","MMddyy"));
    parts = strings(0,1);
    parts(end+1,1) = "betaCosineSim_diff";
    if strlength(figSaveKeyword) > 0, parts(end+1,1) = figSaveKeyword; end
    parts(end+1,1) = dateStr;
    figSaveName = strjoin(parts, "_");
    print(h.fig, fullfile(figSaveDir, figSaveName), '-dpdf', '-painters', '-bestfit');
end

end % function


% ===== helper: simple blue-white-red diverging colormap (no toolbox dependency) =====
function cmap = diverging_redblue_(n)
if nargin < 1, n = 256; end
x = linspace(0, 1, n)';
cmap = zeros(n, 3);
loColor  = [0.05 0.05 0.80];   % blue
midColor = [1.00 1.00 1.00];   % white
hiColor  = [0.80 0.05 0.05];   % red
for k = 1:n
    t = x(k);
    if t <= 0.5
        f = t / 0.5;
        cmap(k,:) = (1-f)*loColor + f*midColor;
    else
        f = (t - 0.5) / 0.5;
        cmap(k,:) = (1-f)*midColor + f*hiColor;
    end
end
end
