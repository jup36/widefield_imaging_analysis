function [S, order, hDendro, hSim] = beta_cosine_map(beta, varargin)
% BETA_COSINE_MAP  Cosine-similarity of motif betas with separate dendrogram & heatmap.
%
% INPUT
%   beta   : [P x K] matrix (P predictors, K motifs) — columns are motifs
%
% NAME-VALUE (optional)
%   'Names'     : 1xK cellstr of motif labels (default {'M1','M2',...})
%   'Linkage'   : linkage method for clustering (default 'average')
%   'Colormap'  : colormap handle or name (default parula)
%   'CLim'      : [min max] for similarity image (default [-1 1])
%   'Title'     : base title string (default 'β cosine similarity')
%
% OUTPUT
%   S        : [K x K] cosine-similarity matrix
%   order    : 1xK index of motifs after hierarchical clustering
%   hDendro  : figure handle for dendrogram
%   hSim     : figure handle for similarity heatmap

% ---- options
p = inputParser;
p.addParameter('Names',    [], @(x) isempty(x) || (iscellstr(x) && isvector(x)));
p.addParameter('Linkage',  'average', @(s) ischar(s) || isstring(s));
p.addParameter('Colormap', parula,   @(c) true);
p.addParameter('CLim',     [-1 1],   @(v) isnumeric(v) && numel(v)==2);
p.addParameter('Title',    'β cosine similarity', @(s) ischar(s) || isstring(s));
p.parse(varargin{:});
opt = p.Results;

% ---- normalize columns (L2) to get cosine via dot products
bnorm = sqrt(sum(beta.^2, 1));
bnorm(bnorm==0) = eps;
Bunit = beta ./ bnorm;                  % [P x K]
S = Bunit' * Bunit;                     % [K x K], cosine similarity

K = size(beta,2);
if isempty(opt.Names)
    names = arrayfun(@(k) sprintf('M%02d', k), 1:K, 'uni', 0);
else
    names = opt.Names(:)';
end

% ---- build a clean distance matrix from cosine similarity ----
% Clamp S to [-1,1] to avoid tiny numerical overshoots
S = max(min(S, 1), -1);

D = 1 - S;                        % similarity -> distance (0 on diag ideally)
D = (D + D.')./2;                 % enforce exact symmetry
K = size(D,1);
D(1:K+1:end) = 0;                 % enforce exact zero diagonal
D(~isfinite(D)) = 0;              % guard against NaN/Inf from weird inputs
D = max(D, 0);                    % no negative distances
% (optional) also cap at 2 (since 1 - (-1) = 2)
D = min(D, 2);

% Handle tiny cases
if K < 2
    error('Need at least 2 motifs to cluster.');
elseif K == 2
    Z = linkage([D(1,2)], opt.Linkage);
else
    % squareform expects a **square**, **hollow**, **symmetric** distance matrix
    Z = linkage(squareform(D), opt.Linkage);
end

order = optimalleaforder(Z, D);         % nicer leaf order

% ---- figure 1: dendrogram
hDendro = figure('Name','β Cosine: Dendrogram','Color','w');
[H, T, perm] = dendrogram(Z, 0, 'Reorder', order);
set(gca, 'TickDir','out', 'Box','off');
if ~isempty(names)
    set(gca, 'XTick', 1:K, 'XTickLabel', names(perm), 'XTickLabelRotation', 45);
end
title(sprintf('%s — dendrogram (%s)', opt.Title, opt.Linkage));

% ---- figure 2: similarity heatmap (reordered)
hSim = figure('Name','β Cosine: Similarity Map','Color','w');
imagesc(S(order, order));
axis image; set(gca, 'TickDir','out', 'Box','off');
colormap(opt.Colormap); caxis(opt.CLim); colorbar;
if ~isempty(names)
    set(gca, 'XTick', 1:K, 'XTickLabel', names(order), 'XTickLabelRotation', 45);
    set(gca, 'YTick', 1:K, 'YTickLabel', names(order));
end
title(sprintf('%s — reordered by clustering', opt.Title));

end
