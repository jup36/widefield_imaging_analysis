function [fh, R, P, order] = plotDesignCorr(X, varargin)
% plotDesignCorr  Visualize column-wise correlations in a design matrix
%
%   [fh, R, P, order] = plotDesignCorr(X, 'Name', Value, ...)
%
% Inputs
%   X           : [N x P] design matrix (rows = samples, cols = predictors)
%
% Name-Value pairs (all optional)
%   'VarNames'  : 1xP cell array of variable names (default: {'x1','x2',...})
%   'Method'    : 'Pearson' (default), 'Spearman', or 'Kendall'
%   'Cluster'   : true/false (default: true) — hierarchical clustering of variables
%   'Mask'      : 'none' (default), 'upper', or 'lower' — show only part of matrix
%   'ShowValues': true/false (default: false) — overlay correlation values
%   'Title'     : char/string for figure title (default: 'Column Correlations')
%
% Outputs
%   fh    : figure handle
%   R     : [P x P] correlation matrix
%   P     : [P x P] p-value matrix from correlation test
%   order : permutation of columns after clustering (identity if Cluster=false)
%
% Example
%   X = randn(35100,20); X(:,3)=X(:,2)*0.9 + 0.1*randn(35100,1);
%   [fh,R,P,ord] = plotDesignCorr(X,'VarNames',compose('x%d',1:20), ...
%                                 'Cluster',true,'Mask','none','ShowValues',false);

% ---------- Parse inputs ----------
p = inputParser;
p.addRequired('X', @(z) isnumeric(z) && ndims(z)==2);
p.addParameter('VarNames', [], @(z) isempty(z) || (iscellstr(z) && numel(z)==size(X,2)));
p.addParameter('Method','Pearson', @(z) any(strcmpi(z,{'Pearson','Spearman','Kendall'})));
p.addParameter('Cluster', true, @(z) islogical(z) && isscalar(z));
p.addParameter('Mask','none', @(z) any(strcmpi(z,{'none','upper','lower'})));
p.addParameter('ShowValues', false, @(z) islogical(z) && isscalar(z));
p.addParameter('Title','Column Correlations', @(z) isstring(z) || ischar(z));
p.parse(X, varargin{:});

varNames   = p.Results.VarNames;
method     = validatestring(p.Results.Method, {'Pearson','Spearman','Kendall'});
doCluster  = p.Results.Cluster;
mask       = validatestring(p.Results.Mask, {'none','upper','lower'});
showVals   = p.Results.ShowValues;
ttl        = p.Results.Title;

[N,P] = size(X);
if isempty(varNames)
    varNames = compose('x%d', 1:P);
end

% ---------- Compute correlations ----------
% Pairwise handling of NaNs so you can pass partially missing columns
switch lower(method)
    case 'pearson'
        [R, Pval] = corrcoef(X, 'Rows','pairwise');
    otherwise
        % corr supports 'Type' for Spearman/Kendall with 'Rows','pairwise'
        [R, Pval] = corr(X, 'Type',method, 'Rows','pairwise');
end

% Guard against numerical issues
R = max(min(R,1),-1);

% ---------- Optional clustering of variables ----------
order = 1:P;
if doCluster && P>1
    % Distance on |R| so positively/negatively correlated variables cluster
    D = 1 - abs(R);
    % Ensure a proper distance vector (upper triangle)
    Dvec = squareform(triu(D,1));
    Z = linkage(Dvec,'average');
    order = optimalleaforder(Z, D);
    R = R(order, order);
    Pval = Pval(order, order);
    varNames = varNames(order);
end

% ---------- Apply mask if requested ----------
M = true(P);
switch lower(mask)
    case 'upper'
        M = tril(true(P));  % keep lower incl diag
    case 'lower'
        M = triu(true(P));  % keep upper incl diag
    case 'none'
        % keep all
end

Rmasked = R;
Rmasked(~M) = NaN;

% ---------- Plot ----------
fh = figure('Color','w'); 
ax = axes(fh); %#ok<LAXES>
imagesc(ax, Rmasked, [-1 1]); axis(ax,'image'); hold(ax,'on');
colormap(ax, divergingMap()); colorbar;
set(ax,'XTick',1:P,'XTickLabel',varNames, 'XTickLabelRotation',45);
set(ax,'YTick',1:P,'YTickLabel',varNames);
title(ax, ttl, 'Interpreter','none');
grid(ax,'on'); ax.GridColor = [0.8 0.8 0.8];

% Overlay value labels if requested
if showVals
    for i = 1:P
        for j = 1:P
            if ~isnan(Rmasked(i,j))
                text(j,i, sprintf('%.2f', R(i,j)), ...
                    'HorizontalAlignment','center','VerticalAlignment','middle', ...
                    'Color', pickTextColor(R(i,j)), 'FontSize',8, 'Parent', ax);
            end
        end
    end
end

% Light separators
xline(ax, 0.5:P+0.5, ':', 'Color',[0.85 0.85 0.85]);
yline(ax, 0.5:P+0.5, ':', 'Color',[0.85 0.85 0.85]);

% Outputs
if nargout > 1
    P = Pval; %#ok<NASGU>  % expose p-values as 3rd output
end

end % function

% ---------- Helpers ----------
function cmap = divergingMap()
% Simple red-white-blue diverging map centered at zero
n = 256;
r = [(0:1/(n/2-1):1)'; ones(n/2,1)];
g = [(0:1/(n/2-1):1)'; (1:-1/(n/2-1):0)'];
b = [ones(n/2,1); (1:-1/(n/2-1):0)'];
cmap = [r g b];
end

function c = pickTextColor(r)
% Dark text on light, light text on dark
if abs(r) > 0.5
    c = [1 1 1]*0.95; % white-ish
else
    c = [0 0 0];      % black
end
end
