function hFig = plotHMMEmissionMeans(model, varargin)
% plotHMMEmissionMeans
% Visualize HMM Gaussian emission means (mu) using imagesc.
%
% Displays model.mu' (motifs x states)
%
% Terminology:
%   - model.mu : emission means (state-conditional motif centroids)
%   - rows     : motifs
%   - columns  : latent states
%
% Usage:
%   hFig = plotHMMEmissionMeans(model)
%   hFig = plotHMMEmissionMeans(model, 'MotifLabels', ...)
%
% Required input:
%   model.mu : S x M matrix (state x motif)
%   model.S  : number of states
%
% Optional name-value pairs:
%   'MotifLabels' : cell array of length M (default = {'M1','M2',...})
%   'FontSize'    : base font size (default = 12)
%   'Title'       : figure title
%   'Colormap'    : colormap (default = parula)
%   'ShowColorbar': true/false (default = true)
%
% Author: Junchol Park (helper: ChatGPT)

% -------------------- parse inputs --------------------
p = inputParser;
p.addRequired('model', @isstruct);
p.addParameter('MotifLabels', {}, @iscell);
p.addParameter('FontSize', 12, @isscalar);
p.addParameter('Title', 'HMM emission means (\mu)', @(x)ischar(x)||isstring(x));
p.addParameter('Colormap', parula, @(x)isnumeric(x));
p.addParameter('ShowColorbar', true, @islogical);
p.parse(model, varargin{:});

fs       = p.Results.FontSize;
titleStr = p.Results.Title;
cmap     = p.Results.Colormap;
showCB   = p.Results.ShowColorbar;

% -------------------- extract params --------------------
mu = model.mu;          % S x M
[S, M] = size(mu);

muT = mu';              % M x S (motifs x states)

% default motif labels
if isempty(p.Results.MotifLabels)
    motifLabels = arrayfun(@(m) sprintf('M%d', m), 1:M, 'UniformOutput', false);
else
    motifLabels = p.Results.MotifLabels;
end

% -------------------- plot --------------------
hFig = figure('Color','w','Name','HMM emission means');
imagesc(muT);
axis tight;

colormap(cmap);
if showCB
    colorbar;
end

xlabel('State', 'FontSize', fs+2, 'FontWeight','bold', 'Interpreter', 'none');
ylabel('Motif', 'FontSize', fs+2, 'FontWeight','bold', 'Interpreter', 'none');

set(gca, ...
    'XTick', 1:S, ...
    'XTickLabel', arrayfun(@(s) sprintf('S%d', s), 1:S, 'UniformOutput', false), ...
    'YTick', 1:M, ...
    'YTickLabel', motifLabels, ...
    'FontSize', fs, ...
    'FontWeight','bold');

title(titleStr, 'FontSize', fs+4, 'FontWeight','bold', 'Interpreter', 'none');

end
