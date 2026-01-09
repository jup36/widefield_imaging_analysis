function hFig = plotHMMTransitionMatrix(model, varargin)
% plotHMMTransitionMatrix
% Visualize HMM state transition matrix using imagesc.
%
% Usage:
%   hFig = plotHMMTransitionMatrix(model)
%   hFig = plotHMMTransitionMatrix(model, 'ShowValues', true)
%
% Required input:
%   model.A : S x S transition matrix
%   model.S : number of states
%
% Optional name-value pairs:
%   'ShowValues'   (true/false, default = true)
%   'ValueFormat' (sprintf format, default = '%.3f')
%   'FontSize'    (default = 12)
%   'Title'       (char or string, default = 'HMM transition matrix')
%
% Author: Junchol Park (helper: ChatGPT)

% -------------------- parse inputs --------------------
p = inputParser;
p.addRequired('model', @isstruct);
p.addParameter('ShowValues', true, @islogical);
p.addParameter('ValueFormat', '%.3f', @ischar);
p.addParameter('FontSize', 12, @isscalar);
p.addParameter('Title', 'HMM transition matrix', @(x)ischar(x)||isstring(x));
p.parse(model, varargin{:});

showValues = p.Results.ShowValues;
valFmt     = p.Results.ValueFormat;
fs         = p.Results.FontSize;
titleStr   = p.Results.Title;

% -------------------- extract model params --------------------
A = model.A;
S = model.S;

% -------------------- plot --------------------
hFig = figure('Color','w','Name','HMM transition matrix');
imagesc(A);
axis square;
colormap(parula);
colorbar;

xlabel('To state j',   'FontSize', fs+2, 'FontWeight','bold');
ylabel('From state i', 'FontSize', fs+2, 'FontWeight','bold');

set(gca, ...
    'XTick', 1:S, ...
    'YTick', 1:S, ...
    'XTickLabel', arrayfun(@(s) sprintf('S%d', s), 1:S, 'UniformOutput', false), ...
    'YTickLabel', arrayfun(@(s) sprintf('S%d', s), 1:S, 'UniformOutput', false), ...
    'FontSize', fs, ...
    'FontWeight','bold');

title(titleStr, 'FontSize', fs+4, 'FontWeight','bold');

% -------------------- annotate values --------------------
if showValues
    for i = 1:S
        for j = 1:S
            text(j, i, sprintf(valFmt, A(i,j)), ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', ...
                'FontSize', fs-1, ...
                'FontWeight','bold', ...
                'Color','w');
        end
    end
end

end
