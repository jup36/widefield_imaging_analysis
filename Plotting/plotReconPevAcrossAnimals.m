function h = plotReconPevAcrossAnimals(rezReconPevC, varargin)
% plotReconPevAcrossAnimals
%
% Scatter plot of mean reconstruction PEV across sessions, grouped by animal.
%
% INPUT
%   rezReconPevC : [nAnimals x nSessions] cell array
%                  each non-empty entry is a struct with field:
%                  - pevChunkFull : [1 x nChunks] reconstruction PEV across chunks
%
% NAME-VALUE OPTIONS
%   'AnimalLabels'  : cell array of animal labels, length nAnimals
%   'MarkerSize'    : default 55
%   'JitterWidth'   : default 0.32
%   'LineWidth'     : default 1.0
%   'ShowMeanLine'  : default true
%   'MeanLineWidth' : default 2.5
%
% EXAMPLE
%   plotReconPevAcrossAnimals(rezReconPevC)
%
%   plotReconPevAcrossAnimals(rezReconPevC, ...
%       'AnimalLabels', {'m1044','m1045','m1092','m1094','m1613','m1859','m1873'});

p = inputParser;
p.addParameter('AnimalLabels', {}, @(x) iscell(x) || isstring(x));
p.addParameter('MarkerSize', 55, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('JitterWidth', 0.32, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('LineWidth', 1.0, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('ShowMeanLine', true, @(x) islogical(x) && isscalar(x));
p.addParameter('MeanLineWidth', 2.5, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.parse(varargin{:});

animalLabels = cellstr(p.Results.AnimalLabels);
markerSize = p.Results.MarkerSize;
jitterWidth = p.Results.JitterWidth;
lineWidth = p.Results.LineWidth;
showMeanLine = p.Results.ShowMeanLine;
meanLineWidth = p.Results.MeanLineWidth;

[nAnimals, nSessMax] = size(rezReconPevC); %#ok<NASGU>

% Default animal labels
if isempty(animalLabels)
    animalLabels = arrayfun(@(x) sprintf('Animal %d', x), 1:nAnimals, 'UniformOutput', false);
elseif numel(animalLabels) ~= nAnimals
    error('AnimalLabels must have length equal to size(rezReconPevC,1).');
end

% Pastel colors for 7 animals (or fewer/more if needed)
pastelColors = localPastelColors(max(nAnimals, 7));
pastelColors = pastelColors(1:nAnimals, :);

% Collect session means per animal
sessMeanC = cell(nAnimals,1);
xCoordC = cell(nAnimals,1);

for a = 1:nAnimals
    y = nan(1, size(rezReconPevC,2));

    for s = 1:size(rezReconPevC,2)
        rez = rezReconPevC{a,s};

        if isempty(rez) || ~isstruct(rez) || ~isfield(rez, 'pevChunkFull') || isempty(rez.pevChunkFull)
            continue
        end

        y(s) = mean(rez.pevChunkFull, 'omitnan');
    end

    y = y(~isnan(y));
    sessMeanC{a} = y;

    n = numel(y);
    if n == 0
        xCoordC{a} = [];
    elseif n == 1
        xCoordC{a} = a;
    else
        % Evenly spread within animal bin for maximal visibility
        xCoordC{a} = linspace(a - jitterWidth, a + jitterWidth, n);
    end
end

% Plot
h = figure('Color', 'w');
hold on;

for a = 1:nAnimals
    x = xCoordC{a};
    y = sessMeanC{a};

    if isempty(y)
        continue
    end

    scatter(x, y, markerSize, ...
        'MarkerFaceColor', pastelColors(a,:), ...
        'MarkerEdgeColor', pastelColors(a,:)*0.65, ...
        'LineWidth', lineWidth);

    if showMeanLine
        mu = mean(y, 'omitnan');
        plot([a-0.18 a+0.18], [mu mu], '-', ...
            'Color', pastelColors(a,:)*0.55, ...
            'LineWidth', meanLineWidth);
    end
end

xlim([0.5 nAnimals+0.5]);
xticks(1:nAnimals);
xticklabels(animalLabels);
xtickangle(0);

ylabel('Mean reconstruction PEV across chunks');
xlabel('Animal');
title('Session-wise mean reconstruction PEV by animal');

box off;
set(gca, 'FontSize', 11, 'TickDir', 'out', 'LineWidth', 1);
grid on;
grid(gca, 'minor');
ax = gca;
ax.GridAlpha = 0.12;
ax.MinorGridAlpha = 0.06;

hold off;

end


function C = localPastelColors(n)
% Soft pastel palette; repeats with interpolation if needed

base = [
    0.80 0.88 1.00  % pastel blue
    1.00 0.82 0.86  % pastel pink
    0.81 0.94 0.84  % pastel green
    0.99 0.90 0.75  % pastel peach
    0.86 0.82 0.96  % pastel lavender
    0.78 0.93 0.93  % pastel cyan
    0.98 0.95 0.74  % pastel yellow
    ];

if n <= size(base,1)
    C = base(1:n,:);
else
    xi = linspace(1, size(base,1), n);
    C = interp1(1:size(base,1), base, xi, 'linear');
end

end