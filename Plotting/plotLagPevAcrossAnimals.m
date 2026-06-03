function h = plotLagPevAcrossAnimals(pevC, varargin)
% plotLagPevAcrossAnimals
%
% Scatter plot of session-wise mean PEV values grouped by animal.
%
% INPUT
%   pevC : [nAnimals x nSessions] cell array
%          each non-empty entry should contain a scalar numeric PEV value.
%
% NAME-VALUE OPTIONS
%   'AnimalLabels'  : cell array of animal labels, length nAnimals
%   'MarkerSize'    : default 55
%   'JitterWidth'   : default 0.32
%   'LineWidth'     : default 1.0
%   'ShowMeanLine'  : default true
%   'MeanLineWidth' : default 2.5
%   'YLabel'        : default 'Mean PEV across chunks'
%   'TitleText'     : default 'Session-wise mean PEV by animal'
%
% EXAMPLE
%   plotLagPevAcrossAnimals(stats_testPevL1C, ...
%       'AnimalLabels', {'m1044','m1045','m1048','m1049','m1613','m1859','m1873'});

p = inputParser;
p.addParameter('AnimalLabels', {}, @(x) iscell(x) || isstring(x));
p.addParameter('MarkerSize', 55, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('JitterWidth', 0.32, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('LineWidth', 1.0, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('ShowMeanLine', true, @(x) islogical(x) && isscalar(x));
p.addParameter('MeanLineWidth', 2.5, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('YLabel', 'Mean PEV across chunks', @(x) ischar(x) || isstring(x));
p.addParameter('TitleText', 'Session-wise mean PEV by animal', @(x) ischar(x) || isstring(x));
p.parse(varargin{:});

animalLabels  = cellstr(p.Results.AnimalLabels);
markerSize    = p.Results.MarkerSize;
jitterWidth   = p.Results.JitterWidth;
lineWidth     = p.Results.LineWidth;
showMeanLine  = p.Results.ShowMeanLine;
meanLineWidth = p.Results.MeanLineWidth;
yLabelText    = char(p.Results.YLabel);
titleText     = char(p.Results.TitleText);

[nAnimals, ~] = size(pevC);

% Default animal labels
if isempty(animalLabels)
    animalLabels = arrayfun(@(x) sprintf('Animal %d', x), ...
        1:nAnimals, 'UniformOutput', false);
elseif numel(animalLabels) ~= nAnimals
    error('AnimalLabels must have length equal to size(pevC,1).');
end

% Pastel colors
pastelColors = localPastelColors(max(nAnimals, 7));
pastelColors = pastelColors(1:nAnimals, :);

% Collect session values per animal
sessValC = cell(nAnimals, 1);
xCoordC  = cell(nAnimals, 1);

for a = 1:nAnimals

    y = nan(1, size(pevC, 2));

    for s = 1:size(pevC, 2)

        val = pevC{a, s};

        if isempty(val)
            continue
        end

        if isnumeric(val) && isscalar(val)
            y(s) = val;
        elseif isnumeric(val)
            % Safety fallback: if a vector sneaks in, average it.
            y(s) = mean(val(:), 'omitnan');
        end
    end

    y = y(~isnan(y));
    sessValC{a} = y;

    n = numel(y);

    if n == 0
        xCoordC{a} = [];
    elseif n == 1
        xCoordC{a} = a;
    else
        xCoordC{a} = linspace(a - jitterWidth, a + jitterWidth, n);
    end
end

% Plot
h = figure('Color', 'w');
hold on;

for a = 1:nAnimals

    x = xCoordC{a};
    y = sessValC{a};

    if isempty(y)
        continue
    end

    scatter(x, y, markerSize, ...
        'MarkerFaceColor', pastelColors(a,:), ...
        'MarkerEdgeColor', pastelColors(a,:) * 0.65, ...
        'LineWidth', lineWidth);

    if showMeanLine
        mu = mean(y, 'omitnan');

        plot([a - 0.18, a + 0.18], [mu, mu], '-', ...
            'Color', pastelColors(a,:) * 0.55, ...
            'LineWidth', meanLineWidth);
    end
end

xlim([0.5, nAnimals + 0.5]);
xticks(1:nAnimals);
xticklabels(animalLabels);
xtickangle(0);

ylabel(yLabelText);
xlabel('Animal');
title(titleText);

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

if n <= size(base, 1)
    C = base(1:n, :);
else
    xi = linspace(1, size(base, 1), n);
    C = interp1(1:size(base, 1), base, xi, 'linear');
end

end