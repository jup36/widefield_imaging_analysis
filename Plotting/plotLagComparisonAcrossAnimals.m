function h = plotLagComparisonAcrossAnimals(pevLagC, varargin)
% plotLagComparisonAcrossAnimals
%
% Plots animal-level mean PEV across multiple lag models.
%
% INPUT
%   pevLagC : cell array of length nLags.
%             Each entry is an [nAnimals x nSessions] cell array.
%             Each cell should contain a scalar session-level PEV.
%
% NAME-VALUE OPTIONS
%   'AnimalLabels'  : cell array of animal labels
%   'LagLabels'     : cell array of lag labels
%   'MarkerSize'    : default 70
%   'LineWidth'     : default 2.0
%   'YLabel'        : default 'Mean PEV across sessions'
%   'TitleText'     : default 'PEV comparison across lags'

p = inputParser;
p.addParameter('AnimalLabels', {}, @(x) iscell(x) || isstring(x));
p.addParameter('LagLabels', {}, @(x) iscell(x) || isstring(x));
p.addParameter('MarkerSize', 70, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('LineWidth', 2.0, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('YLabel', 'Mean PEV across sessions', @(x) ischar(x) || isstring(x));
p.addParameter('TitleText', 'PEV comparison across lags', @(x) ischar(x) || isstring(x));
p.parse(varargin{:});

animalLabels = cellstr(p.Results.AnimalLabels);
lagLabels    = cellstr(p.Results.LagLabels);
markerSize   = p.Results.MarkerSize;
lineWidth    = p.Results.LineWidth;
yLabelText   = char(p.Results.YLabel);
titleText    = char(p.Results.TitleText);

nLags = numel(pevLagC);
[nAnimals, ~] = size(pevLagC{1});

if isempty(animalLabels)
    animalLabels = arrayfun(@(x) sprintf('Animal %d', x), ...
        1:nAnimals, 'UniformOutput', false);
elseif numel(animalLabels) ~= nAnimals
    error('AnimalLabels must have length equal to number of animals.');
end

if isempty(lagLabels)
    lagLabels = arrayfun(@(x) sprintf('Lag %d', x), ...
        1:nLags, 'UniformOutput', false);
elseif numel(lagLabels) ~= nLags
    error('LagLabels must have length equal to number of lag conditions.');
end

% animalMeanMat: nAnimals × nLags
animalMeanMat = nan(nAnimals, nLags);

for l = 1:nLags

    pevC = pevLagC{l};

    if size(pevC, 1) ~= nAnimals
        error('All entries in pevLagC must have the same number of animals.');
    end

    for a = 1:nAnimals

        sessionVals = nan(1, size(pevC, 2));

        for s = 1:size(pevC, 2)

            val = pevC{a, s};

            if isempty(val)
                continue
            end

            if isnumeric(val) && isscalar(val)
                sessionVals(s) = val;
            elseif isnumeric(val)
                sessionVals(s) = mean(val(:), 'omitnan');
            end
        end

        animalMeanMat(a, l) = mean(sessionVals, 'omitnan');
    end
end

% Colors for lags
lagColors = localLagColors(nLags);

% Plot
h = figure('Color', 'w');
hold on;

x = 1:nAnimals;

for l = 1:nLags

    y = animalMeanMat(:, l);

    plot(x, y, '-o', ...
        'Color', lagColors(l,:), ...
        'MarkerFaceColor', lagColors(l,:), ...
        'MarkerEdgeColor', lagColors(l,:) * 0.55, ...
        'MarkerSize', sqrt(markerSize), ...
        'LineWidth', lineWidth, ...
        'DisplayName', lagLabels{l});
end

xlim([0.5, nAnimals + 0.5]);
xticks(1:nAnimals);
xticklabels(animalLabels);
xtickangle(0);

ylabel(yLabelText);
xlabel('Animal');
title(titleText);

legend('Location', 'best', 'Box', 'off');

box off;
set(gca, 'FontSize', 11, 'TickDir', 'out', 'LineWidth', 1);

grid on;
grid(gca, 'minor');

ax = gca;
ax.GridAlpha = 0.12;
ax.MinorGridAlpha = 0.06;

hold off;

end


function C = localLagColors(n)
% Distinct readable colors for lag comparison

base = [
    0.20 0.45 0.80  % blue
    0.90 0.45 0.20  % orange
    0.25 0.65 0.35  % green
    0.55 0.35 0.75  % purple
];

if n <= size(base, 1)
    C = base(1:n, :);
else
    xi = linspace(1, size(base, 1), n);
    C = interp1(1:size(base, 1), base, xi, 'linear');
end

end