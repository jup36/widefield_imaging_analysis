function h = plotDynamicStaticPEVAcrossAnimals(pevMotifC, pevStaticC, varargin)
% plotDynamicStaticPEVAcrossAnimals
%
% Plots animal-level mean PEV across sessions for dynamic/motif-based
% refits versus static-network refits.
%
% INPUT
%   pevMotifC  : [nAnimals x nSessions] cell array.
%                Each cell contains a scalar session-level motif PEV.
%
%   pevStaticC : [nAnimals x nSessions] cell array.
%                Each cell contains a scalar session-level static PEV.
%
% NAME-VALUE OPTIONS
%   'AnimalLabels' : cell array of animal labels
%   'CondLabels'   : default {'Dynamic motif', 'Static motif'}
%   'MarkerSize'   : default 80
%   'LineWidth'    : default 1.5
%   'YLabel'       : default 'Mean PEV across sessions'
%   'TitleText'    : default 'Dynamic vs static motif reconstruction PEV'

p = inputParser;
p.addParameter('AnimalLabels', {}, @(x) iscell(x) || isstring(x));
p.addParameter('CondLabels', {'Dynamic motif', 'Static motif'}, @(x) iscell(x) || isstring(x));
p.addParameter('MarkerSize', 80, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('LineWidth', 1.5, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('YLabel', 'Mean PEV across sessions', @(x) ischar(x) || isstring(x));
p.addParameter('TitleText', 'Dynamic vs static motif reconstruction PEV', @(x) ischar(x) || isstring(x));
p.parse(varargin{:});

animalLabels = cellstr(p.Results.AnimalLabels);
condLabels   = cellstr(p.Results.CondLabels);
markerSize   = p.Results.MarkerSize;
lineWidth    = p.Results.LineWidth;
yLabelText   = char(p.Results.YLabel);
titleText    = char(p.Results.TitleText);

if ~isequal(size(pevMotifC), size(pevStaticC))
    error('pevMotifC and pevStaticC must have the same size.');
end

[nAnimals, ~] = size(pevMotifC);

if isempty(animalLabels)
    animalLabels = arrayfun(@(x) sprintf('Animal %d', x), ...
        1:nAnimals, 'UniformOutput', false);
elseif numel(animalLabels) ~= nAnimals
    error('AnimalLabels must have length equal to number of animals.');
end

if numel(condLabels) ~= 2
    error('CondLabels must have two entries: dynamic/motif and static.');
end

%% Compute animal-level means across sessions

animalMeanMat = nan(nAnimals, 2);
animalMeanMat(:, 1) = localAnimalMeans(pevMotifC);
animalMeanMat(:, 2) = localAnimalMeans(pevStaticC);

%% Plot

h = figure('Color', 'w');
hold on;

x = 1:2;

% Pastel animal colors
animalColors = localPastelColors(nAnimals);

% Horizontal jitter
jitterWidth = 0.08;
rng(1);  % reproducible jitter
jitterA = (rand(nAnimals, 1) - 0.5) * 2 * jitterWidth;

% Plot paired animal lines and points
for a = 1:nAnimals

    y = animalMeanMat(a, :);

    if all(isnan(y))
        continue
    end

    xA = x + jitterA(a);

    % Paired animal line
    plot(xA, y, '-', ...
        'Color', animalColors(a,:) * 0.75, ...
        'LineWidth', 1.2, ...
        'HandleVisibility', 'off');

    % Dynamic point
    scatter(xA(1), y(1), markerSize, ...
        'MarkerFaceColor', animalColors(a,:), ...
        'MarkerEdgeColor', animalColors(a,:) * 0.55, ...
        'LineWidth', 1.0, ...
        'MarkerFaceAlpha', 0.85, ...
        'MarkerEdgeAlpha', 0.9, ...
        'HandleVisibility', 'off');

    % Static point
    scatter(xA(2), y(2), markerSize, ...
        'MarkerFaceColor', animalColors(a,:), ...
        'MarkerEdgeColor', animalColors(a,:) * 0.55, ...
        'LineWidth', 1.0, ...
        'MarkerFaceAlpha', 0.85, ...
        'MarkerEdgeAlpha', 0.9, ...
        'HandleVisibility', 'off');
end

% Group means as horizontal gray bars
groupMean = mean(animalMeanMat, 1, 'omitnan');

barHalfWidth = 0.18;

plot([x(1)-barHalfWidth, x(1)+barHalfWidth], ...
     [groupMean(1), groupMean(1)], ...
     '-', ...
     'Color', [0.35 0.35 0.35], ...
     'LineWidth', 3, ...
     'DisplayName', 'Mean');

plot([x(2)-barHalfWidth, x(2)+barHalfWidth], ...
     [groupMean(2), groupMean(2)], ...
     '-', ...
     'Color', [0.35 0.35 0.35], ...
     'LineWidth', 3, ...
     'HandleVisibility', 'off');

% Optional: animal legend using invisible points
for a = 1:nAnimals
    scatter(nan, nan, markerSize, ...
        'MarkerFaceColor', animalColors(a,:), ...
        'MarkerEdgeColor', animalColors(a,:) * 0.55, ...
        'DisplayName', animalLabels{a});
end

xlim([0.5, 2.5]);
xticks(x);
xticklabels(condLabels);

ylabel(yLabelText);
title(titleText);

legend('Location', 'bestoutside', 'Box', 'off');

box off;
set(gca, 'FontSize', 11, 'TickDir', 'out', 'LineWidth', 1);

grid on;
ax = gca;
ax.GridAlpha = 0.12;

hold off;

%% Print summary

fprintf('\nDynamic/motif mean PEV = %.4f\n', groupMean(1));
fprintf('Static mean PEV        = %.4f\n', groupMean(2));
fprintf('Mean dynamic - static  = %.4f\n', ...
    mean(animalMeanMat(:,1) - animalMeanMat(:,2), 'omitnan'));

end

function animalMeans = localAnimalMeans(pevC)
% Compute mean across sessions for each animal, ignoring empty cells/NaNs.

[nAnimals, nSessions] = size(pevC);
animalMeans = nan(nAnimals, 1);

for a = 1:nAnimals

    sessionVals = nan(1, nSessions);

    for s = 1:nSessions

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

    animalMeans(a) = mean(sessionVals, 'omitnan');
end

end

function C = localPastelColors(n)
% localPastelColors
%
% Generate n readable pastel colors.

base = lines(max(n, 7));

% Mix with white to make pastel
pastelStrength = 0.55;
C = base(1:n, :) * pastelStrength + ones(n, 3) * (1 - pastelStrength);

% Keep values in valid RGB range
C = min(max(C, 0), 1);

end
