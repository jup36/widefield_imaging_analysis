function h = plotDynamicStaticMetricAcrossAnimals(metricMotifC, metricStaticC, varargin)
% plotDynamicStaticMetricAcrossAnimals
%
% Generic paired animal-level plot for Dynamic vs Static metrics.
%
% INPUT
%   metricMotifC  : [nAnimals x nSessions] cell array.
%                   Each cell contains either:
%                       - scalar metric value
%                       - 1-row table containing metricName
%
%   metricStaticC : same format as metricMotifC
%
% NAME-VALUE OPTIONS
%   'MetricName'   : table variable to extract if cells contain tables.
%                    Example: 'mean_H_hoyer_sparsity'
%   'AnimalLabels' : cell array of animal labels
%   'CondLabels'   : default {'Dynamic', 'Static'}
%   'MarkerSize'   : default 80
%   'YLabel'       : y-axis label
%   'TitleText'    : plot title
%   'FigSaveLogic' : true/false
%   'SaveDir'      : directory to save
%   'SaveName'     : base filename

p = inputParser;
p.addParameter('MetricName', '', @(x) ischar(x) || isstring(x));
p.addParameter('AnimalLabels', {}, @(x) iscell(x) || isstring(x));
p.addParameter('CondLabels', {'Dynamic', 'Static'}, @(x) iscell(x) || isstring(x));
p.addParameter('MarkerSize', 80, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('YLabel', 'Metric value', @(x) ischar(x) || isstring(x));
p.addParameter('TitleText', 'Dynamic vs static metric', @(x) ischar(x) || isstring(x));
p.addParameter('FigSaveLogic', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('SaveDir', pwd, @(x) ischar(x) || isstring(x));
p.addParameter('SaveName', 'dynamic_static_metric', @(x) ischar(x) || isstring(x));
p.parse(varargin{:});

metricName   = char(p.Results.MetricName);
animalLabels = cellstr(p.Results.AnimalLabels);
condLabels   = cellstr(p.Results.CondLabels);
markerSize   = p.Results.MarkerSize;
yLabelText   = char(p.Results.YLabel);
titleText    = char(p.Results.TitleText);
figSaveLogic = logical(p.Results.FigSaveLogic);
saveDir      = char(p.Results.SaveDir);
saveName     = char(p.Results.SaveName);

if ~isequal(size(metricMotifC), size(metricStaticC))
    error('metricMotifC and metricStaticC must have the same size.');
end

[nAnimals, ~] = size(metricMotifC);

if isempty(animalLabels)
    animalLabels = arrayfun(@(x) sprintf('Animal %d', x), ...
        1:nAnimals, 'UniformOutput', false);
elseif numel(animalLabels) ~= nAnimals
    error('AnimalLabels must have length equal to number of animals.');
end

if numel(condLabels) ~= 2
    error('CondLabels must have two entries.');
end

%% Compute animal-level means across sessions
animalMeanMat = nan(nAnimals, 2);

animalMeanMat(:, 1) = localAnimalMeansFromCells(metricMotifC, metricName);
animalMeanMat(:, 2) = localAnimalMeansFromCells(metricStaticC, metricName);

%% Plot
h = figure('Color', 'w');
hold on;

x = 1:2;

animalColors = localPastelColors(nAnimals);

jitterWidth = 0.08;
rng(1);
jitterA = (rand(nAnimals, 1) - 0.5) * 2 * jitterWidth;

for a = 1:nAnimals

    y = animalMeanMat(a, :);

    if all(isnan(y))
        continue
    end

    xA = x + jitterA(a);

    plot(xA, y, '-', ...
        'Color', animalColors(a,:) * 0.75, ...
        'LineWidth', 1.2, ...
        'HandleVisibility', 'off');

    scatter(xA(1), y(1), markerSize, ...
        'MarkerFaceColor', animalColors(a,:), ...
        'MarkerEdgeColor', animalColors(a,:) * 0.55, ...
        'LineWidth', 1.0, ...
        'MarkerFaceAlpha', 0.85, ...
        'MarkerEdgeAlpha', 0.9, ...
        'HandleVisibility', 'off');

    scatter(xA(2), y(2), markerSize, ...
        'MarkerFaceColor', animalColors(a,:), ...
        'MarkerEdgeColor', animalColors(a,:) * 0.55, ...
        'LineWidth', 1.0, ...
        'MarkerFaceAlpha', 0.85, ...
        'MarkerEdgeAlpha', 0.9, ...
        'HandleVisibility', 'off');
end

% Group mean bars
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

% Animal legend
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
title(titleText, 'Interpreter', 'none');

legend('Location', 'bestoutside', 'Box', 'off');

box off;
set(gca, 'FontSize', 11, 'TickDir', 'out', 'LineWidth', 1);

grid on;
ax = gca;
ax.GridAlpha = 0.12;

hold off;

%% Print summary
fprintf('\n%s\n', titleText);
fprintf('Dynamic mean = %.4f\n', groupMean(1));
fprintf('Static mean  = %.4f\n', groupMean(2));
fprintf('Dynamic - static = %.4f\n', ...
    mean(animalMeanMat(:,1) - animalMeanMat(:,2), 'omitnan'));

%% Save
if figSaveLogic
    if exist(saveDir, 'dir') ~= 7
        mkdir(saveDir);
    end

    saveNameBase = fullfile(saveDir, ...
        sprintf('%s_%s', saveName, datestr(now, 'mmddyy_HHMMSS')));

    savefig(h, [saveNameBase '.fig']);
    print(h, [saveNameBase '.png'], '-dpng', '-r300');

    fprintf('Saved figure:\n%s.fig\n%s.png\n', saveNameBase, saveNameBase);
end

end

%% ------------------------------------------------------------------------
function animalMeans = localAnimalMeansFromCells(metricC, metricName)

[nAnimals, nSessions] = size(metricC);
animalMeans = nan(nAnimals, 1);

for a = 1:nAnimals

    sessionVals = nan(1, nSessions);

    for s = 1:nSessions

        val = metricC{a, s};

        if isempty(val)
            continue
        end

        if istable(val)

            if isempty(metricName)
                error('MetricName must be provided when cells contain tables.');
            end

            if ~ismember(metricName, val.Properties.VariableNames)
                warning('Metric "%s" not found for animal %d session %d.', ...
                    metricName, a, s);
                continue
            end

            metricVal = val.(metricName);

            if isnumeric(metricVal)
                sessionVals(s) = mean(metricVal(:), 'omitnan');
            end

        elseif isnumeric(val)

            if isscalar(val)
                sessionVals(s) = val;
            else
                sessionVals(s) = mean(val(:), 'omitnan');
            end
        end
    end

    animalMeans(a) = mean(sessionVals, 'omitnan');
end

end

%% ------------------------------------------------------------------------
function C = localPastelColors(n)

base = lines(max(n, 7));

pastelStrength = 0.55;
C = base(1:n, :) * pastelStrength + ones(n, 3) * (1 - pastelStrength);

C = min(max(C, 0), 1);

end