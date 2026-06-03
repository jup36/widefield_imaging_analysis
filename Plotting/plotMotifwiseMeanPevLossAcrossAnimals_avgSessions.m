function [h, pevLossRez] = plotMotifwiseMeanPevLossAcrossAnimals_avgSessions(rezReconPevC, lossField, varargin)
% plotMotifwiseMeanPevLossAcrossAnimals_avgSessions
%
% Plot motif-wise mean PEV loss across animals after averaging across sessions
% within each animal.
%
% INPUTS
%   rezReconPevC : [nAnimals x nSessions] cell array of session result structs
%   lossField    : 'mean_pevTrialLoss' or 'mean_pevChunkLoss'
%
% NAME-VALUE OPTIONS
%   'AnimalLabels' : cell array of animal labels
%   'MarkerSize'   : default 60
%   'Alpha'        : default 0.9
%   'AnimalOffset' : default 0.22
%   'ShowLegend'   : default true
%   'TitleText'    : default auto
%   'YLabel'       : default 'PEV loss (%)'
%   'ConnectLines' : default false
%
% OUTPUTS
%   h          : figure handle
%   pevLossRez : structure containing descriptive statistics across animals
%                for each motif:
%                   .meanLoss         [1 x K]
%                   .stdLoss          [1 x K]
%                   .semLoss          [1 x K]
%                   .meanMat          [nAnimals x K]
%                   .nAnimalsPerMotif [1 x K]
%                   .animalLabels
%                   .lossField
%
% NOTES
%   - Values are converted to percentage by multiplying by 100.
%   - Values > 100 (%) are treated as spurious and set to NaN before
%     averaging across sessions.
%   - One point per animal per motif is plotted after session averaging.
%
% EXAMPLE
%   animalLabels = {'m1044','m1045','m1092','m1094','m1613','m1859','m1873'};
%
%   [h, pevLossRez] = plotMotifwiseMeanPevLossAcrossAnimals_avgSessions( ...
%       rezReconPevC, 'mean_pevTrialLoss', 'AnimalLabels', animalLabels);

p = inputParser;
p.addParameter('AnimalLabels', {}, @(x) iscell(x) || isstring(x));
p.addParameter('MarkerSize', 60, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('Alpha', 0.9, @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 1);
p.addParameter('AnimalOffset', 0.22, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('ShowLegend', true, @(x) islogical(x) && isscalar(x));
p.addParameter('TitleText', '', @(x) ischar(x) || isstring(x));
p.addParameter('YLabel', 'PEV loss (%)', @(x) ischar(x) || isstring(x));
p.addParameter('ConnectLines', false, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});

animalLabels = cellstr(p.Results.AnimalLabels);
markerSize = p.Results.MarkerSize;
alphaVal = p.Results.Alpha;
animalOffset = p.Results.AnimalOffset;
showLegend = p.Results.ShowLegend;
titleText = char(p.Results.TitleText);
yLabelText = char(p.Results.YLabel);
connectLines = p.Results.ConnectLines;

[nAnimals, ~] = size(rezReconPevC);

% Default animal labels
if isempty(animalLabels)
    animalLabels = arrayfun(@(x) sprintf('Animal %d', x), 1:nAnimals, 'UniformOutput', false);
elseif numel(animalLabels) ~= nAnimals
    error('AnimalLabels must have length equal to size(rezReconPevC,1).');
end

% Infer motif count K from first valid entry
K = [];
for a = 1:nAnimals
    for s = 1:size(rezReconPevC,2)
        rez = rezReconPevC{a,s};
        if ~isempty(rez) && isstruct(rez) && isfield(rez, lossField) && ~isempty(rez.(lossField))
            K = numel(rez.(lossField));
            break
        end
    end
    if ~isempty(K)
        break
    end
end

if isempty(K)
    error('Could not infer motif count from rezReconPevC and lossField.');
end

% Same pastel palette style as before
pastelColors = localPastelColors(max(nAnimals, 7));
pastelColors = pastelColors(1:nAnimals, :);

% Fixed animal offsets around each motif center
if nAnimals == 1
    animalOffsets = 0;
else
    animalOffsets = linspace(-animalOffset, animalOffset, nAnimals);
end

% Collect averaged data across sessions for each animal
meanMat = nan(nAnimals, K);

for a = 1:nAnimals
    sessMat = nan(K, size(rezReconPevC,2));

    for s = 1:size(rezReconPevC,2)
        rez = rezReconPevC{a,s};

        if isempty(rez) || ~isstruct(rez) || ~isfield(rez, lossField) || isempty(rez.(lossField))
            continue
        end

        v = rez.(lossField);
        v = v(:);

        if numel(v) ~= K
            warning('Skipping animal %d session %d due to inconsistent motif count.', a, s);
            continue
        end

        % Convert to percent
        v = 100 * v;

        % Treat spurious >100% values as NaN before averaging
        v(v > 100) = NaN;

        sessMat(:, s) = v;
    end

    meanMat(a, :) = mean(sessMat, 2, 'omitnan')';
end

% -------- descriptive stats across animals --------
nAnimalsPerMotif = sum(~isnan(meanMat), 1);
meanLoss = mean(meanMat, 1, 'omitnan');
stdLoss = std(meanMat, 0, 1, 'omitnan');
semLoss = stdLoss ./ sqrt(max(nAnimalsPerMotif, 1));
semLoss(nAnimalsPerMotif == 0) = NaN;

pevLossRez = struct();
pevLossRez.meanLoss = meanLoss;
pevLossRez.stdLoss = stdLoss;
pevLossRez.semLoss = semLoss;
pevLossRez.meanMat = meanMat;
pevLossRez.nAnimalsPerMotif = nAnimalsPerMotif;
pevLossRez.animalLabels = animalLabels;
pevLossRez.lossField = lossField;

% -------- Plot --------
h = figure('Color', 'w');

% Double the current figure size
pos = get(h, 'Position');
set(h, 'Position', [pos(1), pos(2), pos(3)*2, pos(4)*2]);

hold on;

hSc = gobjects(nAnimals, 1);

for a = 1:nAnimals
    x = (1:K) + animalOffsets(a);
    y = meanMat(a, :);

    valid = ~isnan(y);
    if ~any(valid)
        continue
    end

    if connectLines
        plot(x(valid), y(valid), '-', ...
            'Color', pastelColors(a,:) * 0.75, ...
            'LineWidth', 1.2, ...
            'HandleVisibility', 'off');
    end

    hSc(a) = scatter(x(valid), y(valid), markerSize, ...
        'MarkerFaceColor', pastelColors(a,:), ...
        'MarkerEdgeColor', pastelColors(a,:) * 0.6, ...
        'MarkerFaceAlpha', alphaVal, ...
        'MarkerEdgeAlpha', alphaVal, ...
        'LineWidth', 0.9, ...
        'DisplayName', animalLabels{a});
end

xlim([0.5, K + 0.5]);
xticks(1:K);
xlabel('Motif #');
ylabel(yLabelText);

if isempty(titleText)
    switch lossField
        case 'mean_pevTrialLoss'
            title('Motif-wise mean trial PEV loss across animals');
        case 'mean_pevChunkLoss'
            title('Motif-wise mean chunk PEV loss across animals');
        otherwise
            title(strrep(lossField, '_', '\_'));
    end
else
    title(titleText);
end

box off;
set(gca, 'FontSize', 11, 'TickDir', 'out', 'LineWidth', 1);
grid on;
grid(gca, 'minor');
ax = gca;
ax.GridAlpha = 0.12;
ax.MinorGridAlpha = 0.06;

if showLegend
    validHandles = hSc(isgraphics(hSc));
    validLabels = animalLabels(isgraphics(hSc));
    legend(validHandles, validLabels, 'Location', 'eastoutside', 'Box', 'off');
end

pbaspect([2 1 1]);
hold off;

end


function C = localPastelColors(n)
% Soft pastel palette

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
    xi = linspace(1, size(base,1), n);
    C = interp1(1:size(base,1), base, xi, 'linear');
end

end