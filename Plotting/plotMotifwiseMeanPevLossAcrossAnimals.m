function h = plotMotifwiseMeanPevLossAcrossAnimals(rezReconPevC, lossField, varargin)
% plotMotifwiseMeanPevLossAcrossAnimals
%
% Plot motif-wise mean PEV loss across sessions and animals.
%
% INPUTS
%   rezReconPevC : [nAnimals x nSessions] cell array of session result structs
%   lossField    : 'mean_pevTrialLoss' or 'mean_pevChunkLoss'
%
% NAME-VALUE OPTIONS
%   'AnimalLabels' : cell array of animal labels
%   'MarkerSize'   : default 34
%   'Alpha'        : default 0.75
%   'JitterWidth'  : default 0.32
%   'MotifXOffset' : default 0.18
%   'ShowLegend'   : default true
%   'TitleText'    : default auto
%   'YLabel'       : default 'PEV loss'
%
% EXAMPLE
%   plotMotifwiseMeanPevLossAcrossAnimals(rezReconPevC, 'mean_pevTrialLoss')
%   plotMotifwiseMeanPevLossAcrossAnimals(rezReconPevC, 'mean_pevChunkLoss')
%
% Notes
%   - Each session contributes one 24x1 vector.
%   - Within each animal, horizontal offsets are assigned in session order.
%   - Across animals, a fixed small offset separates the color groups.

p = inputParser;
p.addParameter('AnimalLabels', {}, @(x) iscell(x) || isstring(x));
p.addParameter('MarkerSize', 34, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('Alpha', 0.75, @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 1);
p.addParameter('JitterWidth', 0.32, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('MotifXOffset', 0.18, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('ShowLegend', true, @(x) islogical(x) && isscalar(x));
p.addParameter('TitleText', '', @(x) ischar(x) || isstring(x));
p.addParameter('YLabel', 'PEV loss', @(x) ischar(x) || isstring(x));
p.parse(varargin{:});

animalLabels = cellstr(p.Results.AnimalLabels);
markerSize = p.Results.MarkerSize;
alphaVal = p.Results.Alpha;
jitterWidth = p.Results.JitterWidth;
motifXOffset = p.Results.MotifXOffset;
showLegend = p.Results.ShowLegend;
titleText = char(p.Results.TitleText);
yLabelText = char(p.Results.YLabel);

[nAnimals, ~] = size(rezReconPevC);

if isempty(animalLabels)
    animalLabels = arrayfun(@(x) sprintf('Animal %d', x), 1:nAnimals, 'UniformOutput', false);
elseif numel(animalLabels) ~= nAnimals
    error('AnimalLabels must have length equal to size(rezReconPevC,1).');
end

% infer number of motifs from first valid entry
K = [];
for a = 1:nAnimals
    for s = 1:size(rezReconPevC,2)
        rez = rezReconPevC{a,s};
        if ~isempty(rez) && isstruct(rez) && isfield(rez, lossField) && ~isempty(rez.(lossField))
            K = numel(rez.(lossField));
            break
        end
    end
    if ~isempty(K), break; end
end
if isempty(K)
    error('Could not infer motif count from rezReconPevC and lossField.');
end

% same pastel palette idea as previous plot
pastelColors = localPastelColors(max(nAnimals, 7));
pastelColors = pastelColors(1:nAnimals, :);

h = figure('Color', 'w');
hold on;

hSc = gobjects(nAnimals,1);

% fixed animal offset around each motif center
if nAnimals == 1
    animalOffsets = 0;
else
    animalOffsets = linspace(-motifXOffset, motifXOffset, nAnimals);
end

for a = 1:nAnimals
    % collect valid sessions for this animal in session order
    sessVecs = {};
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
        sessVecs{end+1} = v; %#ok<AGROW>
    end

    nSess = numel(sessVecs);
    if nSess == 0
        continue
    end

    % session-order horizontal offsets, left -> right = early -> late
    if nSess == 1
        sessOffsets = 0;
    else
        sessOffsets = linspace(-jitterWidth/2, jitterWidth/2, nSess);
    end

    xAll = [];
    yAll = [];

    for s = 1:nSess
        v = sessVecs{s};

        x = (1:K) + animalOffsets(a) + sessOffsets(s);
        y = 100 * v';

        xAll = [xAll, x]; %#ok<AGROW>
        yAll = [yAll, y]; %#ok<AGROW>
    end

    hSc(a) = scatter(xAll, yAll, markerSize, ...
        'MarkerFaceColor', pastelColors(a,:), ...
        'MarkerEdgeColor', pastelColors(a,:)*0.65, ...
        'MarkerFaceAlpha', alphaVal, ...
        'MarkerEdgeAlpha', alphaVal, ...
        'DisplayName', animalLabels{a});
end

xlim([0.5, K+0.5]);
xticks(1:K);
xlabel('Motif #');
ylabel(yLabelText);
ylim([0 100])

if isempty(titleText)
    switch lossField
        case 'mean_pevTrialLoss'
            title('Motif-wise mean trial PEV loss');
        case 'mean_pevChunkLoss'
            title('Motif-wise mean chunk PEV loss');
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
    legend(hSc(isgraphics(hSc)), animalLabels(isgraphics(hSc)), ...
        'Location', 'eastoutside', 'Box', 'off');
end
pbaspect([1.5 1 1]); 
hold off;

end


function C = localPastelColors(n)
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