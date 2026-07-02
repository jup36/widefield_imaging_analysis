function figH = plot_LK_PEV_vs_Hoyer_scatter(pevSummaryT, varargin)
% plot_LK_PEV_vs_Hoyer_scatter
%
% Scatter plot for choosing Lag/K combinations based on reconstruction
% quality and temporal sparsity.
%
% X-axis:
%   mean_H_hoyer_sparsity
%
% Y-axis:
%   meanPEV
%
% Supported plot levels:
%
%   'animal'
%       One dot per animal x L/K condition.
%       Sessions are averaged within animal first.
%
%   'session'
%       One dot per session x L/K condition.
%
%   'condition'
%       One dot per L/K condition.
%       Data are first averaged across sessions within each animal, then
%       averaged across animals. This is the simplest decision plot.
%
% Color:
%   'animal'/'session'  : animal identity
%   'condition'         : L/K condition
%
% Marker:
%   L/K condition
%
% Expected table columns:
%   mouseID
%   sessionID
%   L
%   K_input
%   meanPEV
%   mean_H_hoyer_sparsity
%
% Optional table column:
%   sourceLabel, e.g. 'LKcombo' or 'legacy_L10_K15'

%% Parse inputs

p = inputParser;

addRequired(p, 'pevSummaryT', @istable);

addParameter(p, 'plotLevel', 'animal', ...
    @(x) ischar(x) || isstring(x));

addParameter(p, 'animalLabels', {}, ...
    @(x) isempty(x) || iscell(x) || isstring(x) || iscategorical(x));

addParameter(p, 'metricX', 'mean_H_hoyer_sparsity', ...
    @(x) ischar(x) || isstring(x));

addParameter(p, 'metricY', 'meanPEV', ...
    @(x) ischar(x) || isstring(x));

addParameter(p, 'xLabel', 'Mean Hoyer sparsity', ...
    @(x) ischar(x) || isstring(x));

addParameter(p, 'yLabel', 'Mean PEV', ...
    @(x) ischar(x) || isstring(x));

addParameter(p, 'xLim', [], ...
    @(x) isempty(x) || (isnumeric(x) && numel(x) == 2));

addParameter(p, 'yLim', [], ...
    @(x) isempty(x) || (isnumeric(x) && numel(x) == 2));

addParameter(p, 'showDiagonal', true, ...
    @(x) islogical(x) || isnumeric(x));

addParameter(p, 'titleText', '', ...
    @(x) ischar(x) || isstring(x));

addParameter(p, 'markerSize', 85, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);

addParameter(p, 'lineWidth', 1.0, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);

addParameter(p, 'figSaveLogic', false, ...
    @(x) islogical(x) || isnumeric(x));

addParameter(p, 'saveDir', pwd, ...
    @(x) ischar(x) || isstring(x));

addParameter(p, 'saveName', 'scatter_PEV_vs_Hoyer', ...
    @(x) ischar(x) || isstring(x));

addParameter(p, 'showPointLabels', true, ...
    @(x) islogical(x) || isnumeric(x));

addParameter(p, 'labelBackground', false, ...
    @(x) islogical(x) || isnumeric(x));

addParameter(p, 'labelFontSize', 8, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);

parse(p, pevSummaryT, varargin{:});

plotLevel    = lower(char(p.Results.plotLevel));
animalLabels = p.Results.animalLabels;
metricX      = char(p.Results.metricX);
metricY      = char(p.Results.metricY);
xLabelText   = char(p.Results.xLabel);
yLabelText   = char(p.Results.yLabel);
xLimVals     = p.Results.xLim;
yLimVals     = p.Results.yLim;
showDiagonal = logical(p.Results.showDiagonal);
titleText    = char(p.Results.titleText);
markerSize   = p.Results.markerSize;
lineWidth    = p.Results.lineWidth;
showPointLabels = logical(p.Results.showPointLabels);
labelBackground = logical(p.Results.labelBackground);
labelFontSize   = p.Results.labelFontSize;
figSaveLogic = logical(p.Results.figSaveLogic);
saveDir      = char(p.Results.saveDir);
saveName     = char(p.Results.saveName);

%% Validate table

requiredVars = {'mouseID', 'sessionID', 'L', 'K_input', metricX, metricY};

for i = 1:numel(requiredVars)
    if ~ismember(requiredVars{i}, pevSummaryT.Properties.VariableNames)
        error('Required variable "%s" is missing from pevSummaryT.', requiredVars{i});
    end
end

if ~ismember(plotLevel, {'animal', 'session', 'condition'})
    error('plotLevel must be "animal", "session", or "condition".');
end

%% Clean table

T = pevSummaryT;

T.mouseID_str   = cellstr(string(T.mouseID));
T.sessionID_str = cellstr(string(T.sessionID));

validI = ~isnan(T.(metricX)) & ~isnan(T.(metricY));
T = T(validI, :);

if isempty(T)
    error('No valid rows found after removing NaN X/Y values.');
end

%% Define L/K labels

T.LK_label = arrayfun(@(L,K) sprintf('L%d/K%d', L, K), ...
    T.L, T.K_input, 'UniformOutput', false);

% Sort by L then K for stable condition ordering
[~, sortI] = sortrows([T.L, T.K_input], [1 2]);
T = T(sortI, :);

condLabels = unique(T.LK_label, 'stable');
nConds = numel(condLabels);

%% Define marker list
% Avoid sideways triangle markers '<' and '>' because they look odd when
% points overlap.

markerList = {'o', 's', 'd', '^', 'v', 'p', 'h', 'o', 's', 'd', '^', 'v'};

if nConds > numel(markerList)
    warning('More L/K conditions than marker types. Marker shapes will repeat.');
end

%% Define animal labels and colors

if isempty(animalLabels)
    animalLabels = unique(T.mouseID_str, 'stable');
else
    animalLabels = cellstr(string(animalLabels));
end

nAnimals = numel(animalLabels);
animalColors = localAnimalPastelColors(animalLabels);

%% Aggregate depending on plot level

switch plotLevel
    
    case 'animal'
        
        plotT = localAggregateAnimalLevel(T, metricX, metricY);
        
        defaultTitle = 'PEV vs Hoyer sparsity: one dot per animal per L/K condition';
        
    case 'session'
        
        plotT = T(:, {'mouseID_str', 'sessionID_str', 'L', 'K_input', ...
            'LK_label', metricX, metricY});
        
        defaultTitle = 'PEV vs Hoyer sparsity: one dot per session per L/K condition';
        
    case 'condition'
        
        animalT = localAggregateAnimalLevel(T, metricX, metricY);
        plotT = localAggregateConditionLevel(animalT, metricX, metricY);
        
        defaultTitle = 'PEV vs Hoyer sparsity: across-animal average per L/K condition';
end

if isempty(titleText)
    titleText = defaultTitle;
end

%% Figure and axes

figH = figure('Color', 'w', 'Position', [100 100 780 760]);
ax = axes(figH);
hold(ax, 'on');

% Manual axes position gives more white space above the title and room for
% the external legend.
switch plotLevel
    case 'condition'
        ax.Position = [0.14 0.14 0.70 0.68];
    otherwise
        ax.Position = [0.12 0.14 0.58 0.68];
end

%% Plot points

switch plotLevel
    
    case {'animal', 'session'}
        
        for a = 1:nAnimals
            
            animalID = animalLabels{a};
            animalI = strcmp(plotT.mouseID_str, animalID);
            
            if ~any(animalI)
                continue
            end
            
            thisColor = animalColors(a, :);
            %% Add optional labels for L/K condition centers
            
            if showPointLabels
                
                for c = 1:nConds
                    
                    condID = condLabels{c};
                    condI = strcmp(plotT.LK_label, condID);
                    
                    if ~any(condI)
                        continue
                    end
                    
                    xCenter = mean(plotT.(metricX)(condI), 'omitnan');
                    yCenter = mean(plotT.(metricY)(condI), 'omitnan');
                    
                    if labelBackground
                        text(xCenter, yCenter, ['  ' condID], ...
                            'FontSize', labelFontSize, ...
                            'FontWeight', 'bold', ...
                            'Color', [0.15 0.15 0.15], ...
                            'BackgroundColor', 'none', ...
                            'Margin', 1, ...
                            'Interpreter', 'none');
                    else
                        % Offset label by a fraction of the current axis range
                        xRange = diff(xlim);
                        yRange = diff(ylim);
                        
                        xOffset = 0.018 * xRange;
                        yOffset = 0.010 * yRange;
                        
                        text(xCenter + xOffset, yCenter + yOffset, condID, ...
                            'FontSize', labelFontSize, ...
                            'FontWeight', 'bold', ...
                            'Color', [0.15 0.15 0.15], ...
                            'BackgroundColor', 'none', ...
                            'HorizontalAlignment', 'left', ...
                            'VerticalAlignment', 'bottom', ...
                            'Interpreter', 'none', ...
                            'Clipping', 'off');
                    end
                end
            end
        end
        
    case 'condition'
        
        condColors = localConditionPastelColors(nConds);
        
        for c = 1:nConds
            
            condID = condLabels{c};
            rowI = strcmp(plotT.LK_label, condID);
            
            if ~any(rowI)
                continue
            end
            
            markerThis = markerList{mod(c-1, numel(markerList)) + 1};
            
            scatter(ax, plotT.(metricX)(rowI), plotT.(metricY)(rowI), ...
                markerSize * 1.25, ...
                'Marker', markerThis, ...
                'MarkerFaceColor', condColors(c, :), ...
                'MarkerEdgeColor', condColors(c, :) .* 0.55, ...
                'LineWidth', lineWidth + 0.3, ...
                'MarkerFaceAlpha', 0.85, ...
                'MarkerEdgeAlpha', 0.95, ...
                'DisplayName', condID);
            
            % Label each condition directly on the plot.
            text(ax, plotT.(metricX)(rowI), plotT.(metricY)(rowI), ...
                ['  ' condID], ...
                'FontSize', 9, ...
                'FontWeight', 'bold', ...
                'BackgroundColor', 'none', ...
                'Color', [0.15 0.15 0.15], ...
                'Margin', 1, ...
                'Interpreter', 'none');
        end
end

%% Add L/K condition center labels for animal/session plots

if ~strcmp(plotLevel, 'condition')
    
    for c = 1:nConds
        
        condID = condLabels{c};
        condI = strcmp(plotT.LK_label, condID);
        
        if ~any(condI)
            continue
        end
        
        xCenter = mean(plotT.(metricX)(condI), 'omitnan');
        yCenter = mean(plotT.(metricY)(condI), 'omitnan');
        
        text(ax, xCenter, yCenter, ['  ' condID], ...
            'FontSize', 9, ...
            'FontWeight', 'bold', ...
            'BackgroundColor', 'none', ...
            'Color', [0.15 0.15 0.15], ...
            'Margin', 1, ...
            'Interpreter', 'none');
    end
end

%% Axes labels and title

xlabel(ax, xLabelText);
ylabel(ax, yLabelText);

titleH = title(ax, titleText, ...
    'Interpreter', 'none', ...
    'FontWeight', 'bold');

% Push title slightly upward inside available white space.
try
    titleH.Units = 'normalized';
    titleH.Position(2) = 1.07;
catch
end

%% Axis limits

if ~isempty(xLimVals)
    xlim(ax, xLimVals);
else
    xlim(ax, [0 1]);
end

if ~isempty(yLimVals)
    ylim(ax, yLimVals);
else
    yVals = plotT.(metricY);
    
    if all(yVals >= 0 & yVals <= 1)
        ylim(ax, [0 1]);
    end
end

%% Diagonal reference line: x = y

if showDiagonal
    
    xl = xlim(ax);
    yl = ylim(ax);
    
    diagMin = max(xl(1), yl(1));
    diagMax = min(xl(2), yl(2));
    
    if diagMax > diagMin
        
        hDiag = plot(ax, [diagMin diagMax], [diagMin diagMax], ...
            '--', ...
            'Color', [0.25 0.25 0.25], ...
            'LineWidth', 1.2, ...
            'HandleVisibility', 'off');
        
        try
            uistack(hDiag, 'bottom');
        catch
        end
        
        text(ax, diagMax, diagMax, '  x = y', ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 9, ...
            'Color', [0.25 0.25 0.25], ...
            'Interpreter', 'none', ...
            'HandleVisibility', 'off');
    end
end

%% Formatting

box(ax, 'off');
grid(ax, 'on');

ax.FontSize = 11;
ax.TickDir = 'out';
ax.LineWidth = 1;
ax.GridAlpha = 0.12;

pbaspect(ax, [1 1 1]);

%% Legends

switch plotLevel
    
    case {'animal', 'session'}
        
        % Animal color legend
        for a = 1:nAnimals
            scatter(ax, nan, nan, markerSize, ...
                'Marker', 'o', ...
                'MarkerFaceColor', animalColors(a,:), ...
                'MarkerEdgeColor', animalColors(a,:) .* 0.55, ...
                'LineWidth', lineWidth, ...
                'DisplayName', animalLabels{a});
        end
        
        % L/K marker legend
        for c = 1:nConds
            
            markerThis = markerList{mod(c-1, numel(markerList)) + 1};
            
            scatter(ax, nan, nan, markerSize, ...
                'Marker', markerThis, ...
                'MarkerFaceColor', [0.80 0.80 0.80], ...
                'MarkerEdgeColor', [0.35 0.35 0.35], ...
                'LineWidth', lineWidth, ...
                'DisplayName', condLabels{c});
        end
        
        legend(ax, 'Location', 'bestoutside', 'Box', 'off');
        
    case 'condition'
        
        legend(ax, 'Location', 'bestoutside', 'Box', 'off');
end

hold(ax, 'off');

%% Print summary

fprintf('\n%s\n', titleText);
fprintf('Plot level: %s\n', plotLevel);
fprintf('Number of plotted rows: %d\n', height(plotT));

condSummaryT = localConditionSummary(plotT, metricX, metricY);
disp(condSummaryT);

%% Save

if figSaveLogic
    
    if exist(saveDir, 'dir') ~= 7
        mkdir(saveDir);
    end
    
    saveBase = fullfile(saveDir, ...
        sprintf('%s_%s', saveName, datestr(now, 'mmddyy_HHMMSS')));
    
    savefig(figH, [saveBase '.fig']);
    print(figH, [saveBase '.png'], '-dpng', '-r300');
    
    fprintf('Saved scatter figure:\n%s.fig\n%s.png\n', saveBase, saveBase);
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotT = localAggregateAnimalLevel(T, metricX, metricY)
% Aggregate rows to one value per animal x L/K condition.
%
% This averages across sessions within each animal.

mouseVals = unique(T.mouseID_str, 'stable');
condVals = unique(T.LK_label, 'stable');

rows = {};

for a = 1:numel(mouseVals)
    
    mouseID = mouseVals{a};
    
    for c = 1:numel(condVals)
        
        condID = condVals{c};
        
        rowI = strcmp(T.mouseID_str, mouseID) & strcmp(T.LK_label, condID);
        
        if ~any(rowI)
            continue
        end
        
        Lvals = T.L(rowI);
        Kvals = T.K_input(rowI);
        
        rows(end+1, :) = { ...
            mouseID, ...
            condID, ...
            Lvals(1), ...
            Kvals(1), ...
            mean(T.(metricX)(rowI), 'omitnan'), ...
            mean(T.(metricY)(rowI), 'omitnan'), ...
            sum(rowI) ...
            };
    end
end

plotT = cell2table(rows, ...
    'VariableNames', { ...
    'mouseID_str', ...
    'LK_label', ...
    'L', ...
    'K_input', ...
    metricX, ...
    metricY, ...
    'nSessions' ...
    });

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotT = localAggregateConditionLevel(animalT, metricX, metricY)
% Aggregate animal-level values to one value per L/K condition.
%
% This avoids overweighting animals with more sessions.

condVals = unique(animalT.LK_label, 'stable');

rows = {};

for c = 1:numel(condVals)
    
    condID = condVals{c};
    rowI = strcmp(animalT.LK_label, condID);
    
    if ~any(rowI)
        continue
    end
    
    Lvals = animalT.L(rowI);
    Kvals = animalT.K_input(rowI);
    
    xVals = animalT.(metricX)(rowI);
    yVals = animalT.(metricY)(rowI);
    
    rows(end+1, :) = { ...
        condID, ...
        Lvals(1), ...
        Kvals(1), ...
        mean(xVals, 'omitnan'), ...
        mean(yVals, 'omitnan'), ...
        std(xVals, 'omitnan') ./ sqrt(sum(~isnan(xVals))), ...
        std(yVals, 'omitnan') ./ sqrt(sum(~isnan(yVals))), ...
        sum(rowI) ...
        };
end

plotT = cell2table(rows, ...
    'VariableNames', { ...
    'LK_label', ...
    'L', ...
    'K_input', ...
    metricX, ...
    metricY, ...
    ['sem_' metricX], ...
    ['sem_' metricY], ...
    'nAnimals' ...
    });

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function condSummaryT = localConditionSummary(plotT, metricX, metricY)
% Summarize each L/K condition across plotted rows.

condVals = unique(plotT.LK_label, 'stable');

rows = {};

for c = 1:numel(condVals)
    
    condID = condVals{c};
    rowI = strcmp(plotT.LK_label, condID);
    
    if ~any(rowI)
        continue
    end
    
    rows(end+1, :) = { ...
        condID, ...
        mean(plotT.(metricX)(rowI), 'omitnan'), ...
        std(plotT.(metricX)(rowI), 'omitnan'), ...
        mean(plotT.(metricY)(rowI), 'omitnan'), ...
        std(plotT.(metricY)(rowI), 'omitnan'), ...
        sum(rowI) ...
        };
end

condSummaryT = cell2table(rows, ...
    'VariableNames', { ...
    'LK_condition', ...
    ['mean_' metricX], ...
    ['std_' metricX], ...
    ['mean_' metricY], ...
    ['std_' metricY], ...
    'nPoints' ...
    });

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C = localAnimalPastelColors(animalLabels)
% Stable pastel colors by animal ID.
%
% These match the usual Mesomapping cohort order.

animalLabels = cellstr(string(animalLabels));
nAnimals = numel(animalLabels);

C = nan(nAnimals, 3);

knownAnimalIDs = { ...
    'm1044', ...
    'm1045', ...
    'm1048', ...
    'm1049', ...
    'm1092', ...
    'm1094', ...
    'm1613', ...
    'm1859', ...
    'm1873' ...
    };

knownBaseColors = [ ...
    0.1216 0.4667 0.7059;  % blue
    1.0000 0.4980 0.0549;  % orange
    0.1725 0.6275 0.1725;  % green
    0.8392 0.1529 0.1569;  % red
    0.5804 0.4039 0.7412;  % purple
    0.5490 0.3373 0.2941;  % brown
    0.8902 0.4667 0.7608;  % pink
    0.4980 0.4980 0.4980;  % gray
    0.7373 0.7412 0.1333   % olive
    ];

for a = 1:nAnimals
    
    idx = find(strcmp(knownAnimalIDs, animalLabels{a}), 1);
    
    if ~isempty(idx)
        C(a, :) = knownBaseColors(idx, :);
    end
end

missingI = any(isnan(C), 2);

if any(missingI)
    fallbackColors = lines(sum(missingI));
    C(missingI, :) = fallbackColors;
end

pastelStrength = 0.55;
C = C .* pastelStrength + ones(size(C)) .* (1 - pastelStrength);

C = min(max(C, 0), 1);

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C = localConditionPastelColors(nConds)
% Pastel colors for L/K condition-level plot.

base = lines(max(nConds, 7));

pastelStrength = 0.60;
C = base(1:nConds, :) .* pastelStrength + ones(nConds, 3) .* (1 - pastelStrength);

C = min(max(C, 0), 1);

end