function [figH, lambdaSummaryByAnimal, lambdaSummaryBySession] = plot_legacy_lambda_byAnimal(pevChunkT_legacy, varargin)
% plot_legacy_lambda_byAnimal
%
% Visualize fitted lambda values from legacy/original motif chunks.
%
% Each dot = one valid chunk-level lambda value.
% X-axis   = animal
% Color    = animal
% Jitter   = session-aware jitter within animal
%
% Required input:
%   pevChunkT_legacy
%       Chunk-level table containing at least:
%           mouseID
%           sessionID
%           lambda_train
%
% Optional name-value inputs:
%   'lambdaVar'
%       Name of lambda variable. Default = 'lambda_train'
%
%   'mouseVar'
%       Name of mouse/animal variable. Default = 'mouseID'
%
%   'sessionVar'
%       Name of session variable. Default = 'sessionID'
%
%   'animalOrder'
%       Optional animal order. Default = unique order in table.
%
%   'lambdaRef'
%       Reference lambda line. Default = 5e-4.
%
%   'showLambdaRef'
%       Whether to show reference line. Default = true.
%
%   'yScale'
%       Y-axis scale. Options: 'log' or 'linear'. Default = 'log'.
%
%   'dotSize'
%       Dot size. Default = 28.
%
%   'jitterWidth'
%       Width of within-animal jitter. Default = 0.28.
%
%   'figSaveLogic'
%       Whether to save figure. Default = false.
%
%   'saveDir'
%       Directory to save figure. Default = pwd.
%
%   'savePrefix'
%       Save filename prefix. Default = 'legacy_lambda_byAnimal'
%
%   'titleText'
%       Figure title. Default = '\lambda by Animal and Session'
%
%   'showLegend'
%       Whether to show legend. Default = true.
%
%   'yLim'
%       Optional y-axis limits. Default = [].
%
%   'rngSeed'
%       Random seed for reproducible jitter. Default = 1.
%
% Outputs:
%   figH
%   lambdaSummaryByAnimal
%   lambdaSummaryBySession

%% ------------------------------------------------------------------------
% Parse inputs
% -------------------------------------------------------------------------

p = inputParser;

addRequired(p, 'pevChunkT_legacy', @istable);

addParameter(p, 'lambdaVar', 'lambda_train', ...
    @(x) ischar(x) || isstring(x));

addParameter(p, 'mouseVar', 'mouseID', ...
    @(x) ischar(x) || isstring(x));

addParameter(p, 'sessionVar', 'sessionID', ...
    @(x) ischar(x) || isstring(x));

addParameter(p, 'animalOrder', {}, ...
    @(x) isempty(x) || iscell(x) || isstring(x) || iscategorical(x));

addParameter(p, 'lambdaRef', 5e-4, ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x)));

addParameter(p, 'showLambdaRef', true, ...
    @(x) islogical(x) || isnumeric(x));

addParameter(p, 'yScale', 'log', ...
    @(x) ischar(x) || isstring(x));

addParameter(p, 'dotSize', 28, ...
    @(x) isnumeric(x) && isscalar(x));

addParameter(p, 'jitterWidth', 0.28, ...
    @(x) isnumeric(x) && isscalar(x));

addParameter(p, 'figSaveLogic', false, ...
    @(x) islogical(x) || isnumeric(x));

addParameter(p, 'saveDir', pwd, ...
    @(x) ischar(x) || isstring(x));

addParameter(p, 'savePrefix', 'legacy_lambda_byAnimal', ...
    @(x) ischar(x) || isstring(x));

addParameter(p, 'titleText', '\lambda by Animal and Session', ...
    @(x) ischar(x) || isstring(x));

addParameter(p, 'showLegend', true, ...
    @(x) islogical(x) || isnumeric(x));

addParameter(p, 'yLim', [], ...
    @(x) isempty(x) || (isnumeric(x) && numel(x) == 2));

addParameter(p, 'rngSeed', 1, ...
    @(x) isnumeric(x) && isscalar(x));

parse(p, pevChunkT_legacy, varargin{:});

lambdaVar      = char(p.Results.lambdaVar);
mouseVar       = char(p.Results.mouseVar);
sessionVar     = char(p.Results.sessionVar);
animalOrder    = p.Results.animalOrder;
lambdaRef      = p.Results.lambdaRef;
showLambdaRef  = logical(p.Results.showLambdaRef);
yScale         = lower(char(p.Results.yScale));
dotSize        = p.Results.dotSize;
jitterWidth    = p.Results.jitterWidth;
figSaveLogic   = logical(p.Results.figSaveLogic);
saveDir        = char(p.Results.saveDir);
savePrefix     = char(p.Results.savePrefix);
titleText      = char(p.Results.titleText);
showLegend     = logical(p.Results.showLegend);
yLimUser       = p.Results.yLim;
rngSeed        = p.Results.rngSeed;

if ~ismember(yScale, {'log', 'linear'})
    error('yScale must be either ''log'' or ''linear''.');
end

%% ------------------------------------------------------------------------
% Validate table
% -------------------------------------------------------------------------

requiredVars = {mouseVar, sessionVar, lambdaVar};
missingVars = setdiff(requiredVars, pevChunkT_legacy.Properties.VariableNames);

if ~isempty(missingVars)
    error('pevChunkT_legacy is missing required variable(s): %s', ...
        strjoin(missingVars, ', '));
end

lambdaT = pevChunkT_legacy;

%% ------------------------------------------------------------------------
% Keep valid lambda values
% -------------------------------------------------------------------------

lambdaValsAll = lambdaT.(lambdaVar);

validI = ~isnan(lambdaValsAll) & ...
          isfinite(lambdaValsAll);

% For log scale, lambda must be positive.
% For linear scale, keep zero if it exists, but lambda should generally be positive.
if strcmpi(yScale, 'log')
    validI = validI & lambdaValsAll > 0;
else
    validI = validI & lambdaValsAll >= 0;
end

lambdaT = lambdaT(validI, :);

if isempty(lambdaT)
    error('No valid %s values found for yScale = %s.', lambdaVar, yScale);
end

%% ------------------------------------------------------------------------
% Convert labels to strings
% -------------------------------------------------------------------------

mouseStr   = string(lambdaT.(mouseVar));
sessionStr = string(lambdaT.(sessionVar));

if isempty(animalOrder)
    animalLabels = unique(mouseStr, 'stable');
else
    if iscategorical(animalOrder)
        animalLabels = string(cellstr(animalOrder));
    elseif iscell(animalOrder)
        animalLabels = string(animalOrder);
    else
        animalLabels = string(animalOrder);
    end
end

% Keep only animals that actually exist in the table
animalLabels = animalLabels(ismember(animalLabels, unique(mouseStr)));

nAnimals = numel(animalLabels);

if nAnimals == 0
    error('No animals found after applying animalOrder.');
end

%% ------------------------------------------------------------------------
% Colors
% -------------------------------------------------------------------------

animalColors = lines(nAnimals);

%% ------------------------------------------------------------------------
% Plot
% -------------------------------------------------------------------------

figH = figure('Color', 'w', ...
    'Name', 'Legacy fitted lambda by animal', ...
    'Position', [200 200 760 620]);

ax = axes(figH);
hold(ax, 'on');

rng(rngSeed);

legendH = gobjects(nAnimals, 1);

for ai = 1:nAnimals
    
    thisAnimal = animalLabels(ai);
    animalI = mouseStr == thisAnimal;
    
    thisLambda  = lambdaT.(lambdaVar)(animalI);
    thisSession = sessionStr(animalI);
    
    sessionLabels = unique(thisSession, 'stable');
    nSess = numel(sessionLabels);
    
    xAll = nan(numel(thisLambda), 1);
    
    for si = 1:nSess
        
        sessI = thisSession == sessionLabels(si);
        nThis = sum(sessI);
        
        if nSess == 1
            sessCenterOffset = 0;
        else
            sessCenterOffsets = linspace(-jitterWidth/2, jitterWidth/2, nSess);
            sessCenterOffset = sessCenterOffsets(si);
        end
        
        localJitter = (rand(nThis, 1) - 0.5) * jitterWidth * 0.45;
        xAll(sessI) = ai + sessCenterOffset + localJitter;
    end
    
    h = scatter(ax, xAll, thisLambda, dotSize, ...
        'MarkerFaceColor', animalColors(ai, :), ...
        'MarkerEdgeColor', animalColors(ai, :), ...
        'MarkerFaceAlpha', 0.65, ...
        'MarkerEdgeAlpha', 0.65);
    
    legendH(ai) = h;
    
    % Overlay animal median
    medLambda = median(thisLambda, 'omitnan');
    
    plot(ax, [ai - 0.32, ai + 0.32], [medLambda, medLambda], ...
        '-', ...
        'Color', animalColors(ai, :), ...
        'LineWidth', 2.5);
end

%% ------------------------------------------------------------------------
% Reference line
% -------------------------------------------------------------------------

if showLambdaRef && ~isempty(lambdaRef) && isfinite(lambdaRef) && lambdaRef >= 0
    
    yline(ax, lambdaRef, 'r:', ...
        sprintf('\\lambda = %.1g', lambdaRef), ...
        'LineWidth', 1.5, ...
        'LabelHorizontalAlignment', 'left', ...
        'LabelVerticalAlignment', 'bottom');
end

%% ------------------------------------------------------------------------
% Axis formatting
% -------------------------------------------------------------------------

set(ax, ...
    'YScale', yScale, ...
    'XTick', 1:nAnimals, ...
    'XTickLabel', cellstr(animalLabels), ...
    'TickDir', 'out', ...
    'Box', 'off', ...
    'FontSize', 11);

xtickangle(ax, 45);

xlabel(ax, 'Animal');
ylabel(ax, '\lambda');
title(ax, titleText, ...
    'FontWeight', 'bold', ...
    'Interpreter', 'tex');

grid(ax, 'on');
ax.GridAlpha = 0.15;

if showLegend
    legend(ax, legendH, cellstr(animalLabels), ...
        'Location', 'eastoutside', ...
        'Interpreter', 'none');
end

%% ------------------------------------------------------------------------
% Y-limits
% -------------------------------------------------------------------------

lambdaVals = lambdaT.(lambdaVar);

if isempty(yLimUser)
    
    yMin = min(lambdaVals, [], 'omitnan');
    yMax = max(lambdaVals, [], 'omitnan');
    
    if showLambdaRef && ~isempty(lambdaRef) && isfinite(lambdaRef) && lambdaRef >= 0
        yMin = min(yMin, lambdaRef);
        yMax = max(yMax, lambdaRef);
    end
    
    if strcmpi(yScale, 'log')
        
        % Guard against zero or negative values for log scale.
        if yMin <= 0
            positiveVals = lambdaVals(lambdaVals > 0);
            
            if showLambdaRef && ~isempty(lambdaRef) && lambdaRef > 0
                positiveVals = [positiveVals; lambdaRef];
            end
            
            if isempty(positiveVals)
                error('No positive lambda values available for log-scale plotting.');
            end
            
            yMin = min(positiveVals, [], 'omitnan');
        end
        
        yLow  = 10^(floor(log10(yMin)) - 0.1);
        yHigh = 10^(ceil(log10(yMax)) + 0.1);
        
    else
        
        yRange = yMax - yMin;
        
        if yRange == 0
            yRange = max(abs(yMax) * 0.2, 1e-4);
        end
        
        yLow  = yMin - 0.10 * yRange;
        yHigh = yMax + 0.10 * yRange;
        
        % Lambda should not go below zero on linear scale
        yLow = max(0, yLow);
    end
    
    ylim(ax, [yLow yHigh]);
    
else
    
    ylim(ax, yLimUser);
end

%% ------------------------------------------------------------------------
% Console summaries
% -------------------------------------------------------------------------

lambdaSummaryByAnimal = groupsummary(lambdaT, mouseVar, ...
    {'mean', 'median', 'std'}, lambdaVar);

lambdaSummaryBySession = groupsummary(lambdaT, {mouseVar, sessionVar}, ...
    {'mean', 'median', 'std'}, lambdaVar);

fprintf('\nLambda summary by animal:\n');
disp(lambdaSummaryByAnimal);

fprintf('\nLambda summary by session:\n');
disp(lambdaSummaryBySession);

%% ------------------------------------------------------------------------
% Save
% -------------------------------------------------------------------------

if figSaveLogic
    
    if exist(saveDir, 'dir') ~= 7
        mkdir(saveDir);
    end
    
    saveBase = fullfile(saveDir, ...
        sprintf('%s_%s', savePrefix, datestr(now, 'mmddyy_HHMMSS')));
    
    savefig(figH, [saveBase '.fig']);
    print(figH, [saveBase '.png'], '-dpng', '-r300');
    print(figH, [saveBase '.pdf'], '-dpdf', '-bestfit', '-painters');
    
    fprintf('\nSaved legacy lambda figure:\n%s.fig\n%s.png\n%s.pdf\n', ...
        saveBase, saveBase, saveBase);
end

end