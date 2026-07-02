function figH = plot_LK_Hmetric_heatmap_perAnimal(summaryT, varargin)
% plot_LK_Hmetric_heatmap_perAnimal
%
% Plot per-animal L x K heatmaps for any scalar metric in the Lag-K summary table.
%
% This is designed for H-overlap metrics such as:
%   - mean_H_overlap_raw
%   - mean_H_overlap_norm
%   - mean_H_multi_active_fraction
%   - mean_H_active_count_mean
%   - mean_H_hoyer_sparsity
%
% Important:
%   For cross-motif H-overlap metrics, K_input = 1 should usually be excluded,
%   because there is no between-motif overlap when only one motif is present.
%
% Inputs:
%   summaryT
%       Table produced by the Lag-K parser.
%
% Name-value inputs:
%   'metricName'
%       Column name in summaryT to plot.
%
%   'metricLabel'
%       Label used in figure titles/colorbar.
%
%   'lags'
%       Lag values to plot as rows.
%
%   'Ks'
%       K_input values to plot as columns.
%
%   'KTickLabels'
%       Optional custom x tick labels.
%       Useful when using display-only K values such as 17.5 for legacy K15.
%       Example:
%           {'K=5','K=10','K=15','K=15 legacy','K=20'}
%
%   'clim'
%       Color limits. Default = [].
%
%   'figSaveLogic'
%       Whether to save figure. Default = false.
%
%   'saveDir'
%       Directory to save figure.
%
%   'saveNamePrefix'
%       Prefix for saved figure name.
%
%   'showColorbar'
%       Whether to show colorbar in each animal subplot. Default = true.
%
%   'logTransform'
%       If true, plots log10(1 + metric). Useful for raw H overlap.
%
%   'animalOrder'
%       Optional cell array/string array/categorical array specifying animal order.
%
%   'xTickAngle'
%       Rotation angle for K tick labels. Default = 25.
%
%   'nanLabel'
%       Label shown in NaN cells. Default = 'NaN'.
%
% Output:
%   figH
%       Figure handle.

%% ------------------------------------------------------------------------
% Parse inputs
% -------------------------------------------------------------------------

p = inputParser;

p.addRequired('summaryT', @istable);

p.addParameter('metricName', 'mean_H_overlap_norm', ...
    @(x) ischar(x) || isstring(x));

p.addParameter('metricLabel', '', ...
    @(x) ischar(x) || isstring(x));

p.addParameter('lags', [1 5 10], ...
    @(x) isnumeric(x) && isvector(x));

p.addParameter('Ks', [5 10], ...
    @(x) isnumeric(x) && isvector(x));

p.addParameter('KTickLabels', {}, ...
    @(x) isempty(x) || iscell(x) || isstring(x));

p.addParameter('clim', [], ...
    @(x) isempty(x) || (isnumeric(x) && numel(x) == 2));

p.addParameter('figSaveLogic', false, ...
    @(x) islogical(x) || isnumeric(x));

p.addParameter('saveDir', pwd, ...
    @(x) ischar(x) || isstring(x));

p.addParameter('saveNamePrefix', '', ...
    @(x) ischar(x) || isstring(x));

p.addParameter('showColorbar', true, ...
    @(x) islogical(x) || isnumeric(x));

p.addParameter('logTransform', false, ...
    @(x) islogical(x) || isnumeric(x));

p.addParameter('animalOrder', {}, ...
    @(x) iscell(x) || isstring(x) || iscategorical(x));

p.addParameter('xTickAngle', 25, ...
    @(x) isnumeric(x) && isscalar(x));

p.addParameter('nanLabel', 'NaN', ...
    @(x) ischar(x) || isstring(x));

p.parse(summaryT, varargin{:});

metricName      = char(p.Results.metricName);
metricLabel     = char(p.Results.metricLabel);
lags            = p.Results.lags(:)';
Ks              = p.Results.Ks(:)';
KTickLabels     = p.Results.KTickLabels;
climVals        = p.Results.clim;
figSaveLogic    = logical(p.Results.figSaveLogic);
saveDir         = char(p.Results.saveDir);
saveNamePrefix  = char(p.Results.saveNamePrefix);
showColorbar    = logical(p.Results.showColorbar);
logTransform    = logical(p.Results.logTransform);
animalOrder     = p.Results.animalOrder;
xTickAngle      = p.Results.xTickAngle;
nanLabel        = char(p.Results.nanLabel);

if isempty(metricLabel)
    metricLabel = strrep(metricName, '_', ' ');
end

if isempty(saveNamePrefix)
    saveNamePrefix = metricName;
end

%% ------------------------------------------------------------------------
% Build K tick labels
% -------------------------------------------------------------------------

if isempty(KTickLabels)
    % Use %g so display-only values like 17.5 are rendered cleanly.
    KTickLabels = arrayfun(@(k) sprintf('K=%g', k), Ks, 'UniformOutput', false);
else
    KTickLabels = cellstr(string(KTickLabels));
    
    if numel(KTickLabels) ~= numel(Ks)
        error('KTickLabels must have the same number of elements as Ks.');
    end
end

%% ------------------------------------------------------------------------
% Validate table
% -------------------------------------------------------------------------

requiredVars = {'mouseID', 'L', 'K_input', metricName};

for rv = 1:numel(requiredVars)
    if ~ismember(requiredVars{rv}, summaryT.Properties.VariableNames)
        error('Required variable "%s" not found in summaryT.', requiredVars{rv});
    end
end

% Convert mouseID to cellstr for reliable indexing
if iscategorical(summaryT.mouseID)
    mouseC = cellstr(summaryT.mouseID);
elseif isstring(summaryT.mouseID)
    mouseC = cellstr(summaryT.mouseID);
elseif iscell(summaryT.mouseID)
    mouseC = summaryT.mouseID;
else
    error('summaryT.mouseID must be categorical, string, or cell array.');
end

if isempty(animalOrder)
    animalLabels = unique(mouseC, 'stable');
else
    if iscategorical(animalOrder)
        animalLabels = cellstr(animalOrder);
    elseif isstring(animalOrder)
        animalLabels = cellstr(animalOrder);
    else
        animalLabels = animalOrder;
    end
end

nAnimals = numel(animalLabels);

if nAnimals == 0
    error('No animals found to plot.');
end

%% ------------------------------------------------------------------------
% Build figure layout
% -------------------------------------------------------------------------

nCols = ceil(sqrt(nAnimals));
nRows = ceil(nAnimals / nCols);

figH = figure('Color', 'w', ...
    'Name', sprintf('%s per animal', metricLabel), ...
    'Position', [100 100 360*nCols 310*nRows]);

tl = tiledlayout(figH, nRows, nCols, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

title(tl, sprintf('%s across L x K', metricLabel), ...
    'Interpreter', 'none', ...
    'FontWeight', 'bold');

%% ------------------------------------------------------------------------
% Plot each animal
% -------------------------------------------------------------------------

for a = 1:nAnimals
    
    ax = nexttile;
    
    animalID = animalLabels{a};
    dataMat = nan(numel(lags), numel(Ks));
    
    for li = 1:numel(lags)
        for ki = 1:numel(Ks)
            
            rowI = strcmp(mouseC, animalID) & ...
                   summaryT.L == lags(li) & ...
                   summaryT.K_input == Ks(ki);
            
            vals = summaryT.(metricName)(rowI);
            
            if ~isempty(vals)
                dataMat(li, ki) = mean(vals, 'omitnan');
            end
        end
    end
    
    if logTransform
        dataMatPlot = log10(1 + dataMat);
    else
        dataMatPlot = dataMat;
    end
    
    imagesc(ax, dataMatPlot);
    axis(ax, 'image');
    
    set(ax, ...
        'XTick', 1:numel(Ks), ...
        'XTickLabel', KTickLabels, ...
        'YTick', 1:numel(lags), ...
        'YTickLabel', compose('L=%d', lags), ...
        'TickDir', 'out', ...
        'Box', 'off', ...
        'FontSize', 10);
    
    try
        xtickangle(ax, xTickAngle);
    catch
        xtickangle(xTickAngle);
    end
    
    xlabel(ax, 'K input');
    ylabel(ax, 'Lag');
    title(ax, animalID, 'Interpreter', 'none');
    
    if ~isempty(climVals)
        caxis(ax, climVals);
    end
    
    if showColorbar
        cb = colorbar(ax);
        ylabel(cb, metricLabel, 'Interpreter', 'none');
    end
    
    %% Add numeric labels inside cells
    
    for li = 1:numel(lags)
        for ki = 1:numel(Ks)
            
            val = dataMatPlot(li, ki);
            
            if isnan(val)
                labelStr = nanLabel;
            else
                if logTransform
                    labelStr = sprintf('%.2f', val);
                elseif abs(val) >= 1000 || abs(val) < 0.01
                    labelStr = sprintf('%.2g', val);
                else
                    labelStr = sprintf('%.3f', val);
                end
            end
            
            text(ax, ki, li, labelStr, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontSize', 9, ...
                'FontWeight', 'bold', ...
                'Color', 'k', ...
                'Interpreter', 'none');
        end
    end
end

%% ------------------------------------------------------------------------
% Save figure
% -------------------------------------------------------------------------

if figSaveLogic
    
    if exist(saveDir, 'dir') ~= 7
        mkdir(saveDir);
    end
    
    logTag = '';
    if logTransform
        logTag = '_log10';
    end
    
    saveBase = sprintf('%s_LxK_heatmap_perAnimal%s_%s', ...
        saveNamePrefix, logTag, datestr(now, 'mmddyy_HHMMSS'));
    
    print(figH, fullfile(saveDir, saveBase), ...
        '-dpng', '-r300');
    
    print(figH, fullfile(saveDir, saveBase), ...
        '-dpdf', '-bestfit', '-painters');
    
    savefig(figH, fullfile(saveDir, [saveBase '.fig']));
    
    fprintf('\nSaved H-metric heatmap figure:\n%s\n', ...
        fullfile(saveDir, saveBase));
end

end