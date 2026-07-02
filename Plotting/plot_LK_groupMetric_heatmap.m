function figH = plot_LK_groupMetric_heatmap(groupLK_T, varargin)
% plot_LK_groupMetric_heatmap
%
% Plot across-animal group summary as L x K_input heatmap.
%
% Each tile shows:
%   group mean
%   +/- SEM
%   n animals
%
% Input:
%   groupLK_T from summarize_LK_acrossAnimals
%
% Optional name-value pairs:
%   'metricName'   : original metric name, e.g. 'meanPEV'
%   'metricLabel'  : display label
%   'lags'         : lags to plot
%   'Ks'           : K_input values to plot
%   'KTickLabels'  : custom K-axis labels.
%                    Useful when using display-only K values such as 17.5
%                    to represent legacy K15.
%                    Example:
%                        {'K=1','K=5','K=10','K=15','K=15 legacy','K=20'}
%   'clim'         : color limits
%   'valueFormat'  : numeric format for mean
%   'semFormat'    : numeric format for SEM
%   'figSaveLogic' : true/false
%   'saveDir'      : save directory
%   'savePrefix'   : save filename prefix
%   'showColorbar' : true/false. Default false due to MATLAB colorbar instability.
%   'xTickAngle'   : x tick label angle. Default = 25.
%
% Output:
%   figH : figure handle

%% ------------------------------------------------------------------------
% Parse inputs
% -------------------------------------------------------------------------

p = inputParser;

addParameter(p, 'metricName', 'meanPEV', ...
    @(x) ischar(x) || isstring(x));

addParameter(p, 'metricLabel', '', ...
    @(x) ischar(x) || isstring(x));

addParameter(p, 'lags', [], ...
    @(x) isempty(x) || isnumeric(x));

addParameter(p, 'Ks', [], ...
    @(x) isempty(x) || isnumeric(x));

addParameter(p, 'KTickLabels', {}, ...
    @(x) isempty(x) || iscell(x) || isstring(x));

addParameter(p, 'clim', [], ...
    @(x) isempty(x) || isnumeric(x));

addParameter(p, 'valueFormat', '%.3f', ...
    @(x) ischar(x) || isstring(x));

addParameter(p, 'semFormat', '%.3f', ...
    @(x) ischar(x) || isstring(x));

addParameter(p, 'figSaveLogic', false, ...
    @(x) islogical(x) || isnumeric(x));

addParameter(p, 'saveDir', pwd, ...
    @(x) ischar(x) || isstring(x));

addParameter(p, 'savePrefix', '', ...
    @(x) ischar(x) || isstring(x));

addParameter(p, 'showColorbar', false, ...
    @(x) islogical(x) || isnumeric(x));

addParameter(p, 'xTickAngle', 25, ...
    @(x) isnumeric(x) && isscalar(x));

parse(p, varargin{:});

metricName    = char(p.Results.metricName);
metricLabel   = char(p.Results.metricLabel);
lags          = p.Results.lags;
Ks            = p.Results.Ks;
KTickLabels   = p.Results.KTickLabels;
climVals      = p.Results.clim;
valueFormat   = char(p.Results.valueFormat);
semFormat     = char(p.Results.semFormat);
figSaveLogic  = logical(p.Results.figSaveLogic);
saveDir       = char(p.Results.saveDir);
savePrefix    = char(p.Results.savePrefix);
showColorbar  = logical(p.Results.showColorbar);
xTickAngle    = p.Results.xTickAngle;

if isempty(metricLabel)
    metricLabel = strrep(metricName, '_', ' ');
end

if isempty(savePrefix)
    savePrefix = ['group_LK_heatmap_' metricName];
end

meanVar = ['groupMean_' metricName];
semVar  = ['groupSem_'  metricName];

%% ------------------------------------------------------------------------
% Validate table
% -------------------------------------------------------------------------

requiredVars = {'L', 'K_input', 'nAnimals', meanVar, semVar};
missingVars = setdiff(requiredVars, groupLK_T.Properties.VariableNames);

if ~isempty(missingVars)
    error('groupLK_T is missing required variable(s): %s', strjoin(missingVars, ', '));
end

if isempty(lags)
    lags = sort(unique(groupLK_T.L(:)))';
end

if isempty(Ks)
    Ks = sort(unique(groupLK_T.K_input(:)))';
end

lags = lags(:)';
Ks   = Ks(:)';

%% ------------------------------------------------------------------------
% Build K tick labels
% -------------------------------------------------------------------------

if isempty(KTickLabels)
    % Use %g so display-only values like 17.5 are shown cleanly.
    KTickLabels = arrayfun(@(k) sprintf('K=%g', k), Ks, 'UniformOutput', false);
else
    KTickLabels = cellstr(string(KTickLabels));

    if numel(KTickLabels) ~= numel(Ks)
        error('KTickLabels must have the same number of elements as Ks.');
    end
end

%% ------------------------------------------------------------------------
% Build matrices
% -------------------------------------------------------------------------

meanMat = nan(numel(lags), numel(Ks));
semMat  = nan(numel(lags), numel(Ks));
nMat    = nan(numel(lags), numel(Ks));

for iL = 1:numel(lags)
    for iK = 1:numel(Ks)

        rowI = groupLK_T.L == lags(iL) & groupLK_T.K_input == Ks(iK);

        if any(rowI)
            % If somehow there are duplicate rows, average them defensively.
            meanMat(iL, iK) = mean(groupLK_T.(meanVar)(rowI), 'omitnan');
            semMat(iL, iK)  = mean(groupLK_T.(semVar)(rowI),  'omitnan');
            nMat(iL, iK)    = max(groupLK_T.nAnimals(rowI), [], 'omitnan');
        end
    end
end

%% ------------------------------------------------------------------------
% Color limits
% -------------------------------------------------------------------------

if isempty(climVals)
    climVals = [min(meanMat(:), [], 'omitnan'), max(meanMat(:), [], 'omitnan')];

    % Protect against all-NaN or flat data
    if any(isnan(climVals))
        climVals = [0 1];
    elseif climVals(1) == climVals(2)
        padVal = max(abs(climVals(1)) * 0.05, 0.01);
        climVals = climVals + [-padVal padVal];
    end
end

if numel(climVals) ~= 2 || climVals(1) >= climVals(2)
    error('clim must be a two-element increasing numeric vector.');
end

%% ------------------------------------------------------------------------
% Plot
% -------------------------------------------------------------------------

figH = figure('Color', 'w', 'Position', [300 300 620 460]);

ax = axes(figH);

imagesc(ax, meanMat);
axis(ax, 'image');

colormap(ax, parula);
caxis(ax, climVals);

xticks(ax, 1:numel(Ks));
xticklabels(ax, KTickLabels);

try
    xtickangle(ax, xTickAngle);
catch
    xtickangle(xTickAngle);
end

xlabel(ax, 'K input');

yticks(ax, 1:numel(lags));
yticklabels(ax, arrayfun(@(x) sprintf('L=%g', x), lags, 'UniformOutput', false));
ylabel(ax, 'Lag');

title(ax, sprintf('Across-animal %s\nColor scale %.3g to %.3g; text = mean \\pm SEM across animals', ...
    metricLabel, climVals(1), climVals(2)), ...
    'Interpreter', 'none', ...
    'FontWeight', 'bold');

hold(ax, 'on');

for iL = 1:numel(lags)
    for iK = 1:numel(Ks)

        val    = meanMat(iL, iK);
        semVal = semMat(iL, iK);
        nVal   = nMat(iL, iK);

        if isnan(val)
            txt = 'n/a';
            txtColor = [0.5 0.5 0.5];
        else
            if isnan(semVal)
                semVal = NaN;
            end

            if isnan(nVal)
                txt = sprintf([valueFormat '\n\\pm ' semFormat '\nn/a'], ...
                    val, semVal);
            else
                txt = sprintf([valueFormat '\n\\pm ' semFormat '\nn=%d'], ...
                    val, semVal, round(nVal));
            end

            normVal = (val - climVals(1)) / (climVals(2) - climVals(1));

            if normVal > 0.55
                txtColor = 'k';
            else
                txtColor = 'w';
            end
        end

        text(ax, iK, iL, txt, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontWeight', 'bold', ...
            'FontSize', 10, ...
            'Color', txtColor, ...
            'Interpreter', 'tex');
    end
end

hold(ax, 'off');

%% ------------------------------------------------------------------------
% Optional colorbar
% -------------------------------------------------------------------------

if showColorbar
    try
        cb = colorbar(ax);
        set(cb, 'Units', 'normalized');
        % Avoid cb.Label.String / ylabel(cb,...) because your MATLAB had
        % colorbar listener issues.
    catch ME
        warning('Colorbar failed: %s. Continuing without colorbar.', ME.message);
    end
end

%% ------------------------------------------------------------------------
% Save
% -------------------------------------------------------------------------

if figSaveLogic

    if exist(saveDir, 'dir') ~= 7
        mkdir(saveDir);
    end

    saveNameBase = fullfile(saveDir, ...
        sprintf('%s_%s', savePrefix, datestr(now, 'mmddyy_HHMMSS')));

    savefig(figH, [saveNameBase '.fig']);
    print(figH, [saveNameBase '.png'], '-dpng', '-r300');

    fprintf('Saved group heatmap:\n%s.fig\n%s.png\n', saveNameBase, saveNameBase);
end

end