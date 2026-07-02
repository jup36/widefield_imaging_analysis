function figH = plot_LK_meanK_heatmap_perAnimal(pevSummaryT, varargin)
% plot_LK_meanK_heatmap_perAnimal
%
% Visualize mean discovered K across L-K_input combinations per animal.
%
% Input:
%   pevSummaryT : table with variables:
%       mouseID, sessionID, L, K_input, meanK_discovered
%
% Optional name-value pairs:
%   'lags'        : vector of lags to plot. Default = sorted unique L.
%   'Ks'          : vector of K_input values to plot. Default = sorted unique K_input.
%   'KTickLabels' : custom K-axis labels. Useful for display-only values
%                   such as K=17.5 representing legacy K15.
%                   Example:
%                       {'K=1','K=5','K=10','K=15','K=15 legacy','K=20'}
%   'clim'        : color limits. Default = [].
%   'figSaveLogic': true/false. Default = false.
%   'saveDir'     : directory to save figures. Default = pwd.
%   'savePrefix'  : filename prefix. Default = 'meanK_LK_heatmap'
%   'valueFormat' : text format for tile values. Default = '%.2f'
%   'showColorbar': true/false. Default = false.
%   'xTickAngle'  : x tick label angle. Default = 25.
%
% Output:
%   figH : figure handle

%% Parse inputs
p = inputParser;

addParameter(p, 'lags', [], @(x) isnumeric(x));
addParameter(p, 'Ks', [], @(x) isnumeric(x));
addParameter(p, 'KTickLabels', {}, @(x) isempty(x) || iscell(x) || isstring(x));
addParameter(p, 'clim', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'figSaveLogic', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'saveDir', pwd, @(x) ischar(x) || isstring(x));
addParameter(p, 'savePrefix', 'meanK_LK_heatmap', @(x) ischar(x) || isstring(x));
addParameter(p, 'valueFormat', '%.2f', @(x) ischar(x) || isstring(x));
addParameter(p, 'showColorbar', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'xTickAngle', 25, @(x) isnumeric(x) && isscalar(x));

parse(p, varargin{:});

lags = p.Results.lags;
Ks = p.Results.Ks;
KTickLabels = p.Results.KTickLabels;
climVals = p.Results.clim;
figSaveLogic = logical(p.Results.figSaveLogic);
saveDir = char(p.Results.saveDir);
savePrefix = char(p.Results.savePrefix);
valueFormat = char(p.Results.valueFormat);
showColorbar = logical(p.Results.showColorbar);
xTickAngle = p.Results.xTickAngle;

%% Basic checks
requiredVars = {'mouseID', 'L', 'K_input', 'meanK_discovered'};
missingVars = setdiff(requiredVars, pevSummaryT.Properties.VariableNames);

if ~isempty(missingVars)
    error('pevSummaryT is missing required variable(s): %s', strjoin(missingVars, ', '));
end

if isempty(lags)
    lags = sort(unique(pevSummaryT.L(:)))';
end

if isempty(Ks)
    Ks = sort(unique(pevSummaryT.K_input(:)))';
end

lags = lags(:)';
Ks = Ks(:)';

%% Build K tick labels
if isempty(KTickLabels)
    % Use %g so display-only values like 17.5 are shown cleanly.
    KTickLabels = arrayfun(@(k) sprintf('K=%g', k), Ks, 'UniformOutput', false);
else
    KTickLabels = cellstr(string(KTickLabels));

    if numel(KTickLabels) ~= numel(Ks)
        error('KTickLabels must have the same number of elements as Ks.');
    end
end

mouseList = unique(pevSummaryT.mouseID, 'stable');

if iscategorical(mouseList)
    mouseListStr = cellstr(mouseList);
else
    mouseListStr = cellstr(string(mouseList));
end

nMice = numel(mouseList);

if nMice == 0
    error('No mice found in pevSummaryT.');
end

%% Global color limits if not provided
if isempty(climVals)
    allVals = pevSummaryT.meanK_discovered;
    climVals = [min(allVals, [], 'omitnan'), max(allVals, [], 'omitnan')];

    % Protect against flat or all-NaN data
    if any(isnan(climVals))
        climVals = [0 1];
    elseif climVals(1) == climVals(2)
        climVals = climVals + [-0.5 0.5];
    end
end

%% Decide subplot layout
nCols = ceil(sqrt(nMice));
nRows = ceil(nMice / nCols);

figH = figure('Color', 'w', 'Position', [100 100 360*nCols 310*nRows]);

tiledlayout(nRows, nCols, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

%% Plot each animal
axC = gobjects(nMice, 1);

for m = 1:nMice

    thisMouse = mouseList(m);

    if iscategorical(pevSummaryT.mouseID)
        thisMouseRows = pevSummaryT.mouseID == thisMouse;
    else
        thisMouseRows = strcmp(string(pevSummaryT.mouseID), string(thisMouse));
    end

    thisT = pevSummaryT(thisMouseRows, :);

    % Matrix: rows = L, columns = K_input
    meanKMat = nan(numel(lags), numel(Ks));

    for iL = 1:numel(lags)
        for iK = 1:numel(Ks)

            rowI = thisT.L == lags(iL) & thisT.K_input == Ks(iK);

            if any(rowI)
                % Average across sessions for this animal and this L/K pair
                meanKMat(iL, iK) = mean(thisT.meanK_discovered(rowI), 'omitnan');
            end
        end
    end

    ax = nexttile;
    axC(m) = ax;

    imagesc(ax, meanKMat);
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

    title(ax, mouseListStr{m}, 'Interpreter', 'none');

    % Overlay values
    hold(ax, 'on');

    for iL = 1:numel(lags)
        for iK = 1:numel(Ks)

            val = meanKMat(iL, iK);

            if isnan(val)
                txt = 'n/a';
                txtColor = [0.5 0.5 0.5];
            else
                txt = sprintf(valueFormat, val);

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
                'FontSize', 11, ...
                'Color', txtColor, ...
                'Interpreter', 'none');
        end
    end

    hold(ax, 'off');
end

%% Apply common color limits
for aa = 1:numel(axC)
    if isgraphics(axC(aa))
        caxis(axC(aa), climVals);
    end
end

%% Optional colorbar
if showColorbar
    try
        validAx = axC(isgraphics(axC));

        if ~isempty(validAx)
            cb = colorbar(validAx(end));
            set(cb, 'Units', 'normalized');
            set(cb, 'Position', [0.93 0.15 0.015 0.70]);
        end

        sgtitle('Mean discovered K across lag-K combinations per animal', ...
            'FontWeight', 'bold');

    catch ME
        warning('Colorbar failed: %s. Continuing without colorbar.', ME.message);

        sgtitle(sprintf('Mean discovered K across lag-K combinations per animal; color scale %.2f to %.2f', ...
            climVals(1), climVals(2)), ...
            'FontWeight', 'bold');
    end
else
    sgtitle(sprintf('Mean discovered K across lag-K combinations per animal; color scale %.2f to %.2f', ...
        climVals(1), climVals(2)), ...
        'FontWeight', 'bold');
end

%% Save
if figSaveLogic

    if exist(saveDir, 'dir') ~= 7
        mkdir(saveDir);
    end

    saveNameBase = fullfile(saveDir, ...
        sprintf('%s_%s', savePrefix, datestr(now, 'mmddyy_HHMMSS')));

    savefig(figH, [saveNameBase '.fig']);
    print(figH, [saveNameBase '.png'], '-dpng', '-r300');

    fprintf('Saved figure:\n%s.fig\n%s.png\n', saveNameBase, saveNameBase);
end

end