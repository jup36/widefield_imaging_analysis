function h = withinSessionOverlap_fastSlow(pairOverlapC, fastIDs, slowIDs, colorMat, varargin)
%WITHINSESSIONOVERLAP_FASTSLOW Plot within-session cross-trial-type overlap across sessions.
%
% h = withinSessionOverlap_fastSlow(pairOverlapC, fastIDs, slowIDs, colorMat, 'Name', value, ...)
%
% INPUT
%   pairOverlapC : cell array, one cell per animal
%                  each cell contains a numeric vector [nSess x 1] or [1 x nSess]
%                  for ONE cross-trial-type pair across sessions
%
%   fastIDs      : logical vector, same length as pairOverlapC
%   slowIDs      : logical vector, same length as pairOverlapC
%   colorMat     : [2 x 3] RGB matrix
%                  row 1 = fast color
%                  row 2 = slow color
%
% NAME-VALUE
%   'leftAlignFirstValid' : logical, trim leading NaNs and align first valid point to session 1 (default = true)
%   'LineWidthGroup'      : group mean line width (default = 2.5)
%   'LineWidthIndiv'      : individual line width (default = 1)
%   'MarkerSizeGroup'     : mean marker size (default = 45)
%   'IndivAlpha'          : fade factor toward white for individual lines (default = 0.75)
%   'XLabel'              : x label (default = 'Session #')
%   'YLabel'              : y label (default = 'Within-session overlap')
%   'TitleStr'            : title string (default = '')
%   'LegendStr'           : legend labels (default = {'Fast','Slow'})
%   'ylim'                : y limits (default = [])
%   'axisTight'           : use axis tight (default = false)
%   'PlotGroupPoints'     : plot mean markers (default = true)
%   'MakeFigure'          : create new figure (default = true)
%   'AxesHandle'          : axes handle to plot into (default = [])
%   'printFigLogic'       : print figure to pdf (default = false)
%   'figSaveName'         : save name without extension (default = [])
%   'figSaveDir'          : save directory
%                           (default = 'Z:\Rodent Data\dualImaging_parkj\collectFigure\motifTDR\glmTDR\glmTDR_subSpaceSim')
%
% OUTPUT
%   h : struct containing handles and summary data

%% -------------------- parse --------------------
ip = inputParser;
ip.FunctionName = mfilename;

addRequired(ip, 'pairOverlapC', @(x) iscell(x) && ~isempty(x));
addRequired(ip, 'fastIDs', @(x) islogical(x) || isnumeric(x));
addRequired(ip, 'slowIDs', @(x) islogical(x) || isnumeric(x));
addRequired(ip, 'colorMat', @(x) isnumeric(x) && isequal(size(x), [2 3]));

addParameter(ip, 'leftAlignFirstValid', true, @(x) islogical(x) && isscalar(x));
addParameter(ip, 'LineWidthGroup', 2.5, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(ip, 'LineWidthIndiv', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(ip, 'MarkerSizeGroup', 45, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(ip, 'IndivAlpha', 0.75, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
addParameter(ip, 'XLabel', 'Session #', @(x) ischar(x) || isstring(x));
addParameter(ip, 'YLabel', 'Within-session overlap', @(x) ischar(x) || isstring(x));
addParameter(ip, 'TitleStr', '', @(x) ischar(x) || isstring(x));
addParameter(ip, 'LegendStr', {'Fast','Slow'}, @(x) iscell(x) || isstring(x));
addParameter(ip, 'ylim', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2 && x(1) < x(2)));
addParameter(ip, 'axisTight', false, @(x) islogical(x) && isscalar(x));
addParameter(ip, 'PlotGroupPoints', true, @(x) islogical(x) && isscalar(x));
addParameter(ip, 'MakeFigure', true, @(x) islogical(x) && isscalar(x));
addParameter(ip, 'AxesHandle', [], @(x) isempty(x) || isgraphics(x, 'axes'));
addParameter(ip, 'printFigLogic', false, @(x) islogical(x) && isscalar(x));
addParameter(ip, 'figSaveName', [], @(x) isempty(x) || ischar(x) || isstring(x));
addParameter(ip, 'figSaveDir', ...
    'Z:\Rodent Data\dualImaging_parkj\collectFigure\motifTDR\glmTDR\glmTDR_subSpaceSim', ...
    @(x) ischar(x) || isstring(x));

parse(ip, pairOverlapC, fastIDs, slowIDs, colorMat, varargin{:});
P = ip.Results;

%% -------------------- validate group IDs --------------------
nAnimals = numel(pairOverlapC);

fastIDs = logical(fastIDs(:));
slowIDs = logical(slowIDs(:));

if numel(fastIDs) ~= nAnimals || numel(slowIDs) ~= nAnimals
    error('fastIDs and slowIDs must match numel(pairOverlapC).');
end

if ~all(xor(fastIDs, slowIDs))
    error('Each animal must belong to exactly one group.');
end

%% -------------------- normalize vectors --------------------
vecC = cell(nAnimals,1);
vecLens = zeros(nAnimals,1);

for i = 1:nAnimals
    y = pairOverlapC{i};

    if isempty(y)
        vecC{i} = nan(0,1);
        vecLens(i) = 0;
        continue;
    end

    if ~isnumeric(y) || ~isvector(y)
        error('Each entry of pairOverlapC must be a numeric vector.');
    end

    y = double(y(:)); % force column

    if P.leftAlignFirstValid
        firstValid = find(~isnan(y), 1, 'first');
        if isempty(firstValid)
            y = nan(0,1);
        else
            y = y(firstValid:end);
        end
    end

    vecC{i} = y;
    vecLens(i) = numel(y);
end

maxLen = max(vecLens);
allMat = nan(nAnimals, maxLen);

for i = 1:nAnimals
    if vecLens(i) > 0
        allMat(i,1:vecLens(i)) = vecC{i}(:)';
    end
end

fastMat = allMat(fastIDs,:);
slowMat = allMat(slowIDs,:);

fastMean = mean(fastMat, 1, 'omitnan');
slowMean = mean(slowMat, 1, 'omitnan');

xAll = 1:maxLen;
fastValid = any(~isnan(fastMat), 1);
slowValid = any(~isnan(slowMat), 1);

if any(fastValid), lastFast = find(fastValid, 1, 'last'); else, lastFast = 0; end
if any(slowValid), lastSlow = find(slowValid, 1, 'last'); else, lastSlow = 0; end

xFast = xAll(1:lastFast);
xSlow = xAll(1:lastSlow);

%% -------------------- figure/axes --------------------
if ~isempty(P.AxesHandle)
    h.ax = P.AxesHandle;
    h.fig = ancestor(h.ax, 'figure');
elseif P.MakeFigure
    h.fig = figure('Color', 'w');
    h.ax = axes('Parent', h.fig);
else
    h.ax = gca;
    h.fig = ancestor(h.ax, 'figure');
end

hold(h.ax, 'on');

%% -------------------- colors --------------------
fastColor = colorMat(1,:);
slowColor = colorMat(2,:);

fastColorFaint = blendToWhite_sub_(fastColor, P.IndivAlpha);
slowColorFaint = blendToWhite_sub_(slowColor, P.IndivAlpha);

%% -------------------- individual traces --------------------
nFast = sum(fastIDs);
nSlow = sum(slowIDs);

h.lineFastIndiv = gobjects(nFast,1);
h.lineSlowIndiv = gobjects(nSlow,1);

fastCounter = 0;
for i = find(fastIDs(:))'
    fastCounter = fastCounter + 1;
    y = vecC{i};
    x = 1:numel(y);

    if ~isempty(y)
        h.lineFastIndiv(fastCounter) = plot(h.ax, x, y, '-', ...
            'Color', fastColorFaint, ...
            'LineWidth', P.LineWidthIndiv);
    end
end

slowCounter = 0;
for i = find(slowIDs(:))'
    slowCounter = slowCounter + 1;
    y = vecC{i};
    x = 1:numel(y);

    if ~isempty(y)
        h.lineSlowIndiv(slowCounter) = plot(h.ax, x, y, '-', ...
            'Color', slowColorFaint, ...
            'LineWidth', P.LineWidthIndiv);
    end
end

%% -------------------- group mean traces --------------------
if ~isempty(xFast)
    h.lineFastMean = plot(h.ax, xFast, fastMean(1:lastFast), '-', ...
        'Color', fastColor, ...
        'LineWidth', P.LineWidthGroup);
else
    h.lineFastMean = gobjects(1);
end

if ~isempty(xSlow)
    h.lineSlowMean = plot(h.ax, xSlow, slowMean(1:lastSlow), '-', ...
        'Color', slowColor, ...
        'LineWidth', P.LineWidthGroup);
else
    h.lineSlowMean = gobjects(1);
end

if P.PlotGroupPoints
    if ~isempty(xFast)
        h.scatterFastMean = scatter(h.ax, xFast, fastMean(1:lastFast), ...
            P.MarkerSizeGroup, fastColor, 'filled');
    else
        h.scatterFastMean = gobjects(1);
    end

    if ~isempty(xSlow)
        h.scatterSlowMean = scatter(h.ax, xSlow, slowMean(1:lastSlow), ...
            P.MarkerSizeGroup, slowColor, 'filled');
    else
        h.scatterSlowMean = gobjects(1);
    end
else
    h.scatterFastMean = gobjects(1);
    h.scatterSlowMean = gobjects(1);
end

%% -------------------- formatting --------------------
xlabel(h.ax, P.XLabel, 'Interpreter', 'none');
ylabel(h.ax, P.YLabel, 'Interpreter', 'none');

if strlength(string(P.TitleStr)) > 0
    title(h.ax, P.TitleStr, 'Interpreter', 'none');
end

if P.axisTight
    axis(h.ax, 'tight');

    % --- add x margin AFTER axis tight ---
    xl = xlim(h.ax);
    xlim(h.ax, [xl(1)-0.5, xl(2)+0.5]);

else
    if ~isempty(P.ylim)
        ylim(h.ax, P.ylim);
    end

    % existing behavior (already has small margin)
    xlim(h.ax, [0.75, max([1, maxLen]) + 0.25]);
end

legendStr = cellstr(string(P.LegendStr));
if numel(legendStr) >= 2
    legend(h.ax, [h.lineFastMean, h.lineSlowMean], legendStr(1:2), ...
        'Location', 'best', 'Interpreter', 'none');
end

box(h.ax, 'off');
set(gca, 'TickDir', 'out')

%% -------------------- outputs --------------------
h.fastMat = fastMat;
h.slowMat = slowMat;
h.fastMean = fastMean;
h.slowMean = slowMean;
h.xFast = xFast;
h.xSlow = xSlow;
h.vecC = vecC;

%% -------------------- optional print --------------------
if P.printFigLogic
    if isempty(P.figSaveName) || strlength(string(P.figSaveName)) == 0
        error('When ''printFigLogic'' is true, provide a non-empty ''figSaveName''.');
    end

    figSaveDir = char(string(P.figSaveDir));
    figSaveName = char(string(P.figSaveName));

    if ~exist(figSaveDir, 'dir')
        mkdir(figSaveDir);
    end

    figPath = fullfile(figSaveDir, [sanitizeFileName_sub_(figSaveName) '.pdf']);
    set(h.fig, 'PaperPositionMode', 'auto');
    print(h.fig, figPath, '-dpdf', '-bestfit', '-painters');
    fprintf('[withinSessionOverlap_fastSlow] saved: %s\n', figPath);
end

hold(h.ax, 'off');

end

%% ========================= local helpers =========================
function cOut = blendToWhite_sub_(cIn, fracToWhite)
cIn = cIn(:)';
fracToWhite = max(min(fracToWhite,1),0);
cOut = (1-fracToWhite)*cIn + fracToWhite*[1 1 1];
cOut = max(min(cOut,1),0);
end

function fn = sanitizeFileName_sub_(fn)
fn = char(fn);
bad = '<>:"/\|?*';
for k = 1:numel(bad)
    fn(fn==bad(k)) = '_';
end
fn = regexprep(fn, '\s+', '_');
end