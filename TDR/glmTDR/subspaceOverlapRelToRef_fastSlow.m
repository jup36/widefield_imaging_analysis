function h = subspaceOverlapRelToRef_fastSlow(overlapC, fastIDs, slowIDs, colorMat, varargin)
%SUBSPACEOVERLAPRELTORREF_FASTSLOW Plot subspace overlap-to-reference trajectories for fast vs slow learners.
%
% h = subspaceOverlapRelToRef_fastSlow(overlapC, fastIDs, slowIDs, colorMat, 'Name', value, ...)
%
% INPUT
%   overlapC : cell array, one cell per animal
%              each cell contains a numeric vector [1 x nSess] or [nSess x 1]
%              of overlap scores relative to the reference
%
%              Example:
%              overlapRef_goToneOnC = cellfun(@(a) a.relativeToExpert.GoToneOn.overlap, ...
%                  rezSubspaceSim.perMouse, 'UniformOutput', false);
%
%   fastIDs  : logical vector, same length as overlapC
%              true for fast learners
%
%   slowIDs  : logical vector, same length as overlapC
%              true for slow learners
%
%   colorMat : [2 x 3] RGB matrix
%              row 1 = fast group color
%              row 2 = slow group color
%
% NAME-VALUE
%   'LineWidthGroup'    : line width for group mean traces (default = 2.5)
%   'LineWidthIndiv'    : line width for individual traces (default = 1)
%   'MarkerSizeGroup'   : marker size for group mean points (default = 45)
%   'IndivAlpha'        : fade factor toward white for individual lines (default = 0.75)
%   'XLabel'            : x-axis label (default = 'Session #')
%   'YLabel'            : y-axis label (default = 'Subspace overlap')
%   'TitleStr'          : title string (default = '')
%   'LegendStr'         : legend entries for fast/slow groups (default = {'Fast','Slow'})
%   'ylim'              : y-axis limits (default = [])
%   'PlotGroupPoints'   : logical, plot group mean scatter points (default = true)
%   'MakeFigure'        : logical, create a new figure (default = true)
%   'AxesHandle'        : axes handle to plot into (default = [])
%   'printFigLogic'     : logical, whether to print the figure to pdf (default = false)
%   'figSaveName'       : char/string, file name without extension (default = [])
%   'figSaveDir'        : char/string, save directory
%                         (default = 'Z:\Rodent Data\dualImaging_parkj\collectFigure\motifTDR\glmTDR\glmTDR_subSpaceSim')
%
% OUTPUT
%   h : struct
%       .fig
%       .ax
%       .lineFastMean
%       .lineSlowMean
%       .scatterFastMean
%       .scatterSlowMean
%       .lineFastIndiv
%       .lineSlowIndiv
%       .fastMat
%       .slowMat
%       .fastMean
%       .slowMean
%       .xFast
%       .xSlow

%% -------------------- parse --------------------
ip = inputParser;
ip.FunctionName = mfilename;

addRequired(ip, 'overlapC', @(x) iscell(x) && ~isempty(x));
addRequired(ip, 'fastIDs', @(x) islogical(x) || isnumeric(x));
addRequired(ip, 'slowIDs', @(x) islogical(x) || isnumeric(x));
addRequired(ip, 'colorMat', @(x) isnumeric(x) && isequal(size(x), [2 3]));

addParameter(ip, 'LineWidthGroup', 2.5, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(ip, 'LineWidthIndiv', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(ip, 'MarkerSizeGroup', 45, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(ip, 'IndivAlpha', 0.75, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
addParameter(ip, 'XLabel', 'Session #', @(x) ischar(x) || isstring(x));
addParameter(ip, 'YLabel', 'Subspace overlap', @(x) ischar(x) || isstring(x));
addParameter(ip, 'TitleStr', '', @(x) ischar(x) || isstring(x));
addParameter(ip, 'LegendStr', {'Fast','Slow'}, @(x) iscell(x) || isstring(x));
addParameter(ip, 'ylim', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2 && x(1) < x(2)));
addParameter(ip, 'axisTight', false, @(x) islogical(x) && isscalar(x));
addParameter(ip, 'PlotGroupPoints', true, @(x) islogical(x) && isscalar(x));
addParameter(ip, 'MakeFigure', true, @(x) islogical(x) && isscalar(x));
addParameter(ip, 'AxesHandle', [], @(x) isempty(x) || isgraphics(x, 'axes'));
addParameter(ip, 'printFigLogic', false, @(x) islogical(x) && isscalar(x));
addParameter(ip, 'figSaveName', [], @(x) isempty(x) || ischar(x) || isstring(x));
addParameter(ip, 'figSaveDir', 'Z:\Rodent Data\dualImaging_parkj\collectFigure\motifTDR\glmTDR\glmTDR_subSpaceSim', ...
    @(x) ischar(x) || isstring(x));

parse(ip, overlapC, fastIDs, slowIDs, colorMat, varargin{:});
P = ip.Results;

%% -------------------- validate group IDs --------------------
nAnimals = numel(overlapC);

fastIDs = logical(fastIDs(:));
slowIDs = logical(slowIDs(:));

if numel(fastIDs) ~= nAnimals || numel(slowIDs) ~= nAnimals
    error('fastIDs and slowIDs must match numel(overlapC).');
end

if ~all(xor(fastIDs, slowIDs))
    error('Each animal must belong to exactly one group (xor of fastIDs and slowIDs must be true).');
end

%% -------------------- normalize overlap vectors --------------------
% Left-align each animal to its FIRST VALID data point.
% Example:
%   [NaN NaN 0.5 0.6 0.7]  ->  [0.5 0.6 0.7]
%   [0.4 0.5 NaN 0.7]      ->  [0.4 0.5 NaN 0.7]
%
% Internal NaNs after the first valid point are preserved.

vecC = cell(nAnimals,1);
vecLens = zeros(nAnimals,1);

for i = 1:nAnimals
    y = overlapC{i};

    if isempty(y)
        vecC{i} = nan(0,1);
        vecLens(i) = 0;
        continue;
    end

    if ~isnumeric(y) || ~isvector(y)
        error('Each entry of overlapC must be a numeric vector.');
    end

    y = double(y(:)); % force column

    firstValid = find(~isnan(y), 1, 'first');

    if isempty(firstValid)
        % all-NaN case
        y = nan(0,1);
    else
        % left-align by trimming leading NaNs / missing sessions
        y = y(firstValid:end);
    end

    vecC{i} = y;
    vecLens(i) = numel(y);
end

maxLen = max(vecLens);

% Left-align into padded matrices
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

% Trim plotting extent per group to where at least one animal contributes
fastValid = any(~isnan(fastMat), 1);
slowValid = any(~isnan(slowMat), 1);

if any(fastValid)
    lastFast = find(fastValid, 1, 'last');
else
    lastFast = 0;
end

if any(slowValid)
    lastSlow = find(slowValid, 1, 'last');
else
    lastSlow = 0;
end

xFast = xAll(1:lastFast);
xSlow = xAll(1:lastSlow);

%% -------------------- set up figure/axes --------------------
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
    y = vecC{i};   % use left-aligned version
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
    y = vecC{i};   % use left-aligned version
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
    % --- tight scaling first ---
    axis(h.ax, 'tight');

    % --- then add margin ONLY to x ---
    xl = xlim(h.ax);
    xlim(h.ax, xl + [-0.5, 0.5]);

else
    % --- manual control ---
    if ~isempty(P.ylim)
        ylim(h.ax, P.ylim);
    end

    xlim(h.ax, [0.75, max([1, maxLen]) + 0.25]);
end

legendStr = cellstr(string(P.LegendStr));
if numel(legendStr) >= 2
    legend(h.ax, [h.lineFastMean, h.lineSlowMean], legendStr(1:2), ...
        'Location', 'best', 'Interpreter', 'none');
end

box(h.ax, 'off');

set(gca, "TickDir", "out")

%% -------------------- outputs --------------------
h.fastMat = fastMat;
h.slowMat = slowMat;
h.fastMean = fastMean;
h.slowMean = slowMean;
h.xFast = xFast;
h.xSlow = xSlow;

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
    fprintf('[subspaceOverlapRelToRef_fastSlow] saved: %s\n', figPath);
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