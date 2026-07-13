function h = plotMeanSemShadeOverTime(meanMat, semMat, x, colorMat, varargin)
%PLOTMEANSEMSHADEOVERTIME Plot group mean +/- SEM over time with shaded error.
%
% h = plotMeanSemShadeOverTime(meanMat, semMat, x, colorMat, 'Name', value, ...)
%
% INPUT
%   meanMat  : [nGroup x nTime] matrix of mean values
%   semMat   : [nGroup x nTime] matrix of SEM values
%   x        : [1 x nTime] or [nTime x 1] time vector
%   colorMat : [nGroup x 3] RGB color matrix
%
% NAME-VALUE PAIRS
%   'LineWidth'    : line width for mean traces (default = 2)
%   'FaceAlpha'    : alpha for shaded SEM patches (default = 0.25)
%   'SmoothN'      : smoothing factor passed to smooth2a(..., 0, SmoothN)
%                    default = 0
%   'XLabel'       : x-axis label (default = 'Time')
%   'YLabel'       : y-axis label (default = 'Value')
%   'TitleStr'     : title string (default = '')
%   'LegendStr'    : cellstr/string of legend entries (default = {})
%   'PlotZeroLine' : logical scalar (default = false)
%   'ZeroLineX'    : x location for vertical line, default = 0
%   'ylim'         : optional y-axis limits, e.g. [0 1].
%                    default = []
%
%   'figSaveLogic' : logical scalar; whether to save the figure (default = false)
%   'figSaveDir'   : figure save directory (default = [])
%   'figName'      : filename stem for saved figure, without extension
%                    (default = [])
%
% OUTPUT
%   h : struct with handles
%       .fig
%       .ax
%       .patch
%       .line
%       .figPath
%

%% Parse inputs
ip = inputParser;
ip.FunctionName = mfilename;

addRequired(ip, 'meanMat', @(x) isnumeric(x) && ismatrix(x) && ~isempty(x));
addRequired(ip, 'semMat',  @(x) isnumeric(x) && ismatrix(x) && ~isempty(x));
addRequired(ip, 'x',       @(x) isnumeric(x) && isvector(x));
addRequired(ip, 'colorMat', @(x) isnumeric(x) && size(x,2)==3);

addParameter(ip, 'LineWidth', 2, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(ip, 'FaceAlpha', 0.25, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
addParameter(ip, 'SmoothN', 0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(ip, 'XLabel', 'Time', @(x) ischar(x) || isstring(x));
addParameter(ip, 'YLabel', 'Value', @(x) ischar(x) || isstring(x));
addParameter(ip, 'TitleStr', '', @(x) ischar(x) || isstring(x));
addParameter(ip, 'LegendStr', {}, @(x) iscell(x) || isstring(x));
addParameter(ip, 'PlotZeroLine', false, @(x) islogical(x) && isscalar(x));
addParameter(ip, 'ZeroLineX', 0, @(x) isnumeric(x) && isscalar(x));

% New optional y-axis limits
addParameter(ip, 'ylim', [], ...
    @(x) isempty(x) || (isnumeric(x) && numel(x)==2 && x(1) < x(2)));

addParameter(ip, 'figSaveLogic', false, @(x) islogical(x) && isscalar(x));
addParameter(ip, 'figSaveDir', [], @(x) isempty(x) || ischar(x) || isstring(x));
addParameter(ip, 'figName', [], @(x) isempty(x) || ischar(x) || isstring(x));

parse(ip, meanMat, semMat, x, colorMat, varargin{:});
P = ip.Results;

%% Validate dimensions
x = x(:)'; % force row vector

if ~isequal(size(meanMat), size(semMat))
    error('meanMat and semMat must be the same size.');
end

[nGroup, nTime] = size(meanMat);

if numel(x) ~= nTime
    error('Length of x must match number of columns in meanMat/semMat.');
end

if size(colorMat,1) ~= nGroup
    error('colorMat must have one RGB row per group.');
end

%% Create figure/axes
h.fig = figure('Color', 'w');
h.ax = axes('Parent', h.fig);
hold(h.ax, 'on');

h.patch = gobjects(nGroup,1);
h.line  = gobjects(nGroup,1);
h.figPath = [];

%% Plot each group
for i = 1:nGroup
    mu  = double(meanMat(i,:));
    sem = double(semMat(i,:));
    col = colorMat(i,:);

    if P.SmoothN > 0
        mu  = smooth2a(mu,  0, P.SmoothN);
        sem = smooth2a(sem, 0, P.SmoothN);
    end

    upper = mu + sem;
    lower = mu - sem;

    h.patch(i) = fill(h.ax, ...
        [x fliplr(x)], ...
        [upper fliplr(lower)], ...
        col, ...
        'FaceAlpha', P.FaceAlpha, ...
        'EdgeColor', 'none');

    h.line(i) = plot(h.ax, x, mu, ...
        'Color', col, ...
        'LineWidth', P.LineWidth);
end

if P.PlotZeroLine
    xline(h.ax, P.ZeroLineX, '--k', 'LineWidth', 1);
end

xlabel(h.ax, P.XLabel, 'Interpreter', 'none');
ylabel(h.ax, P.YLabel, 'Interpreter', 'none');

if strlength(string(P.TitleStr)) > 0
    title(h.ax, P.TitleStr, 'Interpreter', 'none');
end

if ~isempty(P.LegendStr)
    legend(h.ax, h.line, cellstr(string(P.LegendStr)), ...
        'Location', 'best', 'Interpreter', 'none');
end

% Apply optional y-axis limits
if ~isempty(P.ylim)
    ylim(h.ax, P.ylim);
end

box(h.ax, 'off');

%% Save figure if requested
if P.figSaveLogic
    if ~isempty(P.figSaveDir) && ~isempty(P.figName)
        figSaveDir = char(string(P.figSaveDir));
        figName    = char(string(P.figName));

        if ~exist(figSaveDir, 'dir')
            mkdir(figSaveDir);
        end

        figName = sanitize_filename_local(figName);
        figPath = fullfile(figSaveDir, [figName '.pdf']);

        set(h.fig, 'PaperPositionMode', 'auto');
        print(h.fig, figPath, '-dpdf', '-painters', '-bestfit');

        fprintf('[plotMeanSemShadeOverTime] saved: %s\n', figPath);
        h.figPath = figPath;
    end
end

end

%% ---------- Local helper ----------
function fn = sanitize_filename_local(fn)
fn = char(fn);
bad = '<>:"/\|?*';
for k = 1:numel(bad)
    fn(fn==bad(k)) = '_';
end
fn = regexprep(fn, '\s+', '_');
end