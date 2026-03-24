function h = plotCellImagescGrid(dataC, varargin)
%PLOTCELLIMAGESCGRID Plot one imagesc map for each cell entry in a cell array.
%
% h = plotCellImagescGrid(dataC, 'Name', value, ...)
%
% INPUT
%   dataC : cell array, each cell contains a 2D numeric matrix
%           Example: 9x1 cell, each cell is [19 x 117]
%
% NAME-VALUE PAIRS
%   'timeX'      : numeric vector for x-axis values (default = [])
%                  If provided, length must match number of columns.
%                  Typically this is time relative to tone onset, where 0 = tone onset.
%
%   'colormap'   : colormap matrix to apply to all plots
%                  (default = parula(256))
%
%   'axisI'      : logical or numeric index vector selecting rows of each matrix
%                  to plot (default = [])
%
%   'header'     : cell array of labels, one per cell entry
%                  Used as y-axis label for each subplot (default = [])
%
%   'clims'      : 1x2 numeric vector for color limits, e.g. [-1.5 1.5]
%                  (default = [])
%
%   'titleStr'   : overall figure title (default = '')
%
%   'layout'     : [nRows nCols] for subplot layout
%                  (default = auto)
%
%   'showColorbar' : logical scalar, whether to show colorbar on each subplot
%                    (default = false)
%
%   'xlabelStr'  : x-axis label (default = 'Time')
%
%   'ylabelStr'  : fallback y-axis label when header is not provided
%                  (default = 'Axis')
%
%   'figureHandle' : existing figure handle to plot into (default = [])
%
% OUTPUT
%   h : struct containing handles
%       .fig  - figure handle
%       .ax   - axes handles
%       .im   - image handles
%
% EXAMPLE
%   plotCellImagescGrid(muNoGo_slopeC, ...
%       'timeX', -2:0.05:3.8, ...
%       'colormap', turbo(256), ...
%       'axisI', true(19,1), ...
%       'header', {'M1','M2','M3','M4','M5','M6','M7','M8','M9'}, ...
%       'clims', [-1.5 1.5]);
%
% NOTE
%   This function assumes each selected matrix has the same number of columns
%   if 'timeX' is provided.

%% Parse inputs
ip = inputParser;
ip.FunctionName = mfilename;

addRequired(ip, 'dataC', @(x) iscell(x) && ~isempty(x));

addParameter(ip, 'timeX', [], @(x) isempty(x) || isnumeric(x));
addParameter(ip, 'colormap', parula(256), @(x) isnumeric(x) && size(x,2)==3);
addParameter(ip, 'axisI', [], @(x) isempty(x) || isnumeric(x) || islogical(x));
addParameter(ip, 'header', [], @(x) isempty(x) || iscellstr(x) || isstring(x) || iscell(x));
addParameter(ip, 'clims', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
addParameter(ip, 'titleStr', '', @(x) ischar(x) || isstring(x));
addParameter(ip, 'layout', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
addParameter(ip, 'showColorbar', false, @(x) islogical(x) && isscalar(x));
addParameter(ip, 'xlabelStr', 'Time', @(x) ischar(x) || isstring(x));
addParameter(ip, 'ylabelStr', 'Axis', @(x) ischar(x) || isstring(x));
addParameter(ip, 'figureHandle', [], @(x) isempty(x) || isgraphics(x, 'figure'));

parse(ip, dataC, varargin{:});
P = ip.Results;

nPlots = numel(dataC);

%% Validate cell contents
for i = 1:nPlots
    if ~isnumeric(dataC{i}) || ndims(dataC{i}) ~= 2
        error('Each cell entry in dataC must contain a 2D numeric matrix.');
    end
end

%% Determine layout
if isempty(P.layout)
    nCols = ceil(sqrt(nPlots));
    nRows = ceil(nPlots / nCols);
else
    nRows = P.layout(1);
    nCols = P.layout(2);
    if nRows * nCols < nPlots
        error('Specified layout [%d %d] cannot accommodate %d plots.', nRows, nCols, nPlots);
    end
end

%% Create / reuse figure
if isempty(P.figureHandle)
    h.fig = figure('Color', 'w');
else
    h.fig = P.figureHandle;
    figure(h.fig);
end

h.ax = gobjects(nPlots,1);
h.im = gobjects(nPlots,1);

%% Plot each cell entry
for i = 1:nPlots
    M = dataC{i};

    % Apply row selection if requested
    if ~isempty(P.axisI)
        try
            M = M(P.axisI, :);
        catch ME
            error('Failed to index rows for cell %d using axisI.\n%s', i, ME.message);
        end
    end

    % Validate timeX if provided
    if ~isempty(P.timeX)
        if size(M,2) ~= numel(P.timeX)
            error('For cell %d, size(M,2) = %d but numel(timeX) = %d.', ...
                i, size(M,2), numel(P.timeX));
        end
    end

    h.ax(i) = subplot(nRows, nCols, i);

    if isempty(P.timeX)
        h.im(i) = imagesc(M);
    else
        h.im(i) = imagesc(P.timeX, 1:size(M,1), M);
    end

    axis(h.ax(i), 'tight');
    set(h.ax(i), 'YDir', 'normal');

    % Apply colormap and clims
    colormap(h.ax(i), P.colormap);
    if ~isempty(P.clims)
        clim(h.ax(i), P.clims);
    end

    % Labels
    xlabel(h.ax(i), P.xlabelStr);

    if ~isempty(P.header)
        if numel(P.header) >= i
            ylabel(h.ax(i), string(P.header{i}), 'Interpreter', 'none');
        else
            ylabel(h.ax(i), P.ylabelStr);
        end
    else
        ylabel(h.ax(i), P.ylabelStr);
    end

    % Mark tone onset if timeX contains 0
    if ~isempty(P.timeX)
        hold(h.ax(i), 'on');
        if any(P.timeX == 0)
            xline(h.ax(i), 0, '--k', 'LineWidth', 1);
        end
        hold(h.ax(i), 'off');
    end

    % Optional colorbar
    if P.showColorbar
        colorbar(h.ax(i));
    end

    title(h.ax(i), sprintf('Cell %d', i));
end

% Overall title
if strlength(string(P.titleStr)) > 0
    sgtitle(P.titleStr);
end

end