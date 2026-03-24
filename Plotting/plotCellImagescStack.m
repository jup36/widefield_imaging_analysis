function [h, dataC_out] = plotCellImagescStack(dataC, varargin)
%PLOTCELLIMAGESCSTACK Stack matrices from a cell array into one imagesc plot.
%
% [h, dataC_out] = plotCellImagescStack(dataC, 'Name', value, ...)
%
% INPUT
%   dataC : cell array, each cell contains a 2D numeric matrix
%           Example: 9x1 cell, each cell is [19 x 117]
%
% NAME-VALUE PAIRS
%   'timeX'        : numeric vector for x-axis values (default = [])
%                    If provided, length must match number of columns.
%
%   'colormap'     : colormap matrix to apply (default = parula(256))
%
%   'axisI'        : logical or numeric index vector selecting rows from each matrix
%                    (default = [])
%
%   'header'       : cell array of labels, one per cell entry
%                    Used as y tick labels at the vertical center of each mouse block
%                    (default = [])
%
%   'clims'        : 1x2 numeric vector for color limits, e.g. [-1.5 1.5]
%                    (default = [])
%
%   'titleStr'     : figure title (default = '')
%
%   'xlabelStr'    : x-axis label (default = 'Time')
%
%   'ylabelStr'    : y-axis label (default = '')
%
%   'showColorbar' : logical scalar, whether to show colorbar (default = true)
%
%   'sigMask'      : cell array same size as dataC; each cell contains a logical mask
%                    same size as the corresponding matrix in dataC.
%                    Entries where sigMask is false are set to zero before plotting.
%                    (default = [])
%
%   'absSlope'     : logical scalar, whether to take absolute value of each
%                    processed matrix AFTER sigMask and axisI are applied.
%                    If true, both visualization and dataC_out use abs(M).
%                    (default = false)
%
%   'boundaryColor': color for mouse block boundaries (default = 'k')
%
%   'boundaryWidth': line width for mouse block boundaries (default = 1.5)
%
%   'figSaveLogic' : logical scalar, whether to save the figure as PDF
%                    (default = false)
%
%   'figSaveDir'   : figure save directory
%                    (default = 'Z:\Rodent Data\dualImaging_parkj\collectFigure\motifTDR\glmTDR\glmTDR_across_session_slope')
%
%   'figName'      : filename stem for saving (default = [])
%                    Example: "muNoGo_slope_projNoGoAxis"
%
% OUTPUT
%   h : struct containing handles
%       .fig           - figure handle
%       .ax            - axes handle
%       .im            - image handle
%       .stackedMat    - stacked matrix actually plotted
%       .blockStart    - starting row of each mouse block in stacked matrix
%       .blockEnd      - ending row of each mouse block in stacked matrix
%       .blockCenter   - vertical center of each mouse block
%       .figPath       - saved PDF path (if saved), else []
%
%   dataC_out : cell array, same size as dataC
%       Each cell contains the processed matrix after applying:
%       1) sigMask (if provided)
%       2) axisI row selection (if provided)
%       3) abs()      (if absSlope = true)
%
% EXAMPLE
%   [h, dataC_out] = plotCellImagescStack(muNoGo_slopeC, ...
%       'timeX', timeX, ...
%       'colormap', rb, ...
%       'axisI', (goToneOnI | goToneOffI), ...
%       'header', mIdC, ...
%       'clims', [-1.5 1.5], ...
%       'sigMask', sigMaskC, ...
%       'absSlope', true, ...
%       'figSaveLogic', true, ...
%       'figName', 'muNoGo_slope_projNoGoAxis');
%

%% Parse inputs
ip = inputParser;
ip.FunctionName = mfilename;

addRequired(ip, 'dataC', @(x) iscell(x) && ~isempty(x));

addParameter(ip, 'timeX', [], @(x) isempty(x) || isnumeric(x));
addParameter(ip, 'colormap', parula(256), @(x) isnumeric(x) && size(x,2)==3);
addParameter(ip, 'axisI', [], @(x) isempty(x) || isnumeric(x) || islogical(x));
addParameter(ip, 'header', [], @(x) isempty(x) || iscell(x) || isstring(x));
addParameter(ip, 'clims', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
addParameter(ip, 'titleStr', '', @(x) ischar(x) || isstring(x));
addParameter(ip, 'xlabelStr', 'Time', @(x) ischar(x) || isstring(x));
addParameter(ip, 'ylabelStr', '', @(x) ischar(x) || isstring(x));
addParameter(ip, 'showColorbar', true, @(x) islogical(x) && isscalar(x));
addParameter(ip, 'sigMask', [], @(x) isempty(x) || iscell(x));
addParameter(ip, 'absSlope', false, @(x) islogical(x) && isscalar(x));
addParameter(ip, 'boundaryColor', 'k', @(x) ischar(x) || isstring(x) || (isnumeric(x) && numel(x)==3));
addParameter(ip, 'boundaryWidth', 1.5, @(x) isnumeric(x) && isscalar(x) && x > 0);

addParameter(ip, 'figSaveLogic', false, @(x) islogical(x) && isscalar(x));
addParameter(ip, 'figSaveDir', ...
    'Z:\Rodent Data\dualImaging_parkj\collectFigure\motifTDR\glmTDR\glmTDR_across_session_slope', ...
    @(x) ischar(x) || isstring(x));
addParameter(ip, 'figName', [], @(x) isempty(x) || ischar(x) || isstring(x));

parse(ip, dataC, varargin{:});
P = ip.Results;

nCells = numel(dataC);

%% Validate dataC
for i = 1:nCells
    if ~isnumeric(dataC{i}) || ndims(dataC{i}) ~= 2
        error('Each cell entry in dataC must contain a 2D numeric matrix.');
    end
end

%% Validate sigMask
useSigMask = ~isempty(P.sigMask);
if useSigMask
    if ~iscell(P.sigMask) || numel(P.sigMask) ~= nCells
        error('sigMask must be a cell array with the same number of elements as dataC.');
    end
end

%% Build stacked matrix
stackedMat = [];
blockStart = zeros(nCells,1);
blockEnd   = zeros(nCells,1);
blockCenter = zeros(nCells,1);
dataC_out = cell(size(dataC));

currentRow = 0;

for i = 1:nCells
    M = dataC{i};

    % Validate timeX
    if ~isempty(P.timeX)
        if size(M,2) ~= numel(P.timeX)
            error('For cell %d, size(M,2) = %d but numel(timeX) = %d.', ...
                i, size(M,2), numel(P.timeX));
        end
    end

    % Apply significance mask if provided
    if useSigMask
        mask = P.sigMask{i};

        if ~islogical(mask)
            error('sigMask{%d} must be a logical matrix.', i);
        end
        if ~isequal(size(mask), size(M))
            error('sigMask{%d} must be the same size as dataC{%d}.', i, i);
        end

        M(~mask) = 0;
    end

    % Apply row selection if requested
    if ~isempty(P.axisI)
        idx = P.axisI;

        % Convert 0/1 numeric vector to logical selector
        if isnumeric(idx) && numel(idx) == size(M,1) && all(ismember(idx, [0 1]))
            idx = logical(idx);
        end

        if islogical(idx)
            if numel(idx) ~= size(M,1)
                error(['Logical axisI must have one element per row. ' ...
                       'Cell %d has %d rows, but numel(axisI) = %d.'], ...
                       i, size(M,1), numel(idx));
            end
        elseif isnumeric(idx)
            if any(idx < 1) || any(mod(idx,1) ~= 0)
                error(['Numeric axisI must contain positive integers only, ' ...
                       'or be a 0/1 selector vector matching the number of rows.']);
            end
        end

        M = M(idx, :);
    end

    % Apply abs() after all indexing/masking if requested
    if P.absSlope
        M = abs(M);
    end

    % Store processed matrix
    dataC_out{i} = M;

    nRowsThis = size(M,1);

    if nRowsThis == 0
        warning('Cell %d contributed zero rows after processing.', i);
        blockStart(i) = currentRow + 1;
        blockEnd(i)   = currentRow;
        blockCenter(i)= NaN;
        continue;
    end

    blockStart(i) = currentRow + 1;
    blockEnd(i)   = currentRow + nRowsThis;
    blockCenter(i)= mean([blockStart(i), blockEnd(i)]);

    stackedMat = [stackedMat; M]; %#ok<AGROW>
    currentRow = currentRow + nRowsThis;
end

if isempty(stackedMat)
    error('No data remained after processing.');
end

%% Create figure
h.fig = figure('Color', 'w');
h.ax = axes('Parent', h.fig);

if isempty(P.timeX)
    h.im = imagesc(stackedMat, 'Parent', h.ax);
    xlim(h.ax, [0.5 size(stackedMat,2)+0.5]);
else
    h.im = imagesc(P.timeX, 1:size(stackedMat,1), stackedMat, 'Parent', h.ax);
end

% First row at the TOP
set(h.ax, 'YDir', 'reverse');

colormap(h.ax, P.colormap);

if ~isempty(P.clims)
    clim(h.ax, P.clims);
end

xlabel(h.ax, P.xlabelStr);
ylabel(h.ax, P.ylabelStr);

% Title behavior
if strlength(string(P.titleStr)) > 0
    title(h.ax, P.titleStr, 'Interpreter', 'none');
else
    title(h.ax, '');
end

hold(h.ax, 'on');

%% Draw mouse boundaries
for i = 1:nCells-1
    if blockEnd(i) > 0
        y = blockEnd(i) + 0.5;
        if isempty(P.timeX)
            xl = [0.5, size(stackedMat,2)+0.5];
        else
            xl = [P.timeX(1), P.timeX(end)];
        end
        plot(h.ax, xl, [y y], '-', ...
            'Color', P.boundaryColor, ...
            'LineWidth', P.boundaryWidth);
    end
end

%% Set y tick labels at block centers
validBlocks = ~isnan(blockCenter);

if ~isempty(P.header)
    if numel(P.header) ~= nCells
        error('header must have one entry per cell in dataC.');
    end
    set(h.ax, 'YTick', blockCenter(validBlocks), ...
              'YTickLabel', P.header(validBlocks), ...
              'TickLabelInterpreter', 'none');
else
    set(h.ax, 'YTick', blockCenter(validBlocks));
end

%% Tone onset line at x = 0
if ~isempty(P.timeX)
    if any(P.timeX == 0)
        xline(h.ax, 0, '--k', 'LineWidth', 1);
    end
end

if P.showColorbar
    colorbar(h.ax);
end

hold(h.ax, 'off');

%% Save figure if requested
h.figPath = [];

if P.figSaveLogic
    if isempty(P.figName)
        error('When figSaveLogic is true, figName must be provided.');
    end

    figSaveDir = char(string(P.figSaveDir));
    figName    = char(string(P.figName));

    if ~exist(figSaveDir, 'dir')
        mkdir(figSaveDir);
    end

    figName = sanitize_filename_local(figName);
    figPath = fullfile(figSaveDir, [figName '.pdf']);

    set(h.fig, 'PaperPositionMode', 'auto');
    print(h.fig, figPath, '-dpdf', '-painters', '-bestfit');

    fprintf('[plotCellImagescStack] saved: %s\n', figPath);
    h.figPath = figPath;
end

%% Output handles and metadata
h.stackedMat = stackedMat;
h.blockStart = blockStart;
h.blockEnd = blockEnd;
h.blockCenter = blockCenter;

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