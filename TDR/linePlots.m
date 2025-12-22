function h = linePlots(dat, varargin)
% LINEPLOTS  Plot one or more lines with distinct colors and optional labels.
%
%   h = linePlots(dat, 'legend', {...}, 'title', '...', 
%                        'xTickLabel', {...}, 'xLabel', '...', 
%                        'yLabel', '...', 'yLim', [min max],
%                        'figVisible', true/false)
%
% INPUT
%   dat : numeric vector or matrix (N x T)
%
% OUTPUT
%   h.fig, h.ax, h.lines, h.scatter

% -------- Parse inputs --------
p = inputParser;
p.addParameter('legend', {}, @(x) iscell(x) || isempty(x));
p.addParameter('title', '', @(x) ischar(x) || isstring(x));
p.addParameter('xTickLabel', {}, @(x) iscell(x) || isempty(x));
p.addParameter('xLabel', '', @(x) ischar(x) || isstring(x));
p.addParameter('yLabel', '', @(x) ischar(x) || isstring(x));
p.addParameter('yLim', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
p.addParameter('figVisible', true, @(x) islogical(x) || isnumeric(x));  % <--- NEW
p.parse(varargin{:});
opt = p.Results;

% -------- Prepare data --------
dat = squeeze(dat);
if isvector(dat)
    dat = dat(:)';  % single row
end
[nLines, nPts] = size(dat);
x = 1:nPts;

% -------- Create figure --------
figVisibility = ternary(opt.figVisible, 'on', 'off');   % small helper below
h.fig = figure('Color','w', 'Visible', figVisibility);  % <--- UPDATED
h.ax  = axes('Parent',h.fig); 
hold(h.ax,'on');

% Distinct colors
cols = lines(nLines);

% -------- Plot each line --------
for i = 1:nLines
    h.lines(i) = plot(h.ax, x, dat(i,:), '-', ...
        'Color', cols(i,:), 'LineWidth', 1.5, ...
        'DisplayName', getDisplayName(opt.legend, i));

    h.scatter(i) = scatter(h.ax, x, dat(i,:), 50, cols(i,:), ...
        'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', 'none', ...
        'HandleVisibility','off');
end

% -------- Axes and appearance --------
grid(h.ax, 'on');
box(h.ax, 'off');
set(h.ax,'TickDir','out','LineWidth',1,'TickLabelInterpreter','none');

if ~isempty(opt.yLim)
    ylim(h.ax, opt.yLim);
end

if ~isempty(opt.xTickLabel)
    xticks(h.ax, 1:numel(opt.xTickLabel));
    xticklabels(h.ax, opt.xTickLabel);
end

xlabel(h.ax, opt.xLabel, 'Interpreter','none');
ylabel(h.ax, opt.yLabel, 'Interpreter','none');

if ~isempty(opt.title)
    title(h.ax, opt.title, 'Interpreter','none');
end

if ~isempty(opt.legend)
    legend(h.ax, h.lines, opt.legend, ...
        'Location','best', 'Box','off', 'Interpreter','none');
end

end % function


% ===== helper: legend name retrieval =====
function name = getDisplayName(legendList, idx)
if isempty(legendList)
    name = sprintf('Line %d', idx);
elseif numel(legendList) >= idx
    name = legendList{idx};
else
    name = sprintf('Line %d', idx);
end
end

% ===== helper: ternary operator (inline) =====
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
