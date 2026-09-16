function plotPairedLinesAx(ax, V, okMask, isFast, fastColor, slowColor, xLabels, yLabelStr, varargin)
% One line per animal connecting its values across the columns of V
% (animal x level), colored by learning group, with the group mean +/- SEM
% overlaid. Generalized from two levels to any number, so the same helper
% serves the across-class plot and the early-vs-late panels.
%
% Uses the same mask as the corresponding test, so figure and test cover
% identical animals.
%
% OPTIONAL NAME-VALUE PAIRS
%   'YLim'      : [lo hi] explicit y-axis limits, e.g. [-0.1 0.3]. Applied
%                 before ticks, so YTickStep lays ticks inside this range
%                 rather than overriding it. Default [] = automatic.
%   'YTickStep' : scalar tick interval, e.g. 0.05. With no YLim, the
%                 automatic limits are snapped OUTWARD to whole multiples
%                 of the step so ticks land on round numbers and no point
%                 falls outside the axis. Default [] = automatic ticks.
%
% Backward compatible with the old positional ninth argument: a bare
% numeric scalar in that slot is still read as YTickStep.

p = inputParser;
p.addParameter('YLim', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 2 && x(2) > x(1)));
p.addParameter('YTickStep', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));

% Old call style: plotPairedLinesAx(..., yLabelStr, yTickStep)
if numel(varargin) == 1 && (isempty(varargin{1}) || isnumeric(varargin{1}))
    varargin = {'YTickStep', varargin{1}};
end
p.parse(varargin{:});

yLimIn    = p.Results.YLim;
yTickStep = p.Results.YTickStep;

nLevels = size(V, 2);
x = 1:nLevels;

hold(ax, 'on');
for a = 1:size(V, 1)
    if ~okMask(a), continue; end
    c = ternary_local(isFast(a), fastColor, slowColor);
    plot(ax, x, V(a, :), '-', 'Color', [c 0.5], 'LineWidth', 1.1);
    scatter(ax, x, V(a, :), 45, c, 'filled', 'MarkerFaceAlpha', 0.8);
end

gm = mean(V(okMask, :), 1, 'omitnan');
gs = arrayfun(@(k) semLocal(V(okMask, k)), 1:nLevels);
errorbar(ax, x, gm, gs, 'ko-', 'MarkerFaceColor','k', 'MarkerSize', 7, ...
    'LineWidth', 1.6, 'CapSize', 10);

set(ax, 'XTick', x, 'XTickLabel', xLabels, 'XLim', [0.7, nLevels + 0.3]);

% -------- y-axis: explicit limits win; otherwise snap to the tick step --------
if ~isempty(yLimIn)
    % An explicit range can clip data -- say so rather than hiding points.
    plotted = V(okMask, :);
    nOut = sum(plotted(:) < yLimIn(1) | plotted(:) > yLimIn(2));
    if nOut > 0
        warning('plotPairedLinesAx: %d data point(s) fall outside the requested YLim [%g %g] and will be clipped.', ...
            nOut, yLimIn(1), yLimIn(2));
    end
    set(ax, 'YLim', yLimIn);
    if ~isempty(yTickStep)
        firstTick = ceil(yLimIn(1)/yTickStep) * yTickStep;
        set(ax, 'YTick', firstTick:yTickStep:yLimIn(2));
    end
elseif ~isempty(yTickStep)
    yl = ylim(ax);
    yl = [floor(yl(1)/yTickStep), ceil(yl(2)/yTickStep)] * yTickStep;
    set(ax, 'YLim', yl, 'YTick', yl(1):yTickStep:yl(2));
end

ylabel(ax, yLabelStr);
grid(ax, 'on'); box(ax, 'on');
hold(ax, 'off');
end