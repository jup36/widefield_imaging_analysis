function figH = plotPairedLines(v1, v2, okMask, isFast, fastColor, slowColor, ...
    xLabels, yLabelStr, titleStr, yTickStep)
% One line per animal connecting its two paired values, colored by learning
% group, with the group mean +/- SEM overlaid. Uses the same mask as the
% corresponding t-test so figure and test cover identical animals.
%
% yTickStep (optional): y-axis tick interval, e.g. 0.05. The axis limits
% are snapped OUTWARD to whole multiples of the step, so ticks land on
% round numbers and no data point falls outside the axis. Setting YTick
% alone would leave MATLAB's auto-limits in place and leave blank space at
% the ends. Omit or pass [] to keep MATLAB's automatic ticks.

if nargin < 10, yTickStep = []; end

figH = figure('Color','w', 'Position', [200 200 420 460]);
ax = axes(figH);
hold(ax, 'on');

for a = 1:numel(v1)
    if ~okMask(a), continue; end
    c = ternary_local(isFast(a), fastColor, slowColor);
    plot(ax, [1 2], [v1(a), v2(a)], '-', 'Color', [c 0.5], 'LineWidth', 1.1);
    scatter(ax, [1 2], [v1(a), v2(a)], 160, c, 'filled', 'MarkerFaceAlpha', 0.8);
end

gm = [mean(v1(okMask),'omitnan'), mean(v2(okMask),'omitnan')];
gs = [semLocal(v1(okMask)),       semLocal(v2(okMask))];
errorbar(ax, [1 2], gm, gs, 'ko-', 'MarkerFaceColor','k', 'MarkerSize', 14, ...
    'LineWidth', 1.6, 'CapSize', 10);

set(ax, 'XTick', [1 2], 'XTickLabel', xLabels, 'XLim', [0.7 2.3]);

if ~isempty(yTickStep) && isfinite(yTickStep) && yTickStep > 0
    yl = ylim(ax);
    yl = [floor(yl(1)/yTickStep), ceil(yl(2)/yTickStep)] * yTickStep;
    set(ax, 'YLim', yl, 'YTick', yl(1):yTickStep:yl(2));
end

ylabel(ax, yLabelStr);
title(ax, titleStr, 'Interpreter', 'tex');
grid(ax, 'on'); box(ax, 'on');
hold(ax, 'off');
end