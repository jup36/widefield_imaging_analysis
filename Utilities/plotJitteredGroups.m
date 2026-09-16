function plotJitteredGroups(ax, fastVals, slowVals, fastColor, slowColor)
% Jittered strip plot: fast at x=1, slow at x=2, small random horizontal
% jitter for visibility, group mean +/- SEM overlaid as a black errorbar.
jitterWidth = 0.15;

xFast = 1 + (rand(numel(fastVals), 1) - 0.5) * jitterWidth * 2;
xSlow = 2 + (rand(numel(slowVals), 1) - 0.5) * jitterWidth * 2;

hold(ax, 'on');
scatter(ax, xFast, fastVals, 40, fastColor, 'filled', 'MarkerFaceAlpha', 0.7);
scatter(ax, xSlow, slowVals, 40, slowColor, 'filled', 'MarkerFaceAlpha', 0.7);

meanFast = mean(fastVals, 'omitnan'); semFast = std(fastVals, 'omitnan') / sqrt(sum(isfinite(fastVals)));
meanSlow = mean(slowVals, 'omitnan'); semSlow = std(slowVals, 'omitnan') / sqrt(sum(isfinite(slowVals)));

errorbar(ax, 1, meanFast, semFast, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6, 'LineWidth', 1.3, 'CapSize', 8);
errorbar(ax, 2, meanSlow, semSlow, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6, 'LineWidth', 1.3, 'CapSize', 8);

set(ax, 'XTick', [1 2], 'XTickLabel', {'Fast', 'Slow'});
xlim(ax, [0.5 2.5]);
grid(ax, 'on'); box(ax, 'on');
hold(ax, 'off');
end