function plotMaskedXcorrHeatmap(ax, A, motifOrder, climVals, blockBounds)
% Masked xcorr heatmap: NaN entries render as the grey axis background via
% AlphaData rather than as a colormap value.
%
% motifOrder  : the motif indices in plotting order, used for tick labels
%               (pass leafOrderGoNoGo, and index A by it before calling).
% blockBounds : optional vector of positions after which to draw a black
%               divider, marking the Go / unlabeled / NoGo blocks. Omit or
%               pass [] for no dividers.
if nargin < 5, blockBounds = []; end

imAlpha = ones(size(A));
imAlpha(isnan(A)) = 0;
imagesc(ax, A, 'AlphaData', imAlpha, climVals);
set(ax, 'Color', [0.85 0.85 0.85]);
axis(ax, 'square'); axis(ax, 'xy');
pbaspect(ax, [1 1 1]);
colormap(ax, parula);
colorbar(ax);
set(ax, 'XTick', 1:numel(motifOrder), 'XTickLabel', string(motifOrder), ...
        'YTick', 1:numel(motifOrder), 'YTickLabel', string(motifOrder));
xtickangle(ax, 90);
xlabel(ax, 'Motif (Go | unlabeled | NoGo)');
ylabel(ax, 'Motif (Go | unlabeled | NoGo)');

hold(ax, 'on');
for b = blockBounds(:)'
    xline(ax, b + 0.5, 'k-', 'LineWidth', 1.2);
    yline(ax, b + 0.5, 'k-', 'LineWidth', 1.2);
end
hold(ax, 'off');
end
