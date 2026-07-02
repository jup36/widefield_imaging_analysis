function plotTrialLevelXcorrHeatmap(rez, condName)

if isempty(rez.meanR)
    return;
end

figure;
imagesc(rez.lagSec, 1:size(rez.meanR, 1), rez.meanR);
axis xy;
colorbar;
xline(0, 'w--', 'LineWidth', 1.5);
xlabel('Lag, sec');
ylabel('Motif');
title(sprintf('Trial-level H-DA xcorr: %s', strrep(condName, 'I', '')));
box off;

end
