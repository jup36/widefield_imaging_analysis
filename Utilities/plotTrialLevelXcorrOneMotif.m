function plotTrialLevelXcorrOneMotif(rez, motifK, condName)

if isempty(rez.meanR)
    return;
end

if motifK > size(rez.meanR, 1)
    warning('motifK exceeds available motifs.');
    return;
end

m = rez.meanR(motifK, :);
s = rez.semR(motifK, :);

figure;
hold on;

fill([rez.lagSec fliplr(rez.lagSec)], ...
    [m - s fliplr(m + s)], ...
    [0.8 0.8 0.8], ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.5);

plot(rez.lagSec, m, 'k', 'LineWidth', 2);
xline(0, '--');

xlabel('Lag, sec');
ylabel('Trial-averaged correlation');
title(sprintf('Trial-level H-DA xcorr: %s, motif %d', ...
    strrep(condName, 'I', ''), motifK));
box off;

end