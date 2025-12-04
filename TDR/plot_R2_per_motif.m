function plot_R2_per_motif(R2_per)
%PLOT_R2_PER_MOTIF  Visualize cross-validated R² per motif as a bar plot.
%
%   plot_R2_per_motif(R2_per)
%
%   Inputs
%     R2_per : [1 x nMotifs] vector of cross-validated R² per motif

nMotifs = numel(R2_per);

figure('Color','w');
bar(R2_per, 'FaceColor',[0.3 0.5 0.8], 'EdgeColor','none');
hold on;
yline(0, 'k');  % baseline
xlabel('Motif');
ylabel('R^{2}', 'FontWeight','bold', 'Interpreter','tex');
title('Cross-validated R^{2} per motif', 'FontWeight','bold', 'Interpreter','tex');
xlim([0 nMotifs+1]);
ylim([0 max(R2_per)*1.1]);
box off;
grid on;

% Label x-ticks as "Motif 1", "Motif 2", ...
xticklabels(arrayfun(@(k) sprintf('Motif %d', k), 1:nMotifs, 'UniformOutput', false));
xticks(1:nMotifs);
xtickangle(45);

end
