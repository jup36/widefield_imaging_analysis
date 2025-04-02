function plotHwithTaskEventsTrial(saveFigFullPath, hmat, motif_ind, timestamps, trialId, tbytDat, trTypeI)

%cTable = slanCM('set3', 11); 
cTable = slanCM('dark2', 8); 
cTable = [[0 0 0]; cTable];  % Ensuring black for first color

%assert(length(motif_ind)<=5) % Prevent overcrowding the plot
ymax = max(max(hmat(motif_ind, :))); 
evtOnT = tbytDat(trialId).evtOn; 
t_int = timestamps(1):0.01:timestamps(end);

% Iterate through pulse trains
fig = figure; hold on;
legend_entries = []; % Store legend handles
legend_labels = {};  % Store legend labels

for i = 1:length(motif_ind)
    temph = hmat(motif_ind(i), :);
    temph_int = interp1(timestamps, temph, t_int, 'linear'); % Linear interpolation
    label = sprintf("h_motif#%d", motif_ind(i));

    % Assign colors from cTable
    color_i = cTable(motif_ind(i), :);
    
    % Plot with color instead of line styles
    h_plot = plot(t_int - evtOnT, temph_int, 'Color', color_i, 'LineWidth', 1.5, 'DisplayName', label);
    legend_entries(end+1) = h_plot;
    legend_labels{end+1} = label;
end

e1 = 0; % Start time of event
e2 = tbytDat(trialId).evtOff - tbytDat(trialId).evtOn; % Duration of event

if trTypeI.goI(trialId)
    h_fill = fill([e1 e2 e2 e1], [0 0 ymax ymax], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    set(h_fill, 'HandleVisibility', 'off'); % Hide from legend
    % Plot hit licks
    hitLicks = tbytDat(trialId).hitLicks + e1;
    for ii = 1:length(hitLicks)
        plot([hitLicks(ii) hitLicks(ii)], [0 ymax], 'b-', 'LineWidth', 1, 'HandleVisibility', 'off', 'Color', [0 0 1 0.6]);
    end

elseif trTypeI.nogoI(trialId)
    h_fill = fill([e1 e2 e2 e1], [0 0 ymax ymax], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    set(h_fill, 'HandleVisibility', 'off'); % Hide from legend
    % Plot false alarm licks
    faLicks = tbytDat(trialId).faLicks + e1;
    for ii = 1:length(faLicks)
        plot([faLicks(ii) faLicks(ii)], [0 ymax], 'b-', 'LineWidth', 1, 'HandleVisibility', 'off', 'Color', [1 0 0 0.6]);
    end
end

% Labels and legend
xlabel('Time (s)');
ylabel('Temporal Weightings');
legend(legend_entries, legend_labels, 'Interpreter', 'none', 'Location', 'eastoutside');
%grid on;
hold off;
ylim([0 ymax])
axis tight
set(gca, 'TickDir', 'out')

% Save figure 
if ~isempty(saveFigFullPath)
    print(fig, saveFigFullPath, '-dpdf', '-vector', '-bestfit');
end

end
