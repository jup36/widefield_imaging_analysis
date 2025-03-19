function plotHwithTaskEvents(saveFigFullPath, hmat, motif_ind, timestamps, trialId, tbytDat, trTypeI)

cTable = slanCM('set3', 11); 
cTable = [[0 0 0]; cTable];  % Ensuring black for first color

assert(length(motif_ind)<=5) % Prevent overcrowding the plot
ymax = max(max(hmat(motif_ind, :))); 

% Iterate through pulse trains
fig = figure; hold on;
% Get the default figure position
default_pos = get(fig, 'Position'); 
% Double the width while keeping height the same
new_pos = [default_pos(1), default_pos(2), default_pos(3) * 2, default_pos(4)];
set(fig, 'Position', new_pos); % Apply new figure size
legend_entries = []; % Store legend handles
legend_labels = {};  % Store legend labels

for i = 1:length(motif_ind)
    temph = hmat(motif_ind(i), :);
    t_int = timestamps(1):0.01:timestamps(end);
    temph_int = interp1(timestamps, temph, t_int, 'linear'); % Linear interpolation
    label = sprintf("h_motif#%d", motif_ind(i)); 
    
    % Assign colors from cTable
    color_i = cTable(motif_ind(i), :);
    
    % Plot with color instead of line styles
    h_plot = plot(t_int, temph_int, 'Color', color_i, 'LineWidth', 1, 'DisplayName', label);
    legend_entries(end+1) = h_plot;
    legend_labels{end+1} = label;
end

for ii = trialId
    if trTypeI.goI(ii)
        e1 = tbytDat(ii).evtOn;
        e2 = tbytDat(ii).evtOff;
        h_fill = fill([e1 e2 e2 e1], [0 0 ymax ymax], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        set(h_fill, 'HandleVisibility', 'off'); % Hide from legend
        % Plot hit licks
        hitLicks = tbytDat(ii).hitLicks + e1;
        for lick = hitLicks
            plot([lick lick], [0 ymax], 'b-', 'LineWidth', 1, 'HandleVisibility', 'off', 'Color', [0 0 1 0.6]);
        end

    elseif trTypeI.nogoI(ii)
        e1 = tbytDat(ii).evtOn;
        e2 = tbytDat(ii).evtOff;
        h_fill = fill([e1 e2 e2 e1], [0 0 ymax ymax], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        set(h_fill, 'HandleVisibility', 'off'); % Hide from legend
        % Plot false alarm licks
        faLicks = tbytDat(ii).faLicks + e1;
        for lick = faLicks
            plot([lick lick], [0 ymax], 'b-', 'LineWidth', 1, 'HandleVisibility', 'off', 'Color', [1 0 0 0.6]);
        end
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
