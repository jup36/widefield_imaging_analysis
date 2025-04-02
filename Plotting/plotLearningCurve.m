function h = plotLearningCurve(dat, cMat)
    % Function to plot learning curves for multiple animals across days.
    % Each animal has one data point per day, and the data points are connected.
    %
    % Inputs:
    %   dat  - Nx1 cell array, where each cell is a 1xD numeric array 
    %          (one value per day for that animal).
    %   cMat - Nx3 matrix specifying RGB colors for each animal.
    
    if nargin < 2
        error('Both dat and cMat must be provided.');
    end
    
    numMice = length(dat);
    
    % Determine the total number of days (use the maximum across animals)
    D = max(cellfun(@length, dat));
    
    % Create figure and hold on.
    h = figure; hold on;
    set(gca, 'XColor', 'k');  % Ensure x-axis is visible.
    
    % Preallocate legend handles and entries.
    legendHandles = gobjects(numMice, 1);
    legendEntries = cell(numMice, 1);
    
    % Loop through each mouse and plot their data as a continuous line.
    for iMouse = 1:numMice
        mouseData = dat{iMouse}; % Get data for this mouse
        numDays = length(mouseData); % Number of days for this mouse
        xDays = 1:numDays; % X-axis values (Day indices)
        
        % Plot the mouse's data as a connected line
        hPlot = plot(xDays, mouseData, '-o', ...
            'Color', cMat(iMouse, :), ...
            'MarkerFaceColor', cMat(iMouse, :), ...
            'MarkerEdgeColor', 'none', 'LineWidth', 1.5);
        
        % Store the first valid plot handle for legend
        legendHandles(iMouse) = hPlot;
        legendEntries{iMouse} = sprintf('Mouse %d', iMouse);
    end
    
    % Set x-axis ticks to match day numbers
    set(gca, 'XTick', 1:D, ...
             'XTickLabel', arrayfun(@(d) sprintf('Day %d', d), 1:D, 'UniformOutput', false), ...
             'XTickLabelRotation', 45); % Rotate labels for better readability

    % Set limits
    xlim([0.5, D + 0.5]); % Add a small margin around the days
    
    xlabel('Days'); ylabel('Performance');
    title('Learning Curve');
    legend(legendHandles, legendEntries, 'Location', 'best');
    box off;
    hold off;
    grid on; 
    set(gca, 'TickDir', 'out')
end
