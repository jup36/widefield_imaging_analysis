function h = plotLearningCurve(datC, cMat)
    % Function to plot learning curves for multiple animals across days
    % Inputs:
    %   datC - Nx1 cell array where each entry corresponds to a mouse
    %          Each mouse has a 1xD cell array, where each entry is a 1x3 array (data from 3 blocks per session).
    %   cMat - Nx3 matrix specifying RGB colors for each animal
    
    if nargin < 2
        error('Both datC and cMat must be provided.');
    end
    
    numMice = length(datC);  % Number of mice
    h = figure; hold on; % Initialize figure
    set(gca, 'XTick', [], 'XColor', 'w'); % Remove x-axis ticks for clarity
    
    legendEntries = cell(numMice, 1); % Store legend labels
    
    % Loop through each mouse
    for iMouse = 1:numMice
        mouseData = datC{iMouse}; % Get data for this mouse
        numDays = length(mouseData); % Number of days
        
        xPos = []; % X-axis positions for plotting
        yPos = []; % Y-axis values
        
        for iDay = 1:numDays
            if isempty(mouseData{iDay})
                continue; % Skip if there's no data for the day
            end
            % Define x positions for three blocks
            xBase = (iDay - 1) * 4; % Leave a gap of 1 between days (3 blocks + 1 space)
            xPos = [xPos, xBase + (1:3)];
            yPos = [yPos, mouseData{iDay}]; % Append y values
        end
        
        % Plot each mouse's data
        plot(xPos, yPos, '-o', 'Color', cMat(iMouse, :), ...
            'MarkerFaceColor', cMat(iMouse, :), 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        
        % Store legend entry
        legendEntries{iMouse} = sprintf('Mouse %d', iMouse);
    end
    
    xlabel('Days (with gaps)'); ylabel('Performance');
    title('Learning Curve');
    legend(legendEntries, 'Location', 'best');
    box off;
    hold off;
end
