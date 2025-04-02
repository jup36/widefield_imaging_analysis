function h = plotLearningCurveCell(datC, cMat)
    % Function to plot learning curves for multiple animals across days
    % using a common x-axis. The x-axis is divided into days with a gap
    % between them. Each day is allocated the same x-range based on the maximum
    % number of blocks recorded for that day across all animals.
    %
    % Inputs:
    %   datC - Nx1 cell array where each cell corresponds to one animal.
    %          Each animal is a 1xD cell array (D = number of days), where each
    %          cell is a 1xM array of data (M can vary per day and animal).
    %   cMat - Nx3 matrix specifying RGB colors for each animal.
    
    if nargin < 2
        error('Both datC and cMat must be provided.');
    end
    
    numMice = length(datC);
    
    % Determine the total number of days (use the maximum over animals)
    D = max(cellfun(@length, datC));
    
    % For each day, determine the maximum number of blocks across all animals.
    nBlocks = zeros(1, D);
    for d = 1:D
        blocks = [];
        for iMouse = 1:numMice
            if length(datC{iMouse}) >= d && ~isempty(datC{iMouse}{d})
                blocks(end+1) = length(datC{iMouse}{d}); %#ok<AGROW>
            end
        end
        if isempty(blocks)
            nBlocks(d) = 0;
        else
            nBlocks(d) = max(blocks);
        end
    end
    
    % Compute global x positions for each day.
    gapBetweenDays = 1; % Space between days
    globalDayX = cell(1, D);
    dayCenters = zeros(1, D);
    currentX = 1; % Start with an extra gap (previously 0, now 1)

    for d = 1:D
        % Start the day after a gap.
        dayStart = currentX + gapBetweenDays + 1;
        % Allocate positions for this day (using the max number of blocks).
        globalDayX{d} = dayStart:(dayStart + nBlocks(d) - 1);
        dayCenters(d) = mean(globalDayX{d});
        % Update currentX to the end of this day.
        currentX = dayStart + nBlocks(d) - 1;
    end
    
    % Create figure and hold on.
    h = figure; hold on;
    set(gca, 'XColor', 'k');  % Ensure x-axis is visible.
    
    % Preallocate legend handles and entries.
    legendHandles = gobjects(numMice, 1);
    legendEntries = cell(numMice, 1);
    
    % Plot each animal using the global x positions.
    for iMouse = 1:numMice
        mouseData = datC{iMouse};
        firstHandle = [];  % To capture the first plotted line for legend.
        
        % Loop over days that this animal has data.
        for d = 1:length(mouseData)
            if isempty(mouseData{d})
                continue;
            end
            dataDay = mouseData{d};
            m = length(dataDay);
            % Use the first m x positions from the global positions for day d.
            % (If m < nBlocks(d), the animal's data is left-aligned within the day.)
            xPositions = globalDayX{d}(1:m);
            % Plot this day's data as a separate segment.
            hPlot = plot(xPositions, dataDay, '-o', ...
                'Color', cMat(iMouse, :), ...
                'MarkerFaceColor', cMat(iMouse, :), ...
                'MarkerEdgeColor', 'none', 'LineWidth', 1.5);
            
            if isempty(firstHandle)
                firstHandle = hPlot;
            end
        end
        % Store handle and legend entry for this animal.
        legendHandles(iMouse) = firstHandle;
        legendEntries{iMouse} = sprintf('Mouse %d', iMouse);
    end
    
    % Set x-axis ticks at the center of each day and label them.
    set(gca, 'XTick', dayCenters, ...
             'XTickLabel', arrayfun(@(d) sprintf('Day %d', d), 1:D, 'UniformOutput', false), ...
             'XTickLabelRotation', 45); % Ensure labels are readable

    % Adjust x-limits to fully include the last day's label, with space at the left.
    xlim([min(cellfun(@min, globalDayX)) - 1, max(cellfun(@max, globalDayX)) + gapBetweenDays]);

    xlabel('Days'); ylabel('Performance');
    title('Learning Curve');
    legend(legendHandles, legendEntries, 'Location', 'best');
    box off;
    hold off;
    grid on; 
    set(gca, 'TickDir', 'out')
end
