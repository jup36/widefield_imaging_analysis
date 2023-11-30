function insertgrating90(figHandle, img)
    % Ensure we're working with the correct figure
    figure(figHandle);
    hold on; % Ensure that the lines are added to the existing figure
    
    % Determine the extents for the lines based on image size
    yExtent = size(img, 1) * 0.10; % 10% of the image's y dimension
    xSpacing = size(img, 2) * 0.10 / 4; % Divide 10% of x by 4 for each line

    % Draw vertical lines at the top left corner
    for i = 1:4
        xPos = i * xSpacing;
        line([xPos, xPos], [1, yExtent], 'Color', 'white', 'LineWidth', 1.5);
    end
    
    hold off;
end