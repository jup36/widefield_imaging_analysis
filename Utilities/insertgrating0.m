function insertgrating0(figHandle, img)
    % Ensure we're working with the correct figure
    figure(figHandle);
    hold on; % Ensure that the lines are added to the existing figure
    
    % Determine the extents for the lines based on image size
    xExtent = size(img, 2) * 0.10; % 10% of the image's x dimension
    ySpacing = size(img, 1) * 0.10 / 4; % Divide 10% of y by 4 for each line

    % Draw horizontal lines at the top left corner
    for i = 1:4
        yPos = i * ySpacing;
        line([1, xExtent], [yPos, yPos], 'Color', 'white', 'LineWidth', 1.5);
    end
    
    hold off;
end
