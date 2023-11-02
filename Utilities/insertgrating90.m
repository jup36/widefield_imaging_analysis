function insertgrating90(figHandle, img)
    % Ensure we're working with the correct figure
    figure(figHandle);

    % Preserve the current hold state
    holdState = ishold;
    hold on; % Ensure that the lines are added to the existing figure

    % Grating parameters
    lineSpacing = 1.5; % Spacing between lines in pixels
    lineLength = size(img, 1); % Length of lines equals the height of the image
    numLines = ceil(size(img, 2) / lineSpacing); % Number of lines based on width

    % Draw vertical lines
    for i = 1:numLines
        % Calculate x position for the line
        xPos = i * lineSpacing;

        % Draw the line
        line([xPos, xPos], [1, lineLength], 'Color', 'white', 'LineWidth', 1);
    end

    % Return to previous hold state
    if ~holdState
        hold off;
    end
end
