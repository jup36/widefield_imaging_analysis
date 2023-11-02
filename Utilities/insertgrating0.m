function insertgrating0(figHandle, img)
    % Ensure we're working with the correct figure
    figure(figHandle);

    % Preserve the current hold state
    holdState = ishold;
    hold on; % Ensure that the lines are added to the existing figure

    % Grating parameters
    lineSpacing = 1.5; % Spacing between lines in pixels
    lineLength = size(img, 2); % Length of lines equals the width of the image
    numLines = ceil(size(img, 1) / lineSpacing); % Number of lines based on height

    % Draw horizontal lines
    for i = 1:numLines
        % Calculate y position for the line
        yPos = i * lineSpacing;

        % Draw the line
        line([1, lineLength], [yPos, yPos], 'Color', 'white', 'LineWidth', 1);
    end

    % Return to previous hold state
    if ~holdState
        hold off;
    end
end
