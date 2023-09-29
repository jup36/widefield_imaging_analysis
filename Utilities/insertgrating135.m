function insertgrating135(figHandle, img)
% Ensure we're working with the correct figure
figure(figHandle);

% Grating parameters
lineSpacing = 1.5; % Spacing between lines in pixels
lineLength = sqrt((size(img,1)/15)^2 + (size(img,2)/15)^2); % Diagonal length
numLines = ceil(lineLength/lineSpacing); % Number of lines needed

% Draw 45-degree angled lines
for i = 1:numLines
    % Calculate starting point for the line
    startX = i * lineSpacing;
    startY = 1;
    
    % Adjust starting point if it's outside the img's width
    if startX > size(img, 2)
        startY = startX - size(img, 2);
        startX = size(img, 2);
    end
    
    % Calculate ending point for the line
    endX = startX + lineLength;
    endY = startY + lineLength;
    
    % Draw the line
    line([startX, endX], [startY, endY], 'Color', 'white', 'LineWidth', 3);
end

hold off;
end


