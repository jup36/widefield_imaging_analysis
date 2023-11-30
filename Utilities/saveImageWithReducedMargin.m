% Function to save an image with reduced margins
function saveImageWithReducedMargin(imageData, filename)
    fig = figure('Visible', 'off'); % Create a new invisible figure
    ax = axes('Parent', fig, 'Unit', 'normalized', 'Position', [0 0 1 1]);
    imshow(imageData, 'Parent', ax);
    set(fig, 'PaperUnits', 'points', 'PaperPosition', [0 0 size(imageData, 2), size(imageData, 1)]);
    saveas(fig, filename, 'png');
    close(fig); % Close the figure
end