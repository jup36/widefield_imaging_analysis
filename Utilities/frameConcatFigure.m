
filePathFrame = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA008/DA008_101823/DA008_101823_trials/trial_21'; 
[file_list, folder_list] = GrabFiles_sort_trials('FrameWithText', 0, {filePathFrame}); 

frames = 75:94; 
counter = 1; 

for fr = frames
    tempImg = imread(file_list{fr});
    [imgHeight, imgWidth, ~] = size(tempImg);

    % Define margins for cropping
    cutEachX = 120;
    cutEachY = 32;

    % Define the cropping rectangle
    cropRect = [cutEachX, cutEachY, imgWidth - 2*cutEachX, imgHeight - 4*cutEachY];

    % Crop the image
    croppedImage = imcrop(tempImg, cropRect);

    % Store cropped images in the appropriate row based on the counter
    if counter <= 10  % First 10 frames
        row1Images{counter} = croppedImage;
    else  % Next 10 frames
        row2Images{counter - 10} = croppedImage;
    end

    % Increment the counter
    counter = counter + 1;
end

% Concatenate images horizontally for each row
concatRow1 = cat(2, row1Images{:});
concatRow2 = cat(2, row2Images{:});

% Define the file names
firstRowFileName = fullfile(filePathFrame, 'first_row_concat');
secondRowFileName = fullfile(filePathFrame, 'second_row_concat');

% Save the first row with reduced margins
saveImageWithReducedMargin(concatRow1, firstRowFileName);

% Save the second row with reduced margins
saveImageWithReducedMargin(concatRow2, secondRowFileName);



% Function to save an image with reduced margins
function saveImageWithReducedMargin(imageData, filename)
    fig = figure('Visible', 'off'); % Create a new invisible figure
    ax = axes('Parent', fig, 'Unit', 'normalized', 'Position', [0 0 1 1]);
    imshow(imageData, 'Parent', ax);
    set(fig, 'PaperUnits', 'points', 'PaperPosition', [0 0 size(imageData, 2), size(imageData, 1)]);
    saveas(fig, filename, 'png');
    close(fig); % Close the figure
end



