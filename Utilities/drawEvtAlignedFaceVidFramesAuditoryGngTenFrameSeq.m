function drawEvtAlignedFaceVidFramesAuditoryGngTenFrameSeq(pngPath, vidName, evtName, every, pethFrameI, frameTrelEvt, frameLickI, frameWaterI, frameAirpuffI)

figSavePath = fullfile(pngPath, 'frame_sequence'); 
if exist(figSavePath, 'dir')~=7
    mkdir(figSavePath)
end

% Grab all sorted PNG files
pngs_all = GrabFiles_sort_trials('faceDff_frame_', 0, {pngPath});
pngs = pngs_all(pethFrameI);

% Ensure the lengths match
assert(length(pngs) == sum(pethFrameI), 'Mismatch in the number of frames.');

% Define margins and text box dimensions
margin = 20; % Margin from the edge in pixels
textBoxHeight = 20; % Estimate text box height (may need adjustment based on font size)
textBoxWidth = 70; % Estimate text box width (may need adjustment based on font size and text length)

% Initialize counters and arrays for concatenated images
concat = [];
frameCount = 0;
batchCount = 1;

for f = 1:length(pngs) % Loop through frames
    if mod(f, every) == 0
        % Load PNG image
        imgFull = imread(pngs{f});
        [~, imgWidthFull, ~] = size(imgFull); % Get the size of the image

        % Crop the face video
        img = imgFull(:, 1:floor(imgWidthFull/2)-1, :);
        [imgHeight, imgWidth, ~] = size(img);

        %% Label the image
        frameTimePos = [imgWidth - textBoxWidth - margin, imgHeight - textBoxHeight - margin]; % Bottom right
        lickPos = [margin, margin]; % Top left
        waterPos = [margin, 2 * margin + textBoxHeight]; % Below the 'Lick' position

        % Add frame time to the image
        frameTimeStr = sprintf('%.2fs', frameTrelEvt(f));
        imgWithText = insertText(img, frameTimePos, frameTimeStr, 'FontSize', 22, 'TextColor', 'white', 'BoxColor', [255, 0, 255], 'BoxOpacity', 0.3);

        % Add 'Lick' annotation and draw a rectangle around the edge if applicable
        if frameLickI(f)
            imgWithText = insertText(imgWithText, lickPos, 'Lick', 'FontSize', 22, 'TextColor', 'white');
            rectPosition = [1, 1, imgWidth-1, imgHeight-1]; % Rectangle position [x, y, width, height]
            imgWithText = insertShape(imgWithText, 'Rectangle', rectPosition, 'Color', 'cyan', 'LineWidth', 5);
        end

        % Add 'Water' or 'Airpuff' annotation
        if frameWaterI(f)
            imgWithText = insertText(imgWithText, waterPos, 'Water', 'FontSize', 22, 'TextColor', 'white', 'BoxColor', [0, 255, 255], 'BoxOpacity', 0.6);
        elseif frameAirpuffI(f) % Assuming Water and Airpuff don't occur together
            imgWithText = insertText(imgWithText, waterPos, 'Airpuff', 'FontSize', 20, 'TextColor', 'white', 'BoxColor', [255, 0, 0], 'BoxOpacity', 0.6);
        end

        concat = cat(2, concat, imgWithText);

        frameCount = frameCount + 1;

        fprintf('Frame #%d is loaded and written.\n', f);
        % Save concatenated image if 12 frames are reached
        if mod(frameCount/12, 1)==0
            outputFileName = fullfile(figSavePath, sprintf('%s_%d.png', strcat(vidName, '_', evtName), batchCount));
            imwrite(concat, outputFileName);
            fprintf('Concatenated image saved as %s\n', outputFileName);

            % Reset for next batch
            concat = [];
            frameCount = 0;
            batchCount = batchCount + 1;

        end
    end
end

% Save any remaining frames
if frameCount > 0
    outputFileName = fullfile(figSavePath, sprintf('%s_%d.png', strcat(vidName, '_', evtName), batchCount));
    imwrite(concat, outputFileName);
    fprintf('Concatenated image saved as %s\n', outputFileName);
end

end
