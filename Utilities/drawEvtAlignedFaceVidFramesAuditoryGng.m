function drawEvtAlignedFaceVidFramesAuditoryGng(pngPath, vidName, pethFrameI, frameTrelEvt, frameLickI, frameWaterI, frameAirpuffI)
pngs_all = GrabFiles_sort_trials('faceDff_frame_', 0, {pngPath});
pngs = pngs_all(pethFrameI);

assert(length(pngs)==sum(pethFrameI))

% Define margins and estimate text box height based on font size
margin = 20; % Margin from the edge in pixels
textBoxHeight = 20; % Estimate text box height (may need adjustment based on font size)
textBoxWidth = 70; % Estimate text box width (may need adjustment based on font size and text length)

concat_preEvt = [];
concat_postEvt = [];

for f = 1:length(pngs) % frames
    if mod(f, 2) == 0
        % Load png image
        imgFull = imread(pngs{f});
        [~, imgWidthFull, ~] = size(imgFull); % Get the size of the image

        img = imgFull(:, 1:floor(imgWidthFull/2)-1, :); % crop the face video

        [imgHeight, imgWidth, ~] = size(img);

        %% label image
        frameTimePos = [imgWidth - textBoxWidth - margin, imgHeight - textBoxHeight - margin]; % Bottom right
        lickPos = [margin, margin]; % Top left
        waterPos = [margin, 2 * margin + textBoxHeight]; % Below the 'Lick' position

        % Add frametime to the image
        frameTimeStr = sprintf('%.2fs', frameTrelEvt(f));
        imgWithText = insertText(img, frameTimePos, frameTimeStr, 'FontSize', 22, 'TextColor', 'white', 'BoxColor', [255, 0, 255], 'BoxOpacity', 0.3);

        % Add lick for lick frames
        if frameLickI(f)
            imgWithText = insertText(imgWithText, lickPos, 'Lick', 'FontSize', 22, 'TextColor', 'white');

            % Draw a rectangle around the edge of the image
            rectPosition = [1, 1, imgWidth-1, imgHeight-1]; % Rectangle position [x, y, width, height]
            imgWithText = insertShape(imgWithText, 'Rectangle', rectPosition, 'Color', 'cyan', 'LineWidth', 5);
            %imshow(imgWithText)
        end

        % Add 'Water' or 'Airpuff' to the image if applicable
        if frameWaterI(f)
            imgWithText = insertText(imgWithText, waterPos, 'Water', 'FontSize', 22, 'TextColor', 'white',  'BoxColor', [0, 255, 255], 'BoxOpacity', 0.6);
        elseif frameAirpuffI(f) % Assuming Water and Airpuff don't occur together
            imgWithText = insertText(imgWithText, waterPos, 'Airpuff', 'FontSize', 20, 'TextColor', 'white',  'BoxColor', [255, 0, 0], 'BoxOpacity', 0.6);
        end

        % concatenate labeled image
        if frameTrelEvt(f) <= 0
            concat_preEvt = cat(2, concat_preEvt, imgWithText);
        else
            concat_postEvt = cat(2, concat_postEvt, imgWithText);
        end

        fprintf('Frame #%d is loaded and written.\n', f);

    end
end
% Display the image
imshow(concat_postEvt);

end
