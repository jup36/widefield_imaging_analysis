function saveEachFrameAsImage(pngPath, vidName, frameRate, pethFrameI, frameTrel, frameLickI, frameWaterI, frameAirpuffI)
pngs = GrabFiles_sort_trials('frame_', 0, {pngPath});

assert(length(pngs)==length(pethFrameI))

labeledVid = VideoWriter(fullfile(pngPath, [vidName, sprintf('_fr%d', frameRate), '.mp4']), 'MPEG-4');
labeledVid.FrameRate = frameRate;
playspeed = labeledVid.FrameRate/15;
open(labeledVid);

frameCount = 0;

for f = 1:length(pngs) % frames
    if pethFrameI(f)
        frameCount = frameCount+1;

        % Load png image
        img = imread(pngs{f});
        [imgHeight, imgWidth, ~] = size(img); % Get the size of the image

        % Define margins and estimate text box height based on font size
        margin = 20; % Margin from the edge in pixels
        textBoxHeight = 20; % Estimate text box height (may need adjustment based on font size)
        textBoxWidth = 70; % Estimate text box width (may need adjustment based on font size and text length)

        % Calculate positions for the texts
        playspeedPos = [margin, imgHeight - textBoxHeight - margin]; % Bottom left
        frameTimePos = [imgWidth - textBoxWidth - margin, imgHeight - textBoxHeight - margin]; % Bottom right
        lickPos = [margin, margin]; % Top left
        waterPos = [margin, 2 * margin + textBoxHeight]; % Below the 'Lick' position
        % For Airpuff, if it's exclusive with Water, use waterPos, otherwise calculate a new position

        % Add playspeed to the image
        playSpeedStr = sprintf('playspeed x%.1f', playspeed);
        imgWithText = insertText(img, playspeedPos, playSpeedStr, 'FontSize', 22, 'TextColor', 'white', 'BoxColor', [1, 0, 1], 'BoxOpacity', 0.8);

        % Add frametime to the image
        frameTimeStr = sprintf('%.2fs', frameTrel(frameCount));
        imgWithText = insertText(imgWithText, frameTimePos, frameTimeStr, 'FontSize', 22, 'TextColor', 'white');

        % Add 'Lick' to the image if applicable
        if frameLickI(frameCount)
            imgWithText = insertText(imgWithText, lickPos, 'Lick', 'FontSize', 22, 'TextColor', 'white');
        end

        % Add 'Water' or 'Airpuff' to the image if applicable
        if frameWaterI(frameCount)
            imgWithText = insertText(imgWithText, waterPos, 'Water', 'FontSize', 22, 'TextColor', 'white');
        elseif frameAirpuffI(frameCount) % Assuming Water and Airpuff don't occur together
            imgWithText = insertText(imgWithText, waterPos, 'Airpuff', 'FontSize', 22, 'TextColor', 'white');
        end

        % Display the image
        imshow(imgWithText);
    
        % Save the displayed figure as a PDF with higher resolution
        saveFileName = fullfile(pngPath, sprintf('FrameWithText_%d_%d.png', f, frameCount));
        print(gcf, saveFileName, '-dpng', '-r300'); % '-r300' sets the resolution to 300 dpi

        frame = im2frame(imgWithText);

        % write the labeled frame to the output video
        writeVideo(labeledVid, frame);
        %fprintf('Frame #%d is loaded and written.\n', f);
    end
end
close(labeledVid);
end
