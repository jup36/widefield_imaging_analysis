function writeDffVid(pngPath, vidName, frameRate, frameTrel, frameLickI, frameWaterI)
pngs = GrabFiles_sort_trials('frame_', 0, {pngPath});
labeledVid = VideoWriter(fullfile(pngPath, [vidName, sprintf('_fr%d', frameRate), '.mp4']), 'MPEG-4');
labeledVid.FrameRate = frameRate;
playspeed = labeledVid.FrameRate/15;
open(labeledVid);

for f = 1:length(pngs) % frames
    % load png image
    img = imread(pngs{f});

    % Add playspeed to the image
    playSpeedStr = sprintf('playspeed x%.1f', playspeed);
    imgWithText = insertText(img, [10, 607], playSpeedStr, 'FontSize', 22, 'TextColor', 'white', 'BoxColor', [1, 0, 1], 'BoxOpacity', 0.8);

    % Add frametime to the image
    frameTimeStr = sprintf('%.2fs', frameTrel(f));
    imgWithText = insertText(imgWithText, [700, 607], frameTimeStr, 'FontSize', 22, 'TextColor', 'white');

    if frameLickI(f)
        imgWithText = insertText(imgWithText, [10, 60], 'Lick', 'FontSize', 22, 'TextColor', 'white');
    end

    if frameWaterI(f)
        imgWithText = insertText(imgWithText, [10, 20], 'Water', 'FontSize', 22, 'TextColor', 'white');
    end

    imshow(imgWithText);

    frame = im2frame(imgWithText);

    % write the labeled frame to the output video
    writeVideo(labeledVid, frame);
    fprintf('Frame #%d is loaded and written.\n', f);
end
close(labeledVid);
end