function writeSnapshot(pngPath, vidName, frameIndices, maxImageN, frameTrel, frameLickI, frameWaterI, frameAirpuffI)
pngs_all = GrabFiles_sort_trials('faceDffGR_frame_', 0, {pngPath});

numFrames = length(frameIndices);
numImages = min(numFrames, maxImageN);
if numImages == 0
    error('No images found for snapshot.');
end

% Select evenly distributed frames
selectedIndices = round(linspace(1, numFrames, numImages));
pngs = pngs_all(frameIndices(selectedIndices));

images = cell(1, numImages);
for i = 1:numImages
    img = imread(pngs{i});

    % Add frame labels
    frameTimeStr = sprintf('%.2fs', frameTrel(selectedIndices(i)));
    img = insertText(img, [10, 10], frameTimeStr, 'FontSize', 22, 'TextColor', 'white', 'BoxColor', [255, 0, 255], 'BoxOpacity', 0.3);

    if frameLickI(selectedIndices(i))
        img = insertText(img, [10, 40], 'Lick', 'FontSize', 22, 'TextColor', 'white');
    end

    if frameWaterI(selectedIndices(i))
        img = insertText(img, [10, 70], 'Water', 'FontSize', 22, 'TextColor', 'white', 'BoxColor', [0, 255, 255], 'BoxOpacity', 0.6);
    elseif frameAirpuffI(selectedIndices(i))
        img = insertText(img, [10, 70], 'Airpuff', 'FontSize', 22, 'TextColor', 'white', 'BoxColor', [255, 0, 0], 'BoxOpacity', 0.6);
    end

    images{i} = img;
end

[imgHeight, imgWidth, imgChannels] = size(images{1});
gridCols = ceil(sqrt(numImages));
gridRows = ceil(numImages / gridCols);
margin = 10; % Define margin for readability

snapshot = uint8(255 * ones(gridRows * (imgHeight + margin) - margin, gridCols * (imgWidth + margin) - margin, imgChannels));

for idx = 1:numImages
    row = floor((idx-1) / gridCols);
    col = mod((idx-1), gridCols);

    yStart = row * (imgHeight + margin) + 1;
    yEnd = yStart + imgHeight - 1;
    xStart = col * (imgWidth + margin) + 1;
    xEnd = xStart + imgWidth - 1;

    snapshot(yStart:yEnd, xStart:xEnd, :) = images{idx};
end

figure; imshow(snapshot)
imwrite(snapshot, fullfile(pngPath, [vidName, '_snapshot.png']));
fprintf('Snapshot saved as %s\n', fullfile(pngPath, [vidName, '_snapshot.png']));
end