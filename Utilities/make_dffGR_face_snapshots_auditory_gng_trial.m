function make_dffGR_face_snapshots_auditory_gng_trial(filePath, pethTime, maxImageN, trialNum)
% pethTime: peri-event time of the frames to be included in the snapshot

%% load tbytDat
match_header = regexp(filePath, 'm\d{1,4}_\d{6}', 'match');
header = match_header{1};

fileBehG = GrabFiles_sort_trials('green_tbytDat_dff', 0, {fullfile(filePath, 'Matfiles')});
if isempty(fileBehG{1})
    fileBehG = GrabFiles_sort_trials('tbytDat.mat', 1, {filePath});
end
load(fullfile(fileBehG{1}), 'tbytDatG')

filePathTrials = GrabFiles_sort_trials('_trials', 0, {filePath});

%% Process a single trial
filePath_faceDffGR = GrabFiles_sort_trials(['faceDffGR*', sprintf('trial_%d', trialNum)], 0, filePathTrials);
[~, vidName] = fileparts(filePath_faceDffGR{1});
pethStr = convert_pethTime_string(pethTime);
vidName = [vidName, '_', pethStr];

fT = tbytDatG(trialNum).resampledVidFrameTs;
fTrel = fT - tbytDatG(trialNum).evtOn;

% map temporal events to cmosExp pulses
tbytDatG(trialNum).faceFrameLickI = check_timestamp_overlap_adjacentFrames(fTrel, tbytDatG(trialNum).Lick-tbytDatG(trialNum).evtOn, 3);
tbytDatG(trialNum).faceFrameWaterI = check_timestamp_overlap_adjacentFrames(fTrel, tbytDatG(trialNum).water-tbytDatG(trialNum).evtOn, 3);
tbytDatG(trialNum).faceFrameAirpuffI = check_timestamp_overlap_adjacentFrames(fTrel, tbytDatG(trialNum).airpuff-tbytDatG(trialNum).evtOn, 3);

tbytDatG(trialNum).pethFrameI = (fTrel >= min(pethTime)) & (fTrel <= max(pethTime));
fI = find(tbytDatG(trialNum).pethFrameI);
fTrel = fTrel(fI); 

frameLickI = tbytDatG(trialNum).faceFrameLickI(fI);
frameWaterI = tbytDatG(trialNum).faceFrameWaterI(fI);
frameAirpuffI = tbytDatG(trialNum).faceFrameAirpuffI(fI);

writeSnapshot(filePath_faceDffGR{1}, vidName, fI, maxImageN, fTrel, frameLickI, frameWaterI, frameAirpuffI);

fprintf('Snapshot for trial %d is saved.\n', trialNum);
end

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


% function make_dffGR_face_snapshots_auditory_gng_trial(filePath, pethTime, maxImageN, trialNum)
% % pethTime: peri-event time of the frames to be included in the snapshot
% 
% %% load tbytDat
% match_header = regexp(filePath, 'm\d{1,4}_\d{6}', 'match');
% header = match_header{1};
% 
% fileBehG = GrabFiles_sort_trials('green_tbytDat_dff', 0, {fullfile(filePath, 'Matfiles')});
% if isempty(fileBehG{1})
%     fileBehG = GrabFiles_sort_trials('tbytDat.mat', 1, {filePath});
% end
% load(fullfile(fileBehG{1}), 'tbytDatG')
% 
% filePathTrials = GrabFiles_sort_trials('_trials', 0, {filePath});
% 
% %% Process a single trial
% filePath_faceDffGR = GrabFiles_sort_trials(['faceDffGR*', sprintf('trial_%d', trialNum)], 0, filePathTrials);
% [~, vidName] = fileparts(filePath_faceDffGR{1});
% pethStr = convert_pethTime_string(pethTime); 
% vidName = [vidName, '_', pethStr]; 
% 
% fT = tbytDatG(trialNum).resampledVidFrameTs;
% fTrel = fT - tbytDatG(trialNum).evtOn;
% 
% % map temporal events to cmosExp pulses
% tbytDatG(trialNum).faceFrameLickI = check_timestamp_overlap_adjacentFrames(fTrel, tbytDatG(trialNum).Lick-tbytDatG(trialNum).evtOn, 3);
% tbytDatG(trialNum).faceFrameWaterI = check_timestamp_overlap_adjacentFrames(fTrel, tbytDatG(trialNum).water-tbytDatG(trialNum).evtOn, 3);
% tbytDatG(trialNum).faceFrameAirpuffI = check_timestamp_overlap_adjacentFrames(fTrel, tbytDatG(trialNum).airpuff-tbytDatG(trialNum).evtOn, 3);
% 
% tbytDatG(trialNum).pethFrameI = (fTrel >= min(pethTime)) & (fTrel <= max(pethTime));
% fI = find(tbytDatG(trialNum).pethFrameI);
% 
% frameLickI = tbytDatG(trialNum).faceFrameLickI(fI);
% frameWaterI = tbytDatG(trialNum).faceFrameWaterI(fI);
% frameAirpuffI = tbytDatG(trialNum).faceFrameAirpuffI(fI);
% 
% writeSnapshot(filePath_faceDffGR{1}, vidName, fI, maxImageN, fTrel, frameLickI, frameWaterI, frameAirpuffI);
% 
% fprintf('Snapshot for trial %d is saved.\n', trialNum);
% end
% 
% function writeSnapshot(pngPath, vidName, frameIndices, maxImageN, frameTrel, frameLickI, frameWaterI, frameAirpuffI)
% pngs_all = GrabFiles_sort_trials('faceDffGR_frame_', 0, {pngPath});
% 
% numFrames = length(frameIndices);
% numImages = min(numFrames, maxImageN);
% if numImages == 0
%     error('No images found for snapshot.');
% end
% 
% % Select evenly distributed frames
% selectedIndices = round(linspace(1, numFrames, numImages));
% pngs = pngs_all(frameIndices(selectedIndices));
% 
% images = cell(1, numImages);
% for i = 1:numImages
%     img = imread(pngs{i});
% 
%     % Add frame labels
%     frameTimeStr = sprintf('%.2fs', frameTrel(selectedIndices(i)));
%     img = insertText(img, [10, 10], frameTimeStr, 'FontSize', 22, 'TextColor', 'white', 'BoxColor', [255, 0, 255], 'BoxOpacity', 0.3);
% 
%     if frameLickI(selectedIndices(i))
%         img = insertText(img, [10, 40], 'Lick', 'FontSize', 22, 'TextColor', 'white');
%     end
% 
%     if frameWaterI(selectedIndices(i))
%         img = insertText(img, [10, 70], 'Water', 'FontSize', 22, 'TextColor', 'white', 'BoxColor', [0, 255, 255], 'BoxOpacity', 0.6);
%     elseif frameAirpuffI(selectedIndices(i))
%         img = insertText(img, [10, 70], 'Airpuff', 'FontSize', 22, 'TextColor', 'white', 'BoxColor', [255, 0, 0], 'BoxOpacity', 0.6);
%     end
% 
%     images{i} = img;
% end
% 
% [imgHeight, imgWidth, imgChannels] = size(images{1});
% 
% gridCols = ceil(sqrt(numImages));
% gridRows = ceil(numImages / gridCols);
% 
% snapshot = uint8(255 * ones(gridRows * imgHeight, gridCols * imgWidth, imgChannels));
% 
% for idx = 1:numImages
%     row = floor((idx-1) / gridCols);
%     col = mod((idx-1), gridCols);
% 
%     yStart = row * imgHeight + 1;
%     yEnd = (row + 1) * imgHeight;
%     xStart = col * imgWidth + 1;
%     xEnd = (col + 1) * imgWidth;
% 
%     snapshot(yStart:yEnd, xStart:xEnd, :) = images{idx};
% end
% 
% figure; imshow(snapshot)
% imwrite(snapshot, fullfile(pngPath, [vidName, '_snapshot.png']));
% fprintf('Snapshot saved as %s\n', fullfile(pngPath, [vidName, '_snapshot.png']));
% end
