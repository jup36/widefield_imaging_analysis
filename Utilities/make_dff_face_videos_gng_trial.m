function make_dff_face_videos_gng_trial(filePath, frameRate, pethTime, varargin)
% pethTime: peri-event time of the frames to be included in the video

%% load tbytDat
[~, header] = fileparts(filePath);
fileBeh = GrabFiles_sort_trials('tbytDat_dff', 0, {fullfile(filePath, 'Matfiles')});
if isempty(fileBeh{1})
    fileBeh = GrabFiles_sort_trials('tbytDat.mat', 1, {filePath});
end
load(fullfile(fileBeh{1}), 'tbytDat')

filePathTrials =  GrabFiles_sort_trials('trial', 0, {fullfile(filePath, strcat(header, '_trials'))});

%% sort trials to work with
if isempty(varargin)
    trials = 1:length(filePathTrials);
else
    % To match the trial number use this!
    tN = cellfun(@(x) str2double(regexp(x, 'trial_(\d{1,3})', 'tokens', 'once')), filePathTrials);
    trials = [varargin{1}(:)]';
    filePathTrials = filePathTrials(ismember(tN, trials));
end

if isempty(filePathTrials)
    warning('There is no trials ready to be filmed!')
end

%% sort frames to be included in the video
cropFrames = @(a) a >= min(pethTime) & a <= max(pethTime);
pethFrameI = cellfun(cropFrames, {tbytDat(:).faceCamTrelFrames}, 'UniformOutput', false);
[tbytDat.pethFrameI] = deal(pethFrameI{:});

%%
for t = 1:length(filePathTrials) % trials

    [~, vidName] = fileparts(filePathTrials{t});

    % Extract the number after 'trial_'
    matchedNum = regexp(vidName, 'trial_(\d+)', 'tokens');
    % Convert the matched number from cell to double
    trialNum = str2double(matchedNum{1}{1});

    fT = tbytDat(trialNum).faceCamTrelFrames; 
    
    % map temporal events to cmosExp pulses
    tbytDat(trialNum).faceFrameLickI = check_timestamp_overlap_adjacentFrames(fT, tbytDat(trialNum).Lick-tbytDat(trialNum).stimOn, 5);
    tbytDat(trialNum).faceFrameWaterI = check_timestamp_overlap_adjacentFrames(fT, tbytDat(trialNum).water-tbytDat(trialNum).stimOn, 5);
    tbytDat(trialNum).faceFrameAirpuffI = check_timestamp_overlap_adjacentFrames(fT, tbytDat(trialNum).airpuff-tbytDat(trialNum).stimOn, 5); 


    fI = tbytDat(trialNum).pethFrameI; % frame logic for

    frameTrel = fT(fI); % get frame times relative to evt
    frameLickI = tbytDat(trialNum).faceFrameLickI(fI);
    frameWaterI = tbytDat(trialNum).faceFrameWaterI(fI);
    frameAirpuffI = tbytDat(trialNum).faceFrameAirpuffI(fI);

    writeDffFaceVidTrial(filePathTrials{t}, vidName, frameRate, fI, frameTrel, frameLickI, frameWaterI, frameAirpuffI)

    fprintf('Video #%d of total %d videos is writen.\n', t, length(filePathTrials));
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeDffFaceVidTrial(pngPath, vidName, frameRate, pethFrameI, frameTrel, frameLickI, frameWaterI, frameAirpuffI)
pngs_all = GrabFiles_sort_trials('faceDff_frame_', 0, {pngPath});
pngs = pngs_all(pethFrameI); 

assert(length(pngs)==sum(pethFrameI))

labeledVid = VideoWriter(fullfile(pngPath, [vidName, sprintf('_fr%d', frameRate), '.mp4']), 'MPEG-4');
labeledVid.FrameRate = frameRate;
playspeed = labeledVid.FrameRate/200; % face video frame rate
open(labeledVid);

for f = 1:length(pngs) % frames
    if pethFrameI(f)
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
        imgWithText = insertText(img, playspeedPos, playSpeedStr, 'FontSize', 22, 'TextColor', 'white', 'BoxColor', [199, 234, 70], 'BoxOpacity', 0.8);

        % Add frametime to the image
        frameTimeStr = sprintf('%.2fs', frameTrel(f));
        imgWithText = insertText(imgWithText, frameTimePos, frameTimeStr, 'FontSize', 22, 'TextColor', 'white', 'BoxColor', [255, 0, 255], 'BoxOpacity', 0.3);

        % Add 'Lick' to the image if applicable
        if frameLickI(f)
            imgWithText = insertText(imgWithText, lickPos, 'Lick', 'FontSize', 22, 'TextColor', 'white');
        end

        % Add 'Water' or 'Airpuff' to the image if applicable
        if frameWaterI(f)
            imgWithText = insertText(imgWithText, waterPos, 'Water', 'FontSize', 22, 'TextColor', 'white',  'BoxColor', [0, 255, 255], 'BoxOpacity', 0.6);
        elseif frameAirpuffI(f) % Assuming Water and Airpuff don't occur together
            imgWithText = insertText(imgWithText, waterPos, 'Airpuff', 'FontSize', 20, 'TextColor', 'white',  'BoxColor', [255, 0, 0], 'BoxOpacity', 0.6);
        end

        % Display the image
        imshow(imgWithText);

        frame = im2frame(imgWithText);

        % write the labeled frame to the output video
        writeVideo(labeledVid, frame);
        %fprintf('Frame #%d is loaded and written.\n', f);
    end
end
close(labeledVid);
end

end

