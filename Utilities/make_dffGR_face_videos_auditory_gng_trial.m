function make_dffGR_face_videos_auditory_gng_trial(filePath, frameRate, pethTime, varargin)
% pethTime: peri-event time of the frames to be included in the video

%% load tbytDat
match_header = regexp(filePath, 'm\d{1,4}_\d{6}', 'match');
header = match_header{1};

fileBehG = GrabFiles_sort_trials('green_tbytDat_dff', 0, {fullfile(filePath, 'Matfiles')});
if isempty(fileBehG{1})
    fileBehG = GrabFiles_sort_trials('tbytDat.mat', 1, {filePath});
end
load(fullfile(fileBehG{1}), 'tbytDatG')

filePathTrials = GrabFiles_sort_trials('_trials', 0, {filePath});

%% sort trials to work with
if isempty(varargin)
    trials = 1:length(tbytDatG);
else
    trials = [varargin{1}(:)]';
end

%% sort frames to be included in the video
cropFrames = @(a) a >= min(pethTime) & a <= max(pethTime);

%%
for t = trials % trials

    filePath_faceDffGR = GrabFiles_sort_trials(['faceDffGR*', sprintf('trial_%d', t)], 0, filePathTrials); 
    [~, vidName] = fileparts(filePath_faceDffGR{1});
    pethStr = convert_pethTime_string(pethTime);
    vidName = [vidName, '_', pethStr];

    % Extract the number after 'trial_'
    matchedNum = regexp(vidName, 'trial_(\d+)', 'tokens');
    % Convert the matched number from cell to double
    trialNum = str2double(matchedNum{1}{1});

    fT = tbytDatG(trialNum).resampledVidFrameTs; % resampled frame time (likely downsampled) 
    frameResampleRate = 1/mean(diff(fT)); 
    fTrel = fT-tbytDatG(trialNum).evtOn; % resampled frame time relative to event onset (cue onset)

    % map temporal events to cmosExp pulses
    tbytDatG(trialNum).faceFrameLickI = check_timestamp_overlap_adjacentFrames(fTrel, tbytDatG(trialNum).Lick-tbytDatG(trialNum).evtOn, 3);
    tbytDatG(trialNum).faceFrameWaterI = check_timestamp_overlap_adjacentFrames(fTrel, tbytDatG(trialNum).water-tbytDatG(trialNum).evtOn, 3);
    tbytDatG(trialNum).faceFrameAirpuffI = check_timestamp_overlap_adjacentFrames(fTrel, tbytDatG(trialNum).airpuff-tbytDatG(trialNum).evtOn, 3); 

    tbytDatG(trialNum).pethFrameI = cropFrames(fTrel);
    fI = tbytDatG(trialNum).pethFrameI; % frame logic for

    frameTrel = fTrel(fI); % get frame times relative to evt
    frameLickI = tbytDatG(trialNum).faceFrameLickI(fI);
    frameWaterI = tbytDatG(trialNum).faceFrameWaterI(fI);
    frameAirpuffI = tbytDatG(trialNum).faceFrameAirpuffI(fI);

    writeDualDffFaceVidTrialAuditoryGng(filePath_faceDffGR{1}, vidName, frameRate, frameResampleRate, fI, frameTrel, frameLickI, frameWaterI, frameAirpuffI)

    fprintf('Video #%d of total %d videos is writen.\n', t, length(filePathTrials));
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeDualDffFaceVidTrialAuditoryGng(pngPath, vidName, frameRate, frameResampleRate, pethFrameI, frameTrel, frameLickI, frameWaterI, frameAirpuffI)
pngs_all = GrabFiles_sort_trials('faceDffGR_frame_', 0, {pngPath});
pngs = pngs_all(pethFrameI); 

assert(length(pngs)==sum(pethFrameI))

labeledVid = VideoWriter(fullfile(pngPath, [vidName, sprintf('_fr%d', frameRate), '.mp4']), 'MPEG-4');
labeledVid.FrameRate = frameRate;
playspeed = frameRate/frameResampleRate; % face video frame rate
open(labeledVid);

for f = 1:length(pngs) % frames
%    if pethFrameI(f)
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
 %   end
end
close(labeledVid);
end


end

