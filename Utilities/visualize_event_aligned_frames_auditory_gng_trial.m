function visualize_event_aligned_frames_auditory_gng_trial(filePath, numPeriFrames, trial, varargin)
% pethTime: peri-event time of the frames to be included in the video
% numPeriFrames: the # of frames to be included before and after the
%   specified event (e.g., hitLick). For example, if it's 40 it will take
%   40 frames before and 40 frames after the event.

%p = parseInput_visualize_event_aligned_frames(filePath, 10, varargin);
p = parseInput_visualize_event_aligned_frames(filePath, 10, trial, {'eventToAlign', 'hitLick'});

%% load tbytDat
[~, header] = fileparts(filePath);
fileBeh = GrabFiles_sort_trials('tbytDat_dff', 0, {fullfile(filePath, 'Matfiles')});
if isempty(fileBeh{1})
    fileBeh = GrabFiles_sort_trials('tbytDat.mat', 1, {filePath});
end
load(fullfile(fileBeh{1}), 'tbytDat')

fileParseGng = GrabFiles_sort_trials('tbytDat_parseGng', 0, {fullfile(filePath, 'Matfiles')});
tbytDatParseGng = load(fullfile(fileParseGng{1}), 'tbytDat');
tbytDatParseGng = tbytDatParseGng.('tbytDat');

filePathTrials =  GrabFiles_sort_trials('trial', 0, {fullfile(filePath, strcat(header, '_trials'))});

%% sort trials to work with
% To match the trial number use this!
tN = cellfun(@(x) str2double(regexp(x, 'trial(\d{1,3})', 'tokens', 'once')), filePathTrials);
filePathTrialC = filePathTrials(ismember(tN, trial));
filePathTrial = filePathTrialC{1};

if isempty(filePathTrials)
    warning('There is no trials ready to be filmed!')
end

%% sort frames to be included in the video
if isfield(tbytDat, 'resampledVidFrameTs')
    rsfrPeriod = mean(diff(tbytDat(trial).resampledVidFrameTs)); % resampled frame rate (50Hz downsampled from 200Hz)
    %frameResampleRate = 1/rsfrPeriod;
else
    error("resampled video timestamps are missing!!")
end

%% draw frames
[~, vidName] = fileparts(filePathTrial);

% Extract the number after 'trial_'
matchedNum = regexp(vidName, 'trial(\d+)', 'tokens');
% Convert the matched number from cell to double
trialNum = str2double(matchedNum{1}{1});

fT = tbytDat(trialNum).resampledVidFrameTs; % resampled frame time (likely downsampled)

switch p.Results.eventToAlign
    case 'hitLick'
        evtTs = tbytDatParseGng(trial).hitLicks;
        pethTime = evtTs + [-(rsfrPeriod*numPeriFrames) rsfrPeriod*numPeriFrames];
        cropFrames = @(a) a >= min(pethTime) & a <= max(pethTime);
end

fTrel = fT-tbytDat(trialNum).evtOn; % resampled frame time relative to event onset (cue onset)

% map temporal events to cmosExp pulses
tbytDat(trialNum).faceFrameLickI = check_timestamp_overlap_adjacentFrames(fTrel, tbytDat(trialNum).Lick-tbytDat(trialNum).evtOn, 3);
tbytDat(trialNum).faceFrameWaterI = check_timestamp_overlap_adjacentFrames(fTrel, tbytDat(trialNum).water-tbytDat(trialNum).evtOn, 3);
tbytDat(trialNum).faceFrameAirpuffI = check_timestamp_overlap_adjacentFrames(fTrel, tbytDat(trialNum).airpuff-tbytDat(trialNum).evtOn, 3);

tbytDat(trialNum).pethFrameI = cropFrames(fTrel); % mark relevant frames
fI = tbytDat(trialNum).pethFrameI; % frame logic for

frameTrel = fTrel(fI)-evtTs; % get frame times relative to evt
frameLickI = tbytDat(trialNum).faceFrameLickI(fI);
frameWaterI = tbytDat(trialNum).faceFrameWaterI(fI);
frameAirpuffI = tbytDat(trialNum).faceFrameAirpuffI(fI);

% draw peri hitLick frames
drawEvtAlignedFaceVidFramesAuditoryGngTenFrameSeq(filePathTrial, vidName, 'hitLick', 2, fI, frameTrel, frameLickI, frameWaterI, frameAirpuffI);

% draw all frames with cue-aligned timestamps
drawEvtAlignedFaceVidFramesAuditoryGngTenFrameSeq(filePathTrial, vidName, 'cueOn', 5, true(length(fTrel), 1), fTrel, ...,
    tbytDat(trialNum).faceFrameLickI, tbytDat(trialNum).faceFrameWaterI, tbytDat(trialNum).faceFrameAirpuffI);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function p = parseInput_visualize_event_aligned_frames(filePath, numPeriFrames, trial, vargs)
        % parse input, and extract name-value pairs
        default_eventToAlign = 'hitLick'; %

        p = inputParser; % create parser object
        addRequired(p, 'filePath');
        addRequired(p, 'numPeriFrames');
        addRequired(p, 'trial');

        addParameter(p, 'eventToAlign', default_eventToAlign);

        parse(p, filePath, numPeriFrames, trial, vargs{:})
    end

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


end

