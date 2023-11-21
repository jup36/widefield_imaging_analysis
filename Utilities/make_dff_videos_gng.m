function make_dff_videos_gng(filePath, frameRate, pethTime, varargin)
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

%% sort frames to be included in the vide
cropFrames = @(a) a >= min(pethTime) & a <= max(pethTime);
pethFrameI = cellfun(cropFrames, {tbytDat(:).frameTrel}, 'UniformOutput', false);
[tbytDat.pethFrameI] = deal(pethFrameI{:});

%%
for t = 1:length(filePathTrials) % trials

    [~, vidName] = fileparts(filePathTrials{t});

    % Extract the number after 'trial_'
    matchedNum = regexp(vidName, 'trial_(\d+)', 'tokens');
    % Convert the matched number from cell to double
    trialNum = str2double(matchedNum{1}{1});

    fI = tbytDat(trialNum).pethFrameI; % frame logic for

    frameTrel = tbytDat(trialNum).frameTrel(fI); % get frame times relative to evt
    frameLickI = tbytDat(trialNum).frameLickI(fI);
    frameWaterI = tbytDat(trialNum).frameWaterI(fI);
    frameAirpuffI = tbytDat(trialNum).frameAirpuffI(fI);

    writeDffVidTrial(filePathTrials{t}, vidName, frameRate, fI, frameTrel, frameLickI, frameWaterI, frameAirpuffI)

    fprintf('Video #%d of total %d videos is written.\n', t, length(filePathTrials));
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function writeDffVidTrial(pngPath, vidName, frameRate, pethFrameI, frameTrel, frameLickI, frameWaterI, frameAirpuffI)
        pngs = GrabFiles_sort_trials('frame_', 0, {pngPath});

        assert(length(pngs)==length(pethFrameI))

        labeledVid = VideoWriter(fullfile(pngPath, [vidName, sprintf('_fr%d', frameRate), '.mp4']), 'MPEG-4');
        labeledVid.FrameRate = frameRate;
        playspeed = labeledVid.FrameRate/15;
        open(labeledVid);

        for f = 1:length(pngs) % frames
            if pethFrameI(f)
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

                if frameAirpuffI(f)
                    imgWithText = insertText(imgWithText, [10, 20], 'Airpuff', 'FontSize', 22, 'TextColor', 'white');
                end

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

