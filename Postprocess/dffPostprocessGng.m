function dffPostprocessGng(filePath)
%filePath = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA008/DA008_101823';
%tbytDat
%   evtType:
%       1 or 2: visual (common, uncommon)
%           evtOn: photoDiodeOn
%           evtOff: photoDiodeOff
%           periEvtWin: evtOn-1:evtOff+1 (e.g. 4s)
%       3: reward
%           evtOn: water delivery
%           evtOff: 4-s after water delivery
%           periEvtWin: evtOn-1:evtOff (e.g. 5s)
%           Note that 5-s peri-reward window was used -1 to 4s relative to reward

% load trial-by-trial behavior & task data tbytDat.mat
%parentDir = fileparts(filePath);
[~, header] = fileparts(filePath);
fileBeh = GrabFiles_sort_trials('tbytDat', 0, {fullfile(filePath, 'Matfiles')});
if isempty(fileBeh{1})
    fileBeh = GrabFiles_sort_trials('tbytDat.mat', 1, {filePath});
end
load(fullfile(fileBeh{1}), 'tbytDat')

fileStim = GrabFiles_sort_trials('stimInfo', 0, {fullfile(filePath, 'Matfiles')});
load(fileStim{1}, 'stimopts')

% load the preprocessed dffs
%file_list_dff = GrabFiles_sort_trials('dff_combined.mat', 1, {fullfile(filePath, strcat(header, '_img'))}); % use GrabFiles_sort_trials to sort both files and folders

file_list_img_trial = findFoldersWithString(fullfile(filePath, strcat(header, '_img')), header);
file_list_dff = cell(1, length(file_list_img_trial)); 
for ii = 1:length(file_list_img_trial)
    file_list_dff{ii} = findFileWithString(file_list_img_trial{ii}, 'dff_combined'); 
    fprintf("Located dff file #%d\n", ii)
end

valDffFileNumb = sum(~cell2mat(cellfun(@isempty, file_list_dff, 'UniformOutput', false))); 

assert(valDffFileNumb==length(tbytDat))
assert(length(stimopts.rewarded_stim)==length(tbytDat))
[tbytDat(:).rewardTrI] = deal(0); 
[tbytDat(:).punishTrI] = deal(0); 

% file directory for trials
filePathTrials = fullfile(filePath, strcat(header, '_trials'));
if exist(filePathTrials, 'dir') == 0
    mkdir(filePathTrials);
end

%refCmosFrameIdx = 1:2700; % there must be 2700 cmos exposure pulses / frames recorded

dffCell = cell(1, length(file_list_dff));
dffsmCell = cell(1, length(file_list_dff)); 

% take corresponding frames with for each trial with 2D gaussian filtering
for tt = 1:length(tbytDat)
    if ~isempty(tbytDat(tt).cmos)
        trSubDir = fullfile(filePathTrials, sprintf('trial_%d', tt));
        if exist(trSubDir, 'dir') ~= 7
            mkdir(trSubDir);
        end

        if stimopts.rewarded_stim(tt)==1
            tbytDat(tt).rewardTrI = 1; 
        elseif stimopts.punished_stim(tt)==1
            tbytDat(tt).punishTrI = 1; 
        end

        load(file_list_dff{tt}, 'dff'); % load dff

        tbytDat(tt).frameT = tbytDat(tt).cmos(1:2:end); % frame time (needs to alternate due to interleaved violet frames for hemodynamic correction)
        tbytDat(tt).dff = dff(:,:,1:min(length(tbytDat(tt).frameT), size(dff, 3))); % aligned dff
        tbytDat(tt).dffsm = applyImgaussfilt(tbytDat(tt).dff);

        % map temporal events to cmosExp pulses
        tbytDat(tt).frameLickI = check_timestamp_overlap(tbytDat(tt).frameT, tbytDat(tt).Lick);
        tbytDat(tt).frameWaterI = check_timestamp_overlap(tbytDat(tt).frameT, tbytDat(tt).water);
        tbytDat(tt).frameAirpuffI = check_timestamp_overlap(tbytDat(tt).frameT, tbytDat(tt).airpuff); 

        % timestamp each frame relative to the stim onset
        tbytDat(tt).frameStimOnI = check_timestamp_overlap(tbytDat(tt).frameT, tbytDat(tt).stimOn);
        tbytDat(tt).frameStimOffI = check_timestamp_overlap(tbytDat(tt).frameT, tbytDat(tt).stimOff);
        tbytDat(tt).frameStimI = zeros(length(tbytDat(tt).frameT), 1);
        tbytDat(tt).frameStimI(find(tbytDat(tt).frameStimOnI, 1):find(tbytDat(tt).frameStimOffI, 1), 1) = 1;

        tbytDat(tt).frameTrel = tbytDat(tt).frameT-tbytDat(tt).stimOn;
        tbytDat(tt).faceCamTrel = tbytDat(tt).faceCam-tbytDat(tt).stimOn; 

        dffCell{tt} =  tbytDat(tt).dff;
        dffsmCell{tt} =  tbytDat(tt).dffsm;
        
        
        tbytDff =  tbytDat(tt).dff;
        tbytDffsm =  tbytDat(tt).dffsm;
        save(fullfile(trSubDir, 'tbytDff.mat'), 'tbytDff', 'tbytDffsm')

    end
    fprintf('processed dff of trial#%d\n', tt)
end

% record frames save in separate folders
%record_dff_frames(tbytDat, filePathTrials)

% make dff videos
%frameRate = 5;
%make_dff_videos(tbytDat, filePathTrials, frameRate, 37)

% save tbytDat without dffs
tbytDat = rmfield(tbytDat, {'dff', 'dffsm'});
save(fullfile(filePath, 'Matfiles', [header, '_tbytDat_dff']), 'tbytDat')
save(fullfile(filePath, 'Matfiles', [header, '_dffCollect']), 'dffCell', '-v7.3')
save(fullfile(filePath, 'Matfiles', [header, '_dffsmCollect']), 'dffsmCell', '-v7.3')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function record_dff_frames(tbytDat, filePathTrials, varargin)
        % e.g., dff = tbytDat(21).dffsm
        %By default this function iterate through all trials and record frames of
        % each trial as png files. It can also run specific trials designated as varargin.
        if isempty(varargin)
            trials = 1:length(tbytDat);
        else
            trials = [varargin{1}(:)]';
        end

        %% record frames first
        close all;
        for t = trials % 1:length(tbytDat) % trial
            if ~isempty(tbytDat(t).dffsm)
                pngSubDir = fullfile(filePathTrials, sprintf('block_%d_trial_%d', tbytDat(t).cmosExpTrainI, t));

                if exist(pngSubDir, 'dir') == 7

                else
                    mkdir(pngSubDir);
                    for i = 1:size(tbytDat(t).dffsm,3)
                        pngName = sprintf('frame_%d.png', i);
                        figHandle = imageFrameWithNaNs(tbytDat(t).dffsm(:, :, i), [-2 2]); hold on;

                        % Turn off the axes
                        axis off;
                        % Set the figure's background to white
                        set(gcf, 'Color', 'w');

                        if isfield(tbytDat, 'frameStimI')
                            if ~isempty(tbytDat(t).frameStimI) && tbytDat(t).frameStimI(i) % for visual trials
                                if tbytDat(t).evtType == 1 % common visual stim (draw 45 degree lines at the upper left corner)
                                    insertgrating45(figHandle, tbytDat(t).dffsm(:, :, i))
                                elseif tbytDat(t).evtType == 2 % common visual stim (draw 135 degree lines at the upper left corner)
                                    insertgrating135(figHandle, tbytDat(t).dffsm(:, :, i))
                                end
                                hold off;
                            end
                        end
                        % Capture the current figure with dots
                        frameLabled = getframe(gca);
                        frameLabled = frameLabled.cdata;

                        % save the figure with dots

                        imwrite(frameLabled, fullfile(pngSubDir, pngName));

                        fprintf('Frame #%d of trial #%d is labeled and saved.\n', i, t);
                        close all
                    end
                end
            end
        end

    end


    function make_dff_videos(tbytDat, filePathTrials, frameRate)
        close all;
        filePathTrials = GrabFiles_sort_trials('trial', 0, {filePathTrials}); % use GrabFiles_sort_trials to sort both files and folders

        % To match the trial number use this!
        %tN = cellfun(@(x) str2double(regexp(x, 'trial_(\d{1,3})', 'tokens', 'once')), filePathTrials);
        %t = find(tN == 622);

        for t = 1:length(filePathTrials) % trials

            [~, vidName] = fileparts(filePathTrials{t});

            % Extract the number after 'trial_'
            matchedNum = regexp(vidName, 'trial_(\d+)', 'tokens');
            % Convert the matched number from cell to double
            trialNum = str2double(matchedNum{1}{1});
            frameTrel = tbytDat(trialNum).frameTrel; % get frame times relative to evt
            frameLickI = tbytDat(trialNum).frameLickI;
            frameWaterI = tbytDat(trialNum).frameWaterI;

            writeDffVid(filePathTrials{t}, vidName, frameRate, frameTrel, frameLickI, frameWaterI)

        end


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
    end

    function correctTimeIdx = correctFrameI(refArray, timeIdx1, numPoints)
        corrArray = 1:0.5:floor(max(refArray)/2);
        correctTimeIdx1 = floor(corrArray(timeIdx1));
        correctTimeIdx2 = correctTimeIdx1+numPoints-1;
        correctTimeIdx = correctTimeIdx1:correctTimeIdx2;
    end

    function smFrames = applyImgaussfilt(frames)
        smFrames = zeros(size(frames));
        % Loop through each frame and apply smoothing
        for frame = 1:size(frames, 3)
            % Extract the current frame
            current_frame = frames(:, :, frame);

            % Create a binary mask of NaN values
            nan_mask = isnan(current_frame);

            % Temporarily replace NaNs with the mean of non-NaN values in the frame
            current_frame(nan_mask) = nanmean(current_frame(:));

            % Apply Gaussian smoothing
            smFrame = imgaussfilt(current_frame, 1); % Change 1 to desired sigma if needed

            % Restore NaN values using the mask
            smFrame(nan_mask) = NaN;

            % Store the smoothed frame
            smFrames(:,:,frame) = smFrame;
        end
    end

    function figHandle = imageFrameWithNaNs(frame, climits)
        % climits = [-3 3];
        % Create a new figure and store its handle
        figHandle = figure;

        % Create an adjusted colormap with white at the beginning
        cmap = [1, 1, 1; hot(256)];

        % Display the frame
        imagesc(frame);

        % Use the custom colormap
        colormap(cmap);

        % Set NaN values to a value outside the data range for visualization
        clim(climits);

        % Show the colorbar (optional)
        colorbar;
    end

% Example usage:
% frame = rand(68, 68);
% frame(randi([1, 68], 10), randi([1, 68], 10)) = NaN;
% fig = plotFrameWithNans(frame);

    function refTsI = check_timestamp_overlap(refTs, evtTs)
        % Generate example data (replace with your actual data)
        %refTs = faceCam_pt;
        %evtTs = xStampSec_lick;

        % Find the minimum and maximum timestamps of the first set
        minRefTs = min(refTs);
        maxRefTs = max(refTs);

        refTsI = zeros(length(refTs), 1);

        % Select second timestamps that fall between the minimum and maximum of the first set
        selectEvtTs = evtTs(evtTs(:, 1) >= minRefTs & evtTs(:, 1) <= maxRefTs, 1);

        if ~isempty(selectEvtTs)
            for ts = 1:length(selectEvtTs)
                [~, minI] = min(abs(refTs - selectEvtTs(ts, 1)));
                refTsI(minI) = 1;
            end
        end

    end


end




