function dffPostprocess_auditory_gng_dual(filePath, varargin)
%This function
%
% filePath = 'Z:\Rodent Data\dualImaging_parkj\m1045_jRGECO_GRABda\m1045_121324\task';
%tbytDat
%   evtType:
%       1 or 2: auditory (go, no-go)
%           evtOn:
%           evtOff: photoDiodeOff
%           periEvtWin: evtOn-1:evtOff+1 (e.g. 4s)
%       3: No reward was used yet.

% load trial-by-trial behavior & task data tbytDat.mat
%parentDir = fileparts(filePath);

p = parse_postprocessing(filePath, varargin);

%  p = parse_postprocessing(filePath, {'imageToUse', 'dff_combined', 'channelOfInterest', 'green'});

if ~ismember(p.Results.imageToUse, {'ome_stack', 'combinedStack', 'dff_combined'})
    error('Image data type needs to be specified as dff_combined or ome_stack!')
end

%% whereabouts
match_header = regexp(filePath, 'm\d{1,4}_\d{6}', 'match');
header = match_header{1};

fileBeh = GrabFiles_sort_trials('tbytDat_parseGng', 0, {fullfile(filePath, 'Matfiles')});
if isempty(fileBeh{1})
    fileBeh = GrabFiles_sort_trials('tbytDat.mat', 1, {filePath});
end

% directory for imaging data
fileImg = GrabFiles_sort_trials('_img', 0, {fullfile(filePath)});

% directory for nidq behavioral data
% filePath_nidq =  GrabFiles_sort_trials('_g', 0, {filePath});

% file directory for trials
filePathTrials = fullfile(filePath, strcat(header, '_trials'));
if exist(filePathTrials, 'dir') == 0
    mkdir(filePathTrials);
end

% locate and load stimopts
fileStim = GrabFiles_sort_trials('stimInfo', 0, {fullfile(filePath)});
if isempty(fileStim)
    fileStim = uigetdir(filePath, 'Select the stimopts folder!');
    fileStim = GrabFiles_sort_trials('stimInfo', 0, {fileStim});
end
load(fileStim{1}, 'stimopts')
fprintf("finished loading stimopts!\n")

%% load data
% load evtIns for 'cmosExp' essential for temporal alignment of the frames
% load(fullfile(filePath_nidq{1}, 'evtInS'), 'evtInS')
% if exist('evtInS', 'var')~=1
%     [evtInS_file, evtInS_folder] = uigetfile('*.mat', 'Select the evtInS file for cmosExp!', filePath);
%     load(fullfile(evtInS_folder, evtInS_file), 'evtInS')
% end
% fprintf("finished loading evtInS!\n")

% load tbytDat that contains trial-by-trial timestamps (cmosExp, faceCam)
load(fullfile(fileBeh{1}), 'tbytDat')
fprintf("finished loading tbytDat!\n")

%% get trial type info
% trial type info
if ~isfield(tbytDat, 'pos_rwd_tr')
    tt = parse_gng_trialtype(stimopts);
    ttName = fieldnames(tt);
    for ff = 1:length(ttName)
        ttC = num2cell(tt.(ttName{ff}));
        [tbytDat.(ttName{ff})] = deal(ttC{:});
    end
end

%% get trial-by-trial dff
% compute and/or align dff, save each trial's dff
assert(strcmp(p.Results.imageToUse, 'dff_combined'))
% load preprocessed dffs ('dff_combined.mat')
[file_list_img, ~] = GrabFiles_sort_trials('green', 0, fileImg(1));
imgC = cell(1, length(file_list_img));
for ff = 1:length(file_list_img)
    [dff_file, ~] = GrabFiles_sort_trials('dff_combined.mat', 0, file_list_img(ff)); % use GrabFiles_sort_trials to sort both files and folders
    load(dff_file{1}, 'dff')
    imgC{1, ff} = dff;
    clearvars dff
    if exist('prepro_log', 'var')~=1
        load(dff_file{1}, 'opts')
    end
    fprintf("finished loading dff file #%d/%d\n", ff, length(file_list_img))
end

% dffCombinedProcessingDual 
tbytDat = dffCombinedProcessingDual(filePath, imgC, tbytDat, opts, p.Results);

dffCell = {tbytDat.dff};
dffsmCell = {tbytDat.dffsm};

assert(length(dffCell)==length(tbytDat))
assert(length(dffsmCell)==length(tbytDat))

%% do additional alignments
for tt = 1:length(tbytDat)
    % timestamp each frame relative to the stim onset
    tbytDat(tt).frameStimOnI = check_timestamp_overlap(tbytDat(tt).frameT, tbytDat(tt).evtOn);
    tbytDat(tt).frameStimOffI = check_timestamp_overlap(tbytDat(tt).frameT, tbytDat(tt).evtOff);
    tbytDat(tt).frameStimI = zeros(length(tbytDat(tt).frameT), 1);
    tbytDat(tt).frameStimI(find(tbytDat(tt).frameStimOnI, 1):find(tbytDat(tt).frameStimOffI, 1), 1) = 1;
    % map temporal events to cmosExp pulses
    tbytDat(tt).frameLickI = check_timestamp_overlap(tbytDat(tt).frameT, tbytDat(tt).Lick);
    tbytDat(tt).frameWaterI = check_timestamp_overlap(tbytDat(tt).frameT, tbytDat(tt).water);
    tbytDat(tt).frameAirpuffI = check_timestamp_overlap(tbytDat(tt).frameT, tbytDat(tt).airpuff);
end

% record frames save in separate folders
%record_dff_frames_auditory(tbytDat, filePathTrials, 1)

% make dff videos
%frameRate = 3;
%make_dff_videos_auditory(tbytDat, filePathTrials, frameRate, 12)
%make_dff_videos_auditory_selectFrames(tbytDat, filePathTrials, frameRate, 4, [0 1])

% save tbytDat without dffs
tbytDat = rmfield(tbytDat, {'dff', 'dffsm'});
save(fullfile(filePath, 'Matfiles', [header, '_', p.Results.channelOfInterest, '_tbytDat_dff']), 'tbytDat')
save(fullfile(filePath, 'Matfiles', [header, '_', p.Results.channelOfInterest, '_dff_Collect']), 'dffCell', '-v7.3')
save(fullfile(filePath, 'Matfiles', [header, '_', p.Results.channelOfInterest, '_dff_smCollect']), 'dffsmCell', '-v7.3')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function p = parse_postprocessing(filePath, vargs)
        % parse input, and extract name-value pairs
        default_imageToUse = 'ome_stack'; % default image sets to use, it must be 'ome_stack' or 'dff_combined'
        default_dffMethod = 'movingavg'; % default method for dff is moving average
        default_bvCorrectLogic = true; % default logic for hemodynamic correction
        default_tbytBaseWinF = 'all'; % default baseline window for trial-by-trial dff ('all' is to use all frames with negative timestamps relativev to the aligned event), otherwise it should be numeric, e.g., 1 for 1 s
        default_movingWinF = 30; % default window for moving average dff (30 s)
        default_saveTbytDff = true;
        default_combineStacksOfSameSequenceLogic = true;
        default_channelOfInterest = 'green';

        p = inputParser; % create parser object
        addRequired(p, 'filePath')

        addParameter(p, 'imageToUse', default_imageToUse)
        addParameter(p, 'dffMethod', default_dffMethod)
        addParameter(p, 'bvCorrectLogic', default_bvCorrectLogic)
        addParameter(p, 'tbytBaseWinF', default_tbytBaseWinF)
        addParameter(p, 'movingWinF', default_movingWinF)
        addParameter(p, 'saveTbytDff', default_saveTbytDff)
        addParameter(p, 'combineStacksOfSameSequenceLogic', default_combineStacksOfSameSequenceLogic)
        addParameter(p, 'channelOfInterest', default_channelOfInterest)

        parse(p, filePath, vargs{:})
    end


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


    function tt = parse_gng_trialtype(stimopts)
        tt.pos_rwd_tr = zeros(length(stimopts.rewarded_stim), 1);
        tt.pos_oms_tr = zeros(length(stimopts.rewarded_stim), 1);

        if any(find(stimopts.rewarded_stim==1))
            posTr_outcome = [find(stimopts.rewarded_stim==1), stimopts.outcome_positive_stim];
            posTr_outcomeI = posTr_outcome(:, 2) == 1;
            tt.pos_rwd_tr(posTr_outcome(posTr_outcomeI, 1)) = 1;
            tt.pos_oms_tr(posTr_outcome(~posTr_outcomeI, 1)) = 1;
        end

        tt.neg_pns_tr = zeros(length(stimopts.punished_stim), 1);
        tt.neg_oms_tr = zeros(length(stimopts.punished_stim), 1);

        if any(find(stimopts.punished_stim==1))
            negTr_outcome = [find(stimopts.punished_stim==1), stimopts.outcome_negative_stim];
            negTr_outcomeI = negTr_outcome(:, 2) == -1;
            tt.neg_pns_tr(negTr_outcome(negTr_outcomeI, 1)) = 1;
            tt.neg_oms_tr(negTr_outcome(~negTr_outcomeI, 1)) = 1;
        end
    end

end




