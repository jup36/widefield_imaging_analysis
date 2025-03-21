function record_dff_face_frames_auditory_gng_downsample_faceCam(filePath, varargin)
% filePath = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA019/DA019_042224';
% e.g., dff = tbytDat(21).dffsm
%By default this function iterate through all trials and record frames of
% each trial as png files. It can also run specific trials designated as varargin.
% Note that face videos

%% load tbytDat
[~, header] = fileparts(filePath);
filePathTrials = fullfile(filePath, strcat(header, '_trials'));
fileBeh = GrabFiles_sort_trials('tbytDat_dff', 0, {fullfile(filePath, 'Matfiles')});
if isempty(fileBeh{1})
    fileBeh = GrabFiles_sort_trials('tbytDat.mat', 1, {filePath});
end
load(fullfile(fileBeh{1}), 'tbytDat')

if ~isfield(tbytDat, 'rewardTrI')
    tbytDat = addRewardTrI(tbytDat);
end

if ~isfield(tbytDat, 'punishTrI')
    tbytDat = addPunishTrI(tbytDat);
end

% load evtInS
fileNidq = GrabFiles_sort_trials('_g', 0, {filePath});
if ~isempty(fileNidq) && isscalar(fileNidq)
    fileEvt = GrabFiles_sort_trials('evtInS', 0, fileNidq(1));
    load(fullfile(fileEvt{1}), 'evtInS');
else
    [evtInS_file, evtInS_folder] = uigetfile('*.mat', 'Select the evtInS file for faceCam!', filePath);
    load(fullfile(evtInS_folder, evtInS_file), 'evtInS')
end

%% list face videos
faceVidPath = fullfile(filePath, strcat(header, '_vid_facecropped'));
faceVids = GrabFiles_sort_trials('behvid', 0, {faceVidPath});

%% sort trials to work with
if isempty(varargin)
    trials = 1:length(tbytDat);
else
    trials = [varargin{1}(:)]';
end

%% record frames first
downSampleFactor = 4; % default downsample factor to downsample the faceCam frames
close all;
for t = trials
    pngDir = GrabFiles_sort_trials(sprintf('trial%d_', t), 0, {filePathTrials});
    if isempty(pngDir)
        pngDir = fullfile(filePathTrials, sprintf('trial%d', t));
        mkdir(fullfile(filePathTrials, sprintf('trial%d', t)));
    else
        pngDir = pngDir{1};
    end

    % load dff
    load(fullfile(pngDir, 'tbytDff.mat'), 'tbytDffsm');

    % interpolate dff (15Hz) to match the faceCam sampling rate (200Hz)
    faceCamI = tbytDat(t).cmosExpTrainI;
    faceCamTs = evtInS.faceCam(evtInS.faceCam(:, 2)==faceCamI, :);

    targetFaceTsFstI = find(faceCamTs(:, 1) < tbytDat(t).frameT(1), 1, 'last');
    targetFaceTsLastI = find(faceCamTs(:, 1) > tbytDat(t).frameT(end), 1, 'first');
    faceFramesToRead = targetFaceTsFstI:downSampleFactor:targetFaceTsLastI; % downsampled

   v = VideoReader(faceVids{faceCamI});
   totalFaceFrames = v.Duration*v.FrameRate;

   if ~isempty(targetFaceTsFstI) && ~isempty(targetFaceTsLastI) && targetFaceTsLastI <= totalFaceFrames
        targetFaceTs = faceCamTs(faceFramesToRead, 1); % downsampled
        tbytDat(t).resampledVidFrameTs = targetFaceTs; 

       tbytDffsmIntp = temporalIntpImageFrames(tbytDffsm, tbytDat(t).frameT, targetFaceTs);

        % Read, concatenate, and write each frame
        for fr = 1:length(faceFramesToRead)
            if hasFrame(v)
                pngName = sprintf('faceDff_frame_%d.png', faceFramesToRead(fr));

                % Set the CurrentTime property to the time corresponding to the frame
                v.CurrentTime = (faceFramesToRead(fr) - 1) / v.FrameRate;

                % Read the face frame
                faceFrame = readFrame(v);

                % Resize (interpolate) dff to match faceFrame
                tbytDffsmIntpRs = imresize(tbytDffsmIntp(:, :, fr), [size(faceFrame, 1), size(faceFrame, 2)]);

                % Replace NaNs with a value that corresponds to black
                tbytDffsmIntpRs(isnan(tbytDffsmIntpRs)) = -Inf;

                % Create a colormap
                cmap = [0, 0, 0; hot(256)];
                dffImage = ind2rgb(im2uint8(mat2gray(tbytDffsmIntpRs, [-2 2])), cmap);

                % Insert text to indicate tone presentation
                cueOnT = tbytDat(t).frameT(find(tbytDat(t).frameStimI, 1, 'first'));
                cueOffT = tbytDat(t).frameT(find(tbytDat(t).frameStimI, 1, 'last'));
                position = [size(dffImage, 2) / 2, 5]; % [x, y] position

                % Check for tone presentation
                if targetFaceTs(fr) >= cueOnT && targetFaceTs(fr) <= cueOffT
                    if tbytDat(t).rewardTrI == 1 % Reward tone
                        dffImage = insertText(dffImage, position, "Tone1", ...
                            'AnchorPoint', 'CenterTop', 'FontSize', 20, ...
                            'TextColor', 'white', 'BoxColor', [0 0 255], ...
                            'BoxOpacity', 0.4);
                    elseif tbytDat(t).punishTrI == 1 % Punishment tone
                        dffImage = insertText(dffImage, position, "Tone2", ...
                            'AnchorPoint', 'CenterTop', 'FontSize', 20, ...
                            'TextColor', 'white', 'BoxColor', [255 0 0], ...
                            'BoxOpacity', 0.4);
                    end
                end

                % Concatenate dff and face frames horizontally
                dffFaceConcat = cat(2, im2uint8(dffImage), faceFrame);

                % Save the image as png
                imwrite(dffFaceConcat, fullfile(pngDir, pngName));

                fprintf('Frame #%d/%d of trial #%d is labeled and saved.\n', fr, length(faceFramesToRead), t);
            end
        end
    end
    clear v;
    fprintf('Trial #%d is completed.\n', t);
end

save(fullfile(fileBeh{1}), 'tbytDat') % to save tbytDat(t).faceCamTrelFrames

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function middleArrayI = middleIndex(originalArrayLength, numMiddlePoints)
        %This function gets the index to select the middle portion of
        %   the originalArrayLength that corresponds to
        %   numMiddlePoints.

        if originalArrayLength >= numMiddlePoints
            % originalArrayLength = length(targetFaceTs);
            % Calculate start and end indices
            startIndex = ceil((originalArrayLength - numMiddlePoints) / 2) + 1;
            endIndex = startIndex + numMiddlePoints - 1;

            middleArrayI = zeros(originalArrayLength, 1);

            % Select the middle points
            middleArrayI(startIndex:endIndex)=1;
        else
            error("The specified middle points exceeds the original array length!")
        end
    end


    function newMatrix = temporalIntpImageFrames( imgFrameMat, originalTimePoints, targetTimePoints )
        % Assuming tbytDffsm is your original 3D matrix
        % And tbytDat(t).frameTrel and tbytDat(t).fraceCamTrel are your time vectors

        %originalTimePoints = tbytDat(t).frameTrel(valDffFrameMinI:valDffFrameMaxI);
        %targetTimePoints = tbytDat(t).faceCamTrel;

        % Get the size of the original matrix
        [height, width, ~] = size(imgFrameMat);

        % Preallocate the new 3D matrix
        newMatrix = nan(height, width, length(targetTimePoints));

        % Interpolate for each pixel
        for i = 1:height
            for j = 1:width
                pixelTimeSeries = squeeze(imgFrameMat(i, j, :));
                if sum(~isnan(pixelTimeSeries))==length(pixelTimeSeries)
                    newMatrix(i, j, :) = interp1(originalTimePoints, pixelTimeSeries, targetTimePoints, 'linear', 'extrap');
                else
                    newMatrix(i, j, :) = nan(length(targetTimePoints), 1);
                end
            end
        end


    end

end



