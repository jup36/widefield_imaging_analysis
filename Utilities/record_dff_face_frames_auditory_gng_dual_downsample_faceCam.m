function record_dff_face_frames_auditory_gng_dual_downsample_faceCam(filePath, varargin)
% filePath = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA019/DA019_042224';
% e.g., dff = tbytDat(21).dffsm
%By default this function iterate through all trials and record frames of
% each trial as png files. It can also run specific trials designated as varargin.
% Note that face videos

%% load tbytDat
match_header = regexp(filePath, 'm\d{1,4}_\d{6}', 'match');
header = match_header{1};

filePathTrials = fullfile(filePath, strcat(header, '_trials'));

fileBehG= GrabFiles_sort_trials('green_tbytDat_dff', 0, {fullfile(filePath, 'Matfiles')});
if isempty(fileBehG{1})
    fileBehG = GrabFiles_sort_trials('green_tbytDat_dff', 1, {filePath});
end
tbytDatG = load(fullfile(fileBehG{1}), 'tbytDat');
tbytDatG = tbytDatG.('tbytDat'); clearvars tbytDat 

fileBehR= GrabFiles_sort_trials('red_tbytDat_dff', 0, {fullfile(filePath, 'Matfiles')});
if isempty(fileBehR{1})
    fileBehR = GrabFiles_sort_trials('red_tbytDat_dff', 1, {filePath});
end
tbytDatR = load(fullfile(fileBehR{1}), 'tbytDat');
tbytDatR = tbytDatR.('tbytDat'); clearvars tbytDat 

tbytDatG = addRewardTrI(tbytDatG);
tbytDatR = addRewardTrI(tbytDatR);
tbytDatG = addPunishTrI(tbytDatG);
tbytDatR = addPunishTrI(tbytDatR);

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
faceVidPath = GrabFiles_sort_trials([header, '*', '_vid'], 0, {filePath});
faceVids = GrabFiles_sort_trials('behvid', 0, faceVidPath);

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
    pngSubDirG = GrabFiles_sort_trials(['dffG*', sprintf('trial_%d', t)], 0, {filePathTrials});
    pngSubDirR = GrabFiles_sort_trials(['dffR*', sprintf('trial_%d', t)], 0, {filePathTrials});

    % load dff green
    tbytDffsmG = load(fullfile(pngSubDirG{1}, 'tbytDff.mat'), 'tbytDffsm');
    tbytDffsmG = tbytDffsmG.('tbytDffsm');

    tbytDffsmR = load(fullfile(pngSubDirR{1}, 'tbytDff.mat'), 'tbytDffsm');
    tbytDffsmR = tbytDffsmR.('tbytDffsm');

    % interpolate dff (10Hz) to match the faceCam sampling rate (200Hz)
    assert(tbytDatG(t).sideGreenExpTrainI==tbytDatR(t).topRedExpTrainI)
    faceCamI = tbytDatG(t).sideGreenExpTrainI;
    faceCamTs = evtInS.faceCam(evtInS.faceCam(:, 2)==faceCamI, :);

    targetFaceTsFstI = find(faceCamTs(:, 1) < min(tbytDatG(t).frameT(1), tbytDatR(t).frameT(1)), 1, 'last');
    targetFaceTsLastI = find(faceCamTs(:, 1) > max(tbytDatG(t).frameT(end), tbytDatR(t).frameT(end)), 1, 'first');
    faceFramesToRead = targetFaceTsFstI:downSampleFactor:targetFaceTsLastI; % downsampled

    v = VideoReader(faceVids{faceCamI});
    totalFaceFrames = v.Duration*v.FrameRate;

    if ~isempty(targetFaceTsFstI) && ~isempty(targetFaceTsLastI) && targetFaceTsLastI <= totalFaceFrames
        targetFaceTs = faceCamTs(faceFramesToRead, 1); % downsampled
        tbytDatG(t).resampledVidFrameTs = targetFaceTs;
        tbytDatR(t).resampledVidFrameTs = targetFaceTs;

        tbytDffsmGIntp = temporalIntpImageFrames(tbytDffsmG, tbytDatG(t).frameT, targetFaceTs);
        tbytDffsmRIntp = temporalIntpImageFrames(tbytDffsmR, tbytDatR(t).frameT, targetFaceTs);

        % Read, concatenate, and write each frame
        for fr = 1:length(faceFramesToRead)
            if hasFrame(v)
                pngName = sprintf('faceDffGR_frame_%d.png', faceFramesToRead(fr));

                % Set the CurrentTime property to the time corresponding to the frame
                v.CurrentTime = (faceFramesToRead(fr) - 1) / v.FrameRate;

                % Read the face frame
                faceFrame = readFrame(v);

                % Resize (interpolate) dff to match faceFrame
                tbytDffGsmIntpRs = imresize(tbytDffsmGIntp(:, :, fr), [size(faceFrame, 1), size(faceFrame, 2)]);
                tbytDffRsmIntpRs = imresize(tbytDffsmRIntp(:, :, fr), [size(faceFrame, 1), size(faceFrame, 2)]);

                % Replace NaNs with a value that corresponds to black
                tbytDffGsmIntpRs(isnan(tbytDffGsmIntpRs)) = -Inf;
                tbytDffRsmIntpRs(isnan(tbytDffRsmIntpRs)) = -Inf;

                % Create a colormap
                cmap = [0, 0, 0; hot(256)];
                dffImageG = ind2rgb(im2uint8(mat2gray(tbytDffGsmIntpRs, [-2 2])), cmap);
                dffImageR = ind2rgb(im2uint8(mat2gray(tbytDffRsmIntpRs, [-5 5])), cmap);

                % Insert text to indicate tone presentation
                cueOnT = tbytDatG(t).frameT(find(tbytDatG(t).frameStimI, 1, 'first'));
                cueOffT = tbytDatG(t).frameT(find(tbytDatG(t).frameStimI, 1, 'last'));
                position = [size(dffImageG, 2) / 2, 5]; % [x, y] position

                % Check for tone presentation
                if targetFaceTs(fr) >= cueOnT && targetFaceTs(fr) <= cueOffT
                    if tbytDatG(t).rewardTrI == 1 % Reward tone
                        dffImageG = insertText(dffImageG, position, "Tone1", ...
                            'AnchorPoint', 'CenterTop', 'FontSize', 20, ...
                            'TextColor', 'white', 'BoxColor', [0 0 255], ...
                            'BoxOpacity', 0.4);
                    elseif tbytDatG(t).punishTrI == 1 % Punishment tone
                        dffImageG = insertText(dffImageG, position, "Tone2", ...
                            'AnchorPoint', 'CenterTop', 'FontSize', 20, ...
                            'TextColor', 'white', 'BoxColor', [255 0 0], ...
                            'BoxOpacity', 0.4);
                    end
                end

                % Concatenate dff and face frames horizontally
                dffGRFaceConcat = cat(2, im2uint8(dffImageG), im2uint8(dffImageR), faceFrame);

                % Save the image as png
                pngGRdir = fullfile(filePathTrials, ['faceDffGR_', sprintf('trial_%d', t)]); 
                if exist(pngGRdir)~=7
                    mkdir(pngGRdir)
                end
                
                imwrite(dffGRFaceConcat, fullfile(pngGRdir, pngName));

                fprintf('Frame #%d/%d of trial #%d is labeled and saved.\n', fr, length(faceFramesToRead), t);
            end
        end
    end
    clear v;
    fprintf('Trial #%d/%d is completed.\n', t, length(tbytDat));
end

save(fullfile(fileBehG{1}), 'tbytDatG', '-append') % to save tbytDat(t).faceCamTrelFrames
save(fullfile(fileBehR{1}), 'tbytDatR', '-append') % to save tbytDat(t).faceCamTrelFrames

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



