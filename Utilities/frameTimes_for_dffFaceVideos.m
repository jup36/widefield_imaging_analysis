function frameTimes_for_dffFaceVideos(filePath)
% filePath = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA019/DA019_042224';
% e.g., dff = tbytDat(21).dffsm
%By default this function iterate through all trials and record frames of
% each trial as png files. It can also run specific trials designated as varargin.
% Note that face videos

%% load tbytDat
match_header = regexp(filePath, 'm\d{1,4}_\d{6}', 'match');
header = match_header{1};

% filePathTrials = GrabFiles_sort_trials('_trials', 0, {filePath});

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

%% record frames first
downSampleFactor = 4; % default downsample factor to downsample the faceCam frames
close all;
for t = 1:length(tbytDatG)
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
    end
    clear v;
    fprintf('Trial #%d is completed.\n', t);
end

save(fullfile(fileBehG{1}), 'tbytDatG', '-append') % to save tbytDat(t).faceCamTrelFrames
save(fullfile(fileBehR{1}), 'tbytDatR', '-append') % to save tbytDat(t).faceCamTrelFrames

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function renderedImg = renderDffFrame(dffFrame, hemiMaskC, targetSize, colorAxis)
    %This script takes df/f image matrix (dffFrame) and hemi brain masks
    % (hemiMaskC) to render a masked df/f image, which is also resized to match
    % the target size. 
    % 5/29/2025 Junchol Park

        targetHeight = targetSize(1);
        targetWidth = targetSize(2);

        % Set up the figure
        fig = figure('Visible', 'off', ...
            'Units', 'pixels', ...
            'Position', [100, 100, targetWidth, targetHeight], ...
            'Color', 'k');  % figure background

        ax = axes(fig);
        set(ax, 'Units', 'normalized', ...
            'Position', [0 0 1 1], ...
            'Color', 'k');  % axes background
        axis off; hold on;

        % Upsample the image and masks to target size
        upsampledImg = imresize(dffFrame, targetSize, 'bilinear');
        upsampledMask = cellfun(@(x) imresize(double(x), targetSize, 'nearest') > 0.5, hemiMaskC, 'UniformOutput', false);

        % Optionally mask the image data (recommended for aesthetics)
        combinedMask = upsampledMask{1} | upsampledMask{2};
        maskedImg = upsampledImg;
        maskedImg(~combinedMask) = NaN;

        % Display the masked image
        imshow(maskedImg, colorAxis, 'InitialMagnification', 'fit');
        colormap(ax, magma);
        caxis(colorAxis);

        % Overlay just the borders — do NOT fill with black
        for hemi = 1:2
            cc = bwconncomp(upsampledMask{hemi}, 8);
            s = regionprops(cc, 'ConvexHull');
            if ~isempty(s)
                % Only plot white outline — no fill!
                plot(s(1).ConvexHull(:,1), s(1).ConvexHull(:,2), 'w', 'LineWidth', 3);
            end
        end

        % Set axis range and orientation
        axis([1 targetWidth 1 targetHeight]);
        set(gca, 'YDir', 'reverse');
        set(gca, 'Position', [0 0 1 1], 'Units', 'normalized');

        % Capture final image
        frameImg = getframe(gca);
        renderedImg = frameImg.cdata;

        close(fig);
    end


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