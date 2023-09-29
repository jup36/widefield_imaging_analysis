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