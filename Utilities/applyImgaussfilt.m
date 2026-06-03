function smFrames = applyImgaussfilt(frames, varargin)
% applyImgaussfilt applies frame-by-frame Gaussian smoothing while preserving NaNs.
%
% Usage:
%   smFrames = applyImgaussfilt(frames)
%   smFrames = applyImgaussfilt(frames, 'sigma', 2)
%
% Input:
%   frames : X x Y x T image stack
%
% Name-value:
%   'sigma' : Gaussian filter sigma, default = 1
%
% Output:
%   smFrames : smoothed image stack, same size as frames

% --------------------
% Parse inputs
% --------------------
p = inputParser;
p.addParameter('sigma', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.parse(varargin{:});

sigma = p.Results.sigma;

% --------------------
% Initialize output
% --------------------
smFrames = zeros(size(frames), 'like', frames);

% --------------------
% Smooth each frame
% --------------------
for frame = 1:size(frames, 3)

    % Extract current frame
    current_frame = frames(:, :, frame);

    % Create NaN mask
    nan_mask = isnan(current_frame);

    % Temporarily replace NaNs with mean of non-NaN values
    frame_mean = mean(current_frame(~nan_mask), 'omitnan');

    % Safety fallback in case the entire frame is NaN
    if isnan(frame_mean)
        smFrames(:, :, frame) = current_frame;
        continue;
    end

    current_frame(nan_mask) = frame_mean;

    % Apply Gaussian smoothing
    smFrame = imgaussfilt(current_frame, sigma);

    % Restore NaNs
    smFrame(nan_mask) = NaN;

    % Store smoothed frame
    smFrames(:, :, frame) = smFrame;
end

end