function [hfig, hBaseImg, hOverlay] = showImgAndMask(img, maskedImg, varargin)

if nargin < 3 || ~isa(varargin{1}, 'matlab.ui.Figure')
    % Prepare the figure for displaying the images
    hfig = figure;
else
    hfig = varargin{1};
end


% Scale the main image for better visibility
minImg = min(img(:));
maxImg = max(img(:)) / 3; % Enhance contrast by adjusting the max value
scaledImage = (img - minImg) / (maxImg - minImg);

hAxes = axes('Parent', hfig);

hBaseImg = imshow(scaledImage, 'Parent', hAxes); % Display the scaled base image
hold on;

% Assuming maskedImg contains areas from img multiplied by a binary mask
% Convert maskedImg to the same scale as scaledImage for consistent display
minMasked = min(maskedImg(:));
maxMasked = max(maskedImg(:)) / 3;
scaledMaskedImg = (maskedImg - minMasked) / (maxMasked - minMasked);
% Create an RGB overlay where 0 parts of maskedImg are green
%greenOverlay = zeros(size(scaledMaskedImg, 1), size(scaledMaskedImg, 2), 3); % Initialize with zeros
blueOverlay = zeros(size(scaledMaskedImg, 1), size(scaledMaskedImg, 2), 3); % Initialize with zeros
% Set the green channel to 1 (or another value for a different shade of green) where maskedImg is 0
%greenOverlay(:,:,2) = (maskedImg == 0);
blueOverlay(:,:,3) = (maskedImg == 0);
% Overlay maskedImg with transparency
% Note: Ensure scaledMaskedImg is not all zeros; otherwise, it won't display
hOverlay = imshow(blueOverlay, 'Parent', hAxes);
set(hOverlay, 'AlphaData', 0.25); % Adjust alpha for desired transparency

hold off;
title('Image with Transparent Overlay');
end
