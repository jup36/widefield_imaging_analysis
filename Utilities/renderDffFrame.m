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
