function h = overlayMaskOnFrame(frameImg, mask, varargin)
% overlayMaskOnFrame  Display a 2D image with a transparent mask overlay.
%
% h = overlayMaskOnFrame(frameImg, mask, 'Name', value, ...)
%
% INPUT
%   frameImg : 2D image
%   mask     : logical 2D mask, same size as frameImg
%
% NAME-VALUE PAIRS
%   'MaskColor'      : RGB triple for overlay color (default [1 1 0], yellow)
%   'MaskAlpha'      : transparency of mask where mask==true (default 0.35)
%   'DisplayRange'   : two-element vector for imagesc scaling (default [])
%   'TitleStr'       : figure title (default '')
%   'NewFigure'      : true/false, create new figure (default true)
%
% OUTPUT
%   h : struct with graphics handles

    p = inputParser;
    addParameter(p, 'MaskColor', [1 1 0], @(x) isnumeric(x) && numel(x)==3);
    addParameter(p, 'MaskAlpha', 0.35, @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
    addParameter(p, 'DisplayRange', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
    addParameter(p, 'TitleStr', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'NewFigure', true, @(x) islogical(x) || isnumeric(x));
    parse(p, varargin{:});

    if ~isequal(size(frameImg), size(mask))
        error('frameImg and mask must have the same size.');
    end

    if p.Results.NewFigure
        figure;
    end

    ax = gca;

    % base image
    if isempty(p.Results.DisplayRange)
        hImg = imagesc(frameImg, 'Parent', ax);
    else
        hImg = imagesc(frameImg, p.Results.DisplayRange, 'Parent', ax);
    end
    axis image;
    set(ax, 'YDir', 'normal');
    colormap(ax, turbo);
    colorbar;
    hold(ax, 'on');

    % create solid RGB overlay
    overlayRGB = zeros([size(mask), 3]);
    for c = 1:3
        overlayRGB(:,:,c) = p.Results.MaskColor(c);
    end

    % draw overlay with transparency only on masked pixels
    hMask = image(overlayRGB, 'Parent', ax);
    set(hMask, 'AlphaData', double(mask) * p.Results.MaskAlpha);

    if strlength(string(p.Results.TitleStr)) > 0
        title(p.Results.TitleStr);
    end

    h = struct();
    h.ax = ax;
    h.img = hImg;
    h.mask = hMask;
    set(gca, 'YDir', 'reverse')
end

function overlayMaskOutline(img, mask)

    imagesc(img);
    axis image;
    set(gca,'YDir','normal');
    colormap turbo;
    colorbar;
    hold on;

    B = bwboundaries(mask);
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5);
    end

end