%% ========================================================================
function [Xpix, ok, statusMsg] = standardizeDffBlockToPixelsByFrames(X, nPixS, validPixI, blk)
% standardizeDffBlockToPixelsByFrames
%
% Converts a dffC block into pixels x frames format.
%
% Supported input formats:
%   1) height x width x frames, e.g. 64 x 64 x 860
%   2) pixels x frames, e.g. 4096 x 860
%   3) frames x pixels, e.g. 860 x 4096
%   4) validPixels x frames
%   5) frames x validPixels
%
% Output:
%   Xpix : nPixS x frames

Xpix = [];
ok = false;

if isempty(X)
    statusMsg = 'empty block';
    return;
end

X = double(X);

if ndims(X) == 3

    % Expected common case:
    %   64 x 64 x frames
    nFrames = size(X, 3);
    Xpix = reshape(X, [], nFrames);

    if size(Xpix, 1) ~= nPixS
        statusMsg = sprintf( ...
            '3D block reshaped to [%d x %d], but S expects %d pixels', ...
            size(Xpix, 1), size(Xpix, 2), nPixS);
        Xpix = [];
        return;
    end

    ok = true;
    statusMsg = sprintf('block %d: 3D image stack converted to pixels x frames', blk);
    return;
end

if ~ismatrix(X)
    statusMsg = sprintf('block %d: unsupported dimensionality %s', blk, mat2str(size(X)));
    return;
end

if size(X, 1) == nPixS

    % Already pixels x frames.
    Xpix = X;
    ok = true;
    statusMsg = sprintf('block %d: already pixels x frames', blk);
    return;

elseif size(X, 2) == nPixS

    % frames x pixels.
    Xpix = X';
    ok = true;
    statusMsg = sprintf('block %d: transposed frames x pixels to pixels x frames', blk);
    return;

elseif size(X, 1) == sum(validPixI)

    % validPixels x frames.
    Xpix = zeros(nPixS, size(X, 2));
    Xpix(validPixI, :) = X;
    ok = true;
    statusMsg = sprintf('block %d: validPixels x frames inserted into full pixel space', blk);
    return;

elseif size(X, 2) == sum(validPixI)

    % frames x validPixels.
    X = X';
    Xpix = zeros(nPixS, size(X, 2));
    Xpix(validPixI, :) = X;
    ok = true;
    statusMsg = sprintf('block %d: frames x validPixels transposed and inserted into full pixel space', blk);
    return;

else

    statusMsg = sprintf( ...
        'block %d: could not match dffC size %s to S pixel count %d or valid pixel count %d', ...
        blk, mat2str(size(X)), nPixS, sum(validPixI));
    return;
end

end