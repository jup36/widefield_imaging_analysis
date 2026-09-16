function [trace, info] = globalMeanDA(X, imageSize, mask)
% Full arrays use MATLAB column-major pixel order, as in the source pipeline.
info = struct('inputSize', size(X), 'nSelectedPixels', 0, ...
    'nFinitePixelsPerFrame', [], 'validPixelMask', mask, 'layout', 'empty');
trace = [];
if isempty(X), return; end
nPix = prod(imageSize);
if isempty(mask)
    keep = true(nPix,1);
else
    assert(islogical(mask) && numel(mask) == nPix, ...
        'Mask must be logical with %d elements; true means include.', nPix);
    keep = mask(:);
    assert(any(keep), 'Cortical mask selects no pixels.');
end
X = double(X);
if ndims(X) == 3 || isequal(size(X), imageSize)
    assert(size(X,1) == imageSize(1) && size(X,2) == imageSize(2), ...
        'Movie image shape differs from imageSize.');
    X = reshape(X, nPix, []);
    X = X(keep,:);
    info.layout = 'image stack';
elseif ismatrix(X)
    candidates = [size(X,1)==nPix, size(X,2)==nPix, ...
        ~isempty(mask) && size(X,1)==sum(keep), ...
        ~isempty(mask) && size(X,2)==sum(keep)];
    rowMatch = candidates(1) || candidates(3);
    colMatch = candidates(2) || candidates(4);
    assert(xor(rowMatch,colMatch), ...
        'Cannot uniquely determine pixel axis of %s; check imageSize/mask.', mat2str(size(X)));
    if colMatch, X = X.'; end
    if size(X,1) == nPix
        X = X(keep,:);
        info.layout = 'full pixels';
    else
        % Assumes compact rows correspond to find(mask(:)), as in source code.
        info.layout = 'compact selected pixels';
    end
else
    error('Unsupported movie shape: %s.', mat2str(size(X)));
end
X(~isfinite(X)) = NaN;
info.nSelectedPixels = size(X,1);
info.nFinitePixelsPerFrame = sum(isfinite(X),1);
trace = mean(X,1,'omitnan');
trace(info.nFinitePixelsPerFrame == 0) = NaN;
end