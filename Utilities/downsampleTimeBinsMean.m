function [X_ds, timeX_ds] = downsampleTimeBinsMean(X, timeX, newBinSize)
%DOWNSAMPLETIMEBINSMEAN Downsample animal-by-time matrix by averaging time bins.
%
% [X_ds, timeX_ds] = downsampleTimeBinsMean(X, timeX, newBinSize)
%
% INPUT
%   X          : [nAnimals x nTime] numeric matrix
%                rows = animals, columns = time bins
%   timeX      : [1 x nTime] or [nTime x 1] vector of original time-bin centers
%   newBinSize : desired new bin size in same units as timeX
%                Example: 0.1
%
% OUTPUT
%   X_ds       : [nAnimals x nTime_ds] downsampled matrix
%   timeX_ds   : [1 x nTime_ds] new time-bin centers
%
% NOTES
%   - This function assumes uniformly spaced original time bins.
%   - Downsampling is done by averaging consecutive bins.
%   - If the number of bins is not divisible by the downsampling factor,
%     trailing bins are discarded.
%
% EXAMPLE
%   [fasts_ds, timeX_ds] = downsampleTimeBinsMean(muNoGo_absSlope_projNoGoAxis_fasts, timeX, 0.1);
%   [slows_ds, ~]        = downsampleTimeBinsMean(muNoGo_absSlope_projNoGoAxis_slows, timeX, 0.1);
%

%% Validate inputs
if ~isnumeric(X) || ndims(X) ~= 2
    error('X must be a 2D numeric matrix of size [nAnimals x nTime].');
end

if ~isnumeric(timeX) || ~isvector(timeX)
    error('timeX must be a numeric vector of time-bin centers.');
end

timeX = timeX(:)';  % force row vector

if size(X,2) ~= numel(timeX)
    error('size(X,2) must match numel(timeX).');
end

if ~isscalar(newBinSize) || ~isnumeric(newBinSize) || newBinSize <= 0
    error('newBinSize must be a positive numeric scalar.');
end

%% Check original bin size
dt = diff(timeX);
dt0 = median(dt);

tol = 1e-8;
if any(abs(dt - dt0) > tol)
    error('timeX must be uniformly spaced.');
end

%% Compute downsampling factor
factor = newBinSize / dt0;

if abs(factor - round(factor)) > 1e-8
    error('newBinSize must be an integer multiple of the original bin size. Original dt = %.6f, requested newBinSize = %.6f.', dt0, newBinSize);
end

factor = round(factor);

if factor < 1
    error('newBinSize must be >= original bin size.');
elseif factor == 1
    X_ds = X;
    timeX_ds = timeX;
    return;
end

%% Trim trailing bins if necessary
nTime = size(X,2);
nBlocks = floor(nTime / factor);
nKeep = nBlocks * factor;

if nKeep < nTime
    warning('Discarding %d trailing time bins so the data can be grouped into %d-bin blocks.', ...
        nTime - nKeep, factor);
end

X_trim = X(:, 1:nKeep);
timeX_trim = timeX(1:nKeep);

%% Reshape and average across grouped bins
% X_trim: [nAnimals x (nBlocks*factor)]
% reshape to [nAnimals x factor x nBlocks]
X_reshaped = reshape(X_trim, size(X,1), factor, nBlocks);

% average across the 2nd dimension (within each larger bin)
X_ds = squeeze(mean(X_reshaped, 2, 'omitnan'));

% If only one block remains, squeeze may collapse unexpectedly
if nBlocks == 1
    X_ds = reshape(X_ds, size(X,1), 1);
end

%% New time bin centers = mean of original centers within each block
timeX_ds = mean(reshape(timeX_trim, factor, nBlocks), 1, 'omitnan');

end