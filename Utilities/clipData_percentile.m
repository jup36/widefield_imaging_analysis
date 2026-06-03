function data_clean = clipData_percentile(data, lowP, highP)
% clipData_percentile
%
% Clip a numeric data matrix/array based on global percentile thresholds.
%
% Usage:
%   data_clean = clipData_percentile(data)
%   data_clean = clipData_percentile(data, 0.1, 99.9)
%
% Inputs:
%   data  : numeric matrix/array
%   lowP  : lower percentile, default = 0.1
%   highP : upper percentile, default = 99.9
%
% Output:
%   data_clean : data clipped to [lowP, highP] percentile range

if nargin < 2 || isempty(lowP)
    lowP = 0.1;
end

if nargin < 3 || isempty(highP)
    highP = 99.9;
end

if ~isnumeric(data)
    error('Input data must be numeric.');
end

if lowP < 0 || highP > 100 || lowP >= highP
    error('Percentiles must satisfy 0 <= lowP < highP <= 100.');
end

% Collect valid values
allVals = data(:);
allVals = allVals(isfinite(allVals));

if isempty(allVals)
    warning('No finite values found. Returning input unchanged.');
    data_clean = data;
    return
end

% Compute clipping bounds
lo = prctile(allVals, lowP);
hi = prctile(allVals, highP);

fprintf('Clipping range: [%.3f, %.3f]\n', lo, hi);

% Apply clipping while preserving NaNs/Infs as-is unless finite and out of range
data_clean = data;

finiteMask = isfinite(data_clean);

data_clean(finiteMask & data_clean < lo) = lo;
data_clean(finiteMask & data_clean > hi) = hi;

end