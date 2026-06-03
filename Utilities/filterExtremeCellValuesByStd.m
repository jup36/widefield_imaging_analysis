function [valueC_filt, outlierInfo] = filterExtremeCellValuesByStd(valueC, nStdThresh, labelText)
% filterExtremeCellValuesByStd
%
% Removes extreme scalar values from a mouse × session cell array using a
% pooled mean/std threshold.
%
% INPUT
%   valueC
%       Cell array, usually nAnimals × nSessions.
%       Each cell contains a scalar numeric value, vector numeric value,
%       NaN, or empty.
%
%   nStdThresh
%       Threshold in standard deviations, e.g. 3.
%
%   labelText
%       Optional label used for printed diagnostics.
%
% OUTPUT
%   valueC_filt
%       Same size as valueC. Outlier cells are replaced with NaN.
%
%   outlierInfo
%       Struct containing threshold information and outlier indices.

if nargin < 2 || isempty(nStdThresh)
    nStdThresh = 3;
end

if nargin < 3 || isempty(labelText)
    labelText = 'cell values';
end

valueC_filt = valueC;

% Collect all scalar/session-level values
allVals = [];
cellValMat = nan(size(valueC));

for i = 1:numel(valueC)

    val = valueC{i};

    if isempty(val)
        continue
    end

    if isnumeric(val) && isscalar(val)
        thisVal = val;
    elseif isnumeric(val)
        thisVal = mean(val(:), 'omitnan');
    else
        continue
    end

    if ~isnan(thisVal)
        allVals(end+1, 1) = thisVal; %#ok<AGROW>
        cellValMat(i) = thisVal;
    end
end

mu = mean(allVals, 'omitnan');
sigma = std(allVals, 'omitnan');

lowerBound = mu - nStdThresh * sigma;
upperBound = mu + nStdThresh * sigma;

% If sigma is zero/invalid, do not remove anything
if isempty(allVals) || isnan(sigma) || sigma == 0
    outlierInfo = struct();
    outlierInfo.labelText = labelText;
    outlierInfo.nStdThresh = nStdThresh;
    outlierInfo.mean = mu;
    outlierInfo.std = sigma;
    outlierInfo.lowerBound = lowerBound;
    outlierInfo.upperBound = upperBound;
    outlierInfo.outlierLinearIdx = [];
    outlierInfo.outlierSubIdx = [];
    outlierInfo.outlierValues = [];
    outlierInfo.nOutliers = 0;

    fprintf('\n%s | Outlier filtering skipped: no valid variance.\n', labelText);
    return
end

outlierI = cellValMat < lowerBound | cellValMat > upperBound;

outlierLinearIdx = find(outlierI);
[outlierRows, outlierCols] = ind2sub(size(valueC), outlierLinearIdx);

% Replace outliers with NaN
for ii = 1:numel(outlierLinearIdx)
    valueC_filt{outlierLinearIdx(ii)} = NaN;
end

outlierInfo = struct();
outlierInfo.labelText = labelText;
outlierInfo.nStdThresh = nStdThresh;
outlierInfo.mean = mu;
outlierInfo.std = sigma;
outlierInfo.lowerBound = lowerBound;
outlierInfo.upperBound = upperBound;
outlierInfo.outlierLinearIdx = outlierLinearIdx;
outlierInfo.outlierSubIdx = [outlierRows(:), outlierCols(:)];
outlierInfo.outlierValues = cellValMat(outlierLinearIdx);
outlierInfo.nOutliers = numel(outlierLinearIdx);

fprintf('\n%s | Outlier filtering\n', labelText);
fprintf('mean = %.4f | std = %.4f | bounds = [%.4f, %.4f] | removed %d/%d values\n', ...
    mu, sigma, lowerBound, upperBound, outlierInfo.nOutliers, numel(allVals));

if outlierInfo.nOutliers > 0
    fprintf('Removed outliers at row/col indices:\n');
    disp(array2table([outlierRows(:), outlierCols(:), outlierInfo.outlierValues(:)], ...
        'VariableNames', {'AnimalRow', 'SessionCol', 'Value'}));
end

end