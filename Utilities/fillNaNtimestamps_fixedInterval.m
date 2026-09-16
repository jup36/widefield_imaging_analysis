function filledTimestamps = fillNaNtimestamps_fixedInterval(timestamps, expectedInterval)
% fillNaNtimestamps_fixedInterval
%
% Fill NaN values in a timestamp vector assuming a fixed frame interval.
%
% Behavior:
%   1) Leading NaN chunks are extrapolated backward from the first valid timestamp.
%   2) Interior NaN chunks are filled backward from the next valid timestamp.
%   3) Trailing NaN chunks are extrapolated forward from the last valid timestamp.
%
% This intentionally mirrors the original fillNaNtimestamps behavior, but
% keeps the expected interval explicit.
%
% Input:
%   timestamps       : numeric vector containing timestamps with NaNs
%   expectedInterval : frame interval in seconds, e.g. 0.1
%
% Output:
%   filledTimestamps : same size as input, with NaNs filled where possible
if nargin < 2 || isempty(expectedInterval)
    expectedInterval = 0.1;
end
% Preserve input shape.
originalSize = size(timestamps);
filledTimestamps = timestamps(:)';  % work as row vector
n = numel(filledTimestamps);
% If all NaN, cannot infer anything.
if all(isnan(filledTimestamps))
    warning('fillNaNtimestamps_fixedInterval:AllNaN', ...
        'Input timestamps are all NaN. Returning unchanged.');
    filledTimestamps = reshape(filledTimestamps, originalSize);
    return;
end
%% 1) Fill leading NaNs using first valid timestamp
firstValidIdx = find(~isnan(filledTimestamps), 1, 'first');
if firstValidIdx > 1
    numLeading = firstValidIdx - 1;
    firstValidTime = filledTimestamps(firstValidIdx);
    filledTimestamps(1:firstValidIdx-1) = ...
        firstValidTime - (numLeading:-1:1) * expectedInterval;
end
%% 2) Fill interior NaN chunks using next valid timestamp
idx = firstValidIdx + 1;
while idx <= n
    if isnan(filledTimestamps(idx))
        startNaN = idx;
        % Move until first non-NaN after this chunk.
        while idx <= n && isnan(filledTimestamps(idx))
            idx = idx + 1;
        end
        if idx <= n
            % Interior chunk: fill backward from next valid time.
            nextValidTime = filledTimestamps(idx);
            numMissing = idx - startNaN;
            filledTimestamps(startNaN:idx-1) = ...
                nextValidTime - (numMissing:-1:1) * expectedInterval;
        else
            % Trailing chunk: handled below.
            break;
        end
    else
        idx = idx + 1;
    end
end
%% 3) Fill trailing NaNs using last valid timestamp
lastValidIdx = find(~isnan(filledTimestamps), 1, 'last');
if lastValidIdx < n
    numTrailing = n - lastValidIdx;
    lastValidTime = filledTimestamps(lastValidIdx);
    filledTimestamps(lastValidIdx+1:n) = ...
        lastValidTime + (1:numTrailing) * expectedInterval;
end
% Restore original shape.
filledTimestamps = reshape(filledTimestamps, originalSize);
end