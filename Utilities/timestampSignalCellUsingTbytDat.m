function [sigCell, timingInfo] = timestampSignalCellUsingTbytDat(sigCell, tbytDat, ledPrefix, expectedInterval)
% timestampSignalCellUsingTbytDat
%
% Fills row 2 of a 2-row block signal cell with frame timestamps.
%
% Input:
%   sigCell{1, blk}
%       K x frames signal matrix.
%
%   sigCell{2, blk}
%       Empty or NaN 1 x frames timestamp vector.
%
% Output:
%   sigCell{2, blk}
%       1 x frames timestamp vector.
%
% Strategy:
%   1) Assign LED-derived timestamps directly to frame indices.
%   2) If pulseTime length does not match frame range length, do NOT stretch
%      or interpolate pulseTime. Assign only matched samples conservatively.
%   3) Fill remaining NaNs using a fixed expected interval, default 0.1 sec.
%
% This mirrors the earlier hC-style timestamping logic and avoids the
% final-frame stretching artifact caused by interpolating pulseTime.
if nargin < 3 || isempty(ledPrefix)
    ledPrefix = 'blue';
end
if nargin < 4 || isempty(expectedInterval)
    expectedInterval = 0.1;
end
[trainField, pulseEdgeField, pulseTimeField, resolvedPrefix] = ...
    resolveLedFields(tbytDat, ledPrefix);
blockIds = getNumericFieldColumn(tbytDat, trainField);
nBlocks = size(sigCell, 2);
% Initialize timingInfo.
timingInfo = struct();
timingInfo.ledPrefixRequested = ledPrefix;
timingInfo.ledPrefixResolved = resolvedPrefix;
timingInfo.trainField = trainField;
timingInfo.pulseEdgeField = pulseEdgeField;
timingInfo.pulseTimeField = pulseTimeField;
timingInfo.expectedInterval = expectedInterval;
timingInfo.timestampFillMethod = 'direct_LED_assignment_then_fixed_interval_fill';
timingInfo.block = repmat(struct( ...
    'nFrames', NaN, ...
    'nTrialsAssigned', 0, ...
    'nFramesAssigned', 0, ...
    'nMismatchedTrials', 0, ...
    'dtMedianAfterFill', NaN, ...
    'nBadDtAfterFill', NaN), ...
    1, nBlocks);
% Preallocate timestamp rows for every non-empty block.
for blk = 1:nBlocks
    if ~isempty(sigCell{1, blk})
        nFrames = size(sigCell{1, blk}, 2);
        sigCell{2, blk} = NaN(1, nFrames);
        timingInfo.block(blk).nFrames = nFrames;
    else
        sigCell{2, blk} = [];
        timingInfo.block(blk).nFrames = 0;
    end
end
% Assign LED-derived timestamps directly.
for tr = 1:numel(tbytDat)
    blk = blockIds(tr);
    if ~isfinite(blk) || blk ~= fix(blk) || blk < 1 || blk > nBlocks
        continue;
    end
    if isempty(sigCell{1, blk})
        continue;
    end
    nFrames = size(sigCell{1, blk}, 2);
    pulseEdge = getNumericVector(tbytDat(tr).(pulseEdgeField));
    pulseTime = getNumericVector(tbytDat(tr).(pulseTimeField));
    if isempty(pulseEdge) || isempty(pulseTime)
        continue;
    end
    pulseEdge = round(pulseEdge(:)');
    pulseTime = pulseTime(:)';
    frameIdxAll = pulseEdge(1):pulseEdge(end);
    if isempty(frameIdxAll)
        continue;
    end
    % Do not interpolate or stretch.
    % Just assign the matched portion.
    nAssign = min(numel(frameIdxAll), numel(pulseTime));
    if numel(frameIdxAll) ~= numel(pulseTime)
        timingInfo.block(blk).nMismatchedTrials = ...
            timingInfo.block(blk).nMismatchedTrials + 1;
    end
    frameIdxUse = frameIdxAll(1:nAssign);
    pulseTimeUse = pulseTime(1:nAssign);
    % Keep only valid frame indices.
    validFrameI = frameIdxUse >= 1 & frameIdxUse <= nFrames;
    frameIdxUse = frameIdxUse(validFrameI);
    pulseTimeUse = pulseTimeUse(validFrameI);
    if isempty(frameIdxUse)
        continue;
    end
    sigCell{2, blk}(frameIdxUse) = pulseTimeUse;
    timingInfo.block(blk).nTrialsAssigned = ...
        timingInfo.block(blk).nTrialsAssigned + 1;
    timingInfo.block(blk).nFramesAssigned = ...
        timingInfo.block(blk).nFramesAssigned + numel(frameIdxUse);
end
% Fill remaining NaNs using fixed frame interval.
for blk = 1:nBlocks
    if isempty(sigCell{1, blk})
        continue;
    end
    nFrames = size(sigCell{1, blk}, 2);
    sigCell{2, blk} = sigCell{2, blk}(1, 1:nFrames);
    sigCell{2, blk} = fillNaNtimestamps_fixedInterval( ...
        sigCell{2, blk}, expectedInterval);
    % Lightweight diagnostic.
    dt = diff(sigCell{2, blk});
    assert(~any(isfinite(dt) & dt <= 0), ...
        'Non-increasing DA timestamps in block %d; inspect LED metadata before alignment.', blk);
    dtFinite = dt(isfinite(dt));
    if ~isempty(dtFinite)
        timingInfo.block(blk).dtMedianAfterFill = median(dtFinite, 'omitnan');
        badDtI = isfinite(dt) & ...
            (dt <= 0 | abs(dt - expectedInterval) > 0.02);
        timingInfo.block(blk).nBadDtAfterFill = sum(badDtI);
    end
    fprintf(['Timestamped block #%d/%d using %s LED fields; ' ...
        'fixed %.3f sec NaN fill; mismatched trials = %d; bad dt steps = %d.\n'], ...
        blk, nBlocks, resolvedPrefix, expectedInterval, ...
        timingInfo.block(blk).nMismatchedTrials, ...
        timingInfo.block(blk).nBadDtAfterFill);
end
end
