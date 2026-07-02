function [sigCell, timingInfo] = timestampSignalCellUsingTbytDat(sigCell, tbytDat, ledPrefix)
% timestampSignalCellUsingTbytDat
%
% Fills row 2 of a 2-row block signal cell with frame timestamps.
%
% Strategy:
%   1) Assign LED-derived timestamps directly to frames.
%   2) Fill remaining NaNs using fixed 100 ms interval via fillNaNtimestamps.
%
% This mirrors the original hC timestamping logic.

if nargin < 3 || isempty(ledPrefix)
    ledPrefix = 'lime';
end

[trainField, pulseEdgeField, pulseTimeField, resolvedPrefix] = ...
    resolveLedFields(tbytDat, ledPrefix);

blockIds = getNumericFieldColumn(tbytDat, trainField);

nBlocks = size(sigCell, 2);

% Preallocate timestamp rows for every non-empty block.
for blk = 1:nBlocks
    if ~isempty(sigCell{1, blk})
        nFrames = size(sigCell{1, blk}, 2);
        sigCell{2, blk} = NaN(1, nFrames);
    else
        sigCell{2, blk} = [];
    end
end

for tr = 1:numel(tbytDat)

    blk = blockIds(tr);

    if isnan(blk) || blk < 1 || blk > nBlocks
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

    % Original-style frame range.
    frameIdxAll = pulseEdge(1):pulseEdge(end);

    % Conservative length matching.
    % No interpolation or stretching.
    nAssign = min(numel(frameIdxAll), numel(pulseTime));

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
end

% Fill all remaining NaNs using the fixed 100 ms frame clock.
for blk = 1:nBlocks

    if isempty(sigCell{1, blk})
        continue;
    end

    sigCell{2, blk} = sigCell{2, blk}(1, 1:size(sigCell{1, blk}, 2));

    sigCell{2, blk} = fillNaNtimestamps(sigCell{2, blk});

    fprintf('Timestamped block #%d/%d using %s LED fields and fixed 100 ms NaN fill.\n', ...
        blk, nBlocks, resolvedPrefix);
end

timingInfo = struct();
timingInfo.ledPrefixRequested = ledPrefix;
timingInfo.ledPrefixResolved = resolvedPrefix;
timingInfo.trainField = trainField;
timingInfo.pulseEdgeField = pulseEdgeField;
timingInfo.pulseTimeField = pulseTimeField;
timingInfo.timestampFillMethod = 'fillNaNtimestamps_fixed_100ms';

end


%% ========================================================================
function [trainField, pulseEdgeField, pulseTimeField, resolvedPrefix] = resolveLedFields(tbytDat, ledPrefix)
% resolveLedFields
%
% Resolves LED field names from a requested prefix.
%
% For ledPrefix = 'lime':
%   limeLEDTrainI
%   limeLEDPulsesOfTrain
%   limeLED
%
% For ledPrefix = 'blue':
%   blueLEDTrainI
%   blueLEDPulsesOfTrain
%   blueLED

if nargin < 2 || isempty(ledPrefix)
    ledPrefix = 'lime';
end

prefixCandidates = {};
rawCandidates = {ledPrefix, 'lime', 'blue'};

for i = 1:numel(rawCandidates)
    if isempty(rawCandidates{i})
        continue;
    end

    if ~any(strcmp(prefixCandidates, rawCandidates{i}))
        prefixCandidates{end+1} = rawCandidates{i}; %#ok<AGROW>
    end
end

for i = 1:numel(prefixCandidates)

    p = prefixCandidates{i};

    trainFieldTmp = sprintf('%sLEDTrainI', p);
    pulseEdgeFieldTmp = sprintf('%sLEDPulsesOfTrain', p);
    pulseTimeFieldTmp = sprintf('%sLED', p);

    if isfield(tbytDat, trainFieldTmp) && ...
            isfield(tbytDat, pulseEdgeFieldTmp) && ...
            isfield(tbytDat, pulseTimeFieldTmp)

        trainField = trainFieldTmp;
        pulseEdgeField = pulseEdgeFieldTmp;
        pulseTimeField = pulseTimeFieldTmp;
        resolvedPrefix = p;
        return;
    end
end

error(['Could not find compatible LED fields in tbytDat. ' ...
    'Tried prefixes: %s'], strjoin(prefixCandidates, ', '));

end


%% ========================================================================
function vals = getNumericFieldColumn(s, fieldName)
% getNumericFieldColumn
%
% Extracts a scalar numeric field from a struct array as a column vector.

vals = NaN(numel(s), 1);

for i = 1:numel(s)

    if ~isfield(s, fieldName)
        continue;
    end

    v = s(i).(fieldName);

    if iscell(v)
        v = cell2mat(v);
    end

    if isempty(v)
        continue;
    end

    vals(i) = double(v(1));
end

end


%% ========================================================================
function v = getNumericVector(x)
% getNumericVector
%
% Converts numeric or cell-wrapped numeric content to a row vector.

if iscell(x)
    x = cell2mat(x);
end

if isempty(x)
    v = [];
else
    v = double(x(:)');
end

end