%% ========================================================================
function sigIntAlignedC = tbytSignalCellPsth_2row(sigCell, tbytDat, timeWin, timeStep, ledPrefix)
% tbytSignalCellPsth_2row
%
% Trial-aligns a 2-row block signal cell.
%
% New DAsC convention:
%   sigCell{1, blk} = projected DA, K x frames
%   sigCell{2, blk} = timestamps, 1 x frames
%
% Output:
%   sigIntAlignedC{1, tr}
%       K x aligned time bins signal matrix.
%
%   sigIntAlignedC{2, tr}
%       1 x aligned time bins time vector relative to tbytDat(tr).evtOn.
%
% Notes:
%   - Output uses a fixed time axis for all trials.
%   - Bins outside available data are left as NaN.
%   - No signal extrapolation is performed.

if nargin < 5 || isempty(ledPrefix)
    ledPrefix = 'lime';
end

eventField = 'evtOn';

[trainField, pulseEdgeField, ~, resolvedPrefix] = ...
    resolveLedFields(tbytDat, ledPrefix);

assert(isfield(tbytDat, eventField), ...
    'tbytDat must contain field %s.', eventField);

blockIds = getNumericFieldColumn(tbytDat, trainField);

tint = timeWin(1):timeStep:timeWin(2);
nTime = numel(tint);

sigIntAlignedC = cell(2, numel(tbytDat));

for tr = 1:numel(tbytDat)

    blk = blockIds(tr);

    if isnan(blk) || blk < 1 || blk > size(sigCell, 2)
        continue;
    end

    if isempty(sigCell{1, blk}) || isempty(sigCell{2, blk})
        continue;
    end

    sigMat = sigCell{1, blk};      % K x frames
    sigTime = sigCell{2, blk};     % 1 x frames

    K = size(sigMat, 1);

    % Preallocate full-length trial output.
    sigMatTrInt = NaN(K, nTime);
    sigIntAlignedC{2, tr} = tint;

    framesOfTrial = getNumericVector(tbytDat(tr).(pulseEdgeField));

    if isempty(framesOfTrial)
        sigIntAlignedC{1, tr} = sigMatTrInt;
        continue;
    end

    framesOfTrial = round(framesOfTrial(:)');

    f1 = max(1, framesOfTrial(1));
    f2 = min(size(sigMat, 2), framesOfTrial(end));

    if f2 <= f1
        sigIntAlignedC{1, tr} = sigMatTrInt;
        continue;
    end

    frameIdx = f1:f2;

    eventTime = getNumericScalar(tbytDat(tr).(eventField));

    if isnan(eventTime)
        sigIntAlignedC{1, tr} = sigMatTrInt;
        continue;
    end

    sigMatTr = sigMat(:, frameIdx);
    sigTimeTr = sigTime(1, frameIdx) - eventTime;

    % Keep only frames with valid timestamps.
    validTimeI = isfinite(sigTimeTr);

    if sum(validTimeI) < 2
        sigIntAlignedC{1, tr} = sigMatTrInt;
        continue;
    end

    sigMatTr = sigMatTr(:, validTimeI);
    sigTimeTr = sigTimeTr(validTimeI);

    % interp1 requires monotonic x values.
    [sigTimeTr, sortI] = sort(sigTimeTr, 'ascend');
    sigMatTr = sigMatTr(:, sortI);

    % Remove duplicate timestamps if any.
    [sigTimeTr, uniqueI] = unique(sigTimeTr, 'stable');
    sigMatTr = sigMatTr(:, uniqueI);

    if numel(sigTimeTr) < 2
        sigIntAlignedC{1, tr} = sigMatTrInt;
        continue;
    end

    % Only interpolate where the requested time axis is covered by the trial.
    tintI = sigTimeTr(1) <= tint & tint <= sigTimeTr(end);

    if ~any(tintI)
        sigIntAlignedC{1, tr} = sigMatTrInt;
        continue;
    end

    targetTime = tint(tintI);

    % Interpolate each motif-specific DA trace.
    % No extrapolation: unavailable bins remain NaN.
    for k = 1:K

        y = sigMatTr(k, :);
        validYI = isfinite(y) & isfinite(sigTimeTr);

        if sum(validYI) < 2
            continue;
        end

        sigMatTrInt(k, tintI) = interp1( ...
            sigTimeTr(validYI), ...
            y(validYI), ...
            targetTime, ...
            'linear');
    end

    sigIntAlignedC{1, tr} = sigMatTrInt;
end

fprintf('Aligned projected DA signals using %s LED fields: %d trials.\n', ...
    resolvedPrefix, numel(tbytDat));

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


%% ========================================================================
function x = getNumericScalar(x)
% getNumericScalar
%
% Converts numeric or cell-wrapped numeric content to a scalar.

if iscell(x)
    x = cell2mat(x);
end

if isempty(x)
    x = NaN;
else
    x = double(x(1));
end

end


%% ========================================================================
function t = fillNaNtimestamps_safe(t)
% fillNaNtimestamps_safe
%
% Uses your existing fillNaNtimestamps function if available.
% Otherwise falls back to linear interpolation/extrapolation.

if exist('fillNaNtimestamps', 'file') == 2
    try
        t = fillNaNtimestamps(t);
        return;
    catch ME
        warning('fillNaNtimestamps failed; using local fallback. Error: %s', ME.message);
    end
end

t = double(t(:)');

badI = ~isfinite(t);
goodI = isfinite(t);

if ~any(badI)
    return;
end

x = 1:numel(t);

if sum(goodI) >= 2
    t(badI) = interp1(x(goodI), t(goodI), x(badI), 'linear', 'extrap');

elseif sum(goodI) == 1
    % Not ideal, but prevents downstream crashes.
    t(badI) = t(goodI);

else
    warning('Timestamp vector contains no finite values. Leaving as NaN.');
end

end