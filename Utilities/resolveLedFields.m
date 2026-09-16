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
    ledPrefix = 'blue';
end
prefixCandidates = {ledPrefix}; % Require the configured acquisition channel.
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