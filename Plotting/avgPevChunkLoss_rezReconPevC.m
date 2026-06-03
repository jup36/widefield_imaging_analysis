function [pevLoss_final, pevLoss_mouseSessAvg, pevLoss_sessC] = avgPevChunkLoss_rezReconPevC(rezReconPevC)
% avgPevChunkLoss_rezReconPevC
%
% Computes PEV loss averaged in the following order:
%   1) across chunks within each session
%   2) across sessions within each mouse
%   3) across mice
%
% INPUT
%   rezReconPevC : nMice x nSessions cell array
%                  Each non-empty cell should contain a struct with field:
%                  .pevChunkLoss [nMotifs x nChunks]
%
% OUTPUT
%   pevLoss_final        : [nMotifs x 1] final average across mice
%   pevLoss_mouseSessAvg : [nMotifs x nMice] mouse-level averages
%   pevLoss_sessC        : nMice x 1 cell array.
%                          Each cell is [nMotifs x nValidSessions] session-level averages

% ---------------------------------------------------------
% Get dimensions
% ---------------------------------------------------------
[nMice, nSessMax] = size(rezReconPevC);

% ---------------------------------------------------------
% Find number of motifs from the first valid entry
% ---------------------------------------------------------
nMotifs = [];

for iM = 1:nMice
    for iS = 1:nSessMax
        if ~isempty(rezReconPevC{iM, iS}) && ...
                isstruct(rezReconPevC{iM, iS}) && ...
                isfield(rezReconPevC{iM, iS}, 'pevChunkLoss') && ...
                ~isempty(rezReconPevC{iM, iS}.pevChunkLoss)

            nMotifs = size(rezReconPevC{iM, iS}.pevChunkLoss, 1);
            break
        end
    end

    if ~isempty(nMotifs)
        break
    end
end

if isempty(nMotifs)
    error('No valid pevChunkLoss field found in rezReconPevC.');
end

% ---------------------------------------------------------
% Preallocate
% ---------------------------------------------------------
pevLoss_mouseSessAvg = nan(nMotifs, nMice);
pevLoss_sessC = cell(nMice, 1);

% ---------------------------------------------------------
% Loop through mice and sessions
% ---------------------------------------------------------
for iM = 1:nMice

    pevLoss_sess = [];

    for iS = 1:nSessMax

        % Skip empty entries
        if isempty(rezReconPevC{iM, iS})
            continue
        end

        % Skip invalid entries
        if ~isstruct(rezReconPevC{iM, iS}) || ...
                ~isfield(rezReconPevC{iM, iS}, 'pevChunkLoss') || ...
                isempty(rezReconPevC{iM, iS}.pevChunkLoss)
            continue
        end

        pevChunkLoss = rezReconPevC{iM, iS}.pevChunkLoss;  % [nMotifs x nChunks]

        % Safety check
        if size(pevChunkLoss, 1) ~= nMotifs
            error('Motif number mismatch at mouse %d, session %d.', iM, iS);
        end

        % Average across chunks
        pevLoss_thisSess = mean(pevChunkLoss, 2, 'omitnan');  % [nMotifs x 1]

        % Append session vector
        pevLoss_sess = [pevLoss_sess, pevLoss_thisSess]; %#ok<AGROW>

    end

    % Store session-level vectors for this mouse
    pevLoss_sessC{iM} = pevLoss_sess;

    % Average across sessions within this mouse
    if ~isempty(pevLoss_sess)
        pevLoss_mouseSessAvg(:, iM) = mean(pevLoss_sess, 2, 'omitnan');
    end

end

% ---------------------------------------------------------
% Average across mice
% ---------------------------------------------------------
pevLoss_final = mean(pevLoss_mouseSessAvg, 2, 'omitnan');  % [nMotifs x 1]

end