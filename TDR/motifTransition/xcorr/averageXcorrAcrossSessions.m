function Manimal = averageXcorrAcrossSessions(xcorrRezC, rowIdx, fieldName)
% Average one animal's (rowIdx) per-session [K x K] xcorr matrices across
% its valid (non-empty) sessions, ignoring NaNs (e.g. sessions with zero
% Go or NoGo trials, which return an all-NaN matrix for that field).
stack = [];
nSessions = size(xcorrRezC, 2);
for j = 1:nSessions
    entry = xcorrRezC{rowIdx, j};
    if isempty(entry) || ~isstruct(entry) || ~isfield(entry, 'obs') || ~isfield(entry.obs, fieldName)
        continue;
    end
    M = entry.obs.(fieldName);
    if isempty(M)
        continue;
    end
    stack = cat(3, stack, M);
end

if isempty(stack)
    Manimal = [];
else
    Manimal = mean(stack, 3, 'omitnan');
end
end
