function validJ = findValidSessions(xcorrRezC, rowIdx, extractFn)
% Returns the column indices of xcorrRezC(rowIdx, :) that hold a usable
% (non-empty struct, extractable non-empty matrix) session result. Used so
% the early/late split is defined over USABLE sessions rather than raw
% column indices.
validJ = [];
nSessions = size(xcorrRezC, 2);
for j = 1:nSessions
    entry = xcorrRezC{rowIdx, j};
    if isempty(entry) || ~isstruct(entry)
        continue;
    end
    try
        M = extractFn(entry);
    catch
        continue;   % missing field
    end
    if isempty(M)
        continue;
    end
    validJ(end+1) = j;   %#ok<AGROW>
end
end