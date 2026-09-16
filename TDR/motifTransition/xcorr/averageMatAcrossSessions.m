function Manimal = averageMatAcrossSessions(xcorrRezC, rowIdx, extractFn)
% Generalized version of the earlier averageXcorrAcrossSessions: extractFn
% is a function handle taking one session's result struct and returning a
% [K x K] matrix -- either the observed matrix (@(e) e.obs.(fieldName)) or
% one shuffle slice (@(e) e.shuf.(fieldName)(:,:,s)). Using one function
% for both keeps the observed and null computations mechanically
% identical, which is the whole point: the null must be built by the exact
% same averaging procedure as the observed statistic.
stack = [];
nSessions = size(xcorrRezC, 2);
for j = 1:nSessions
    entry = xcorrRezC{rowIdx, j};
    if isempty(entry) || ~isstruct(entry)
        continue;
    end
    try
        M = extractFn(entry);
    catch
        continue;   % missing field (e.g. .shuf absent if doTimeShuffle was false) -- skip
    end
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
