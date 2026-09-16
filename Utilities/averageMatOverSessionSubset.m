function Manimal = averageMatOverSessionSubset(xcorrRezC, rowIdx, sessionJ, extractFn)
% Same averaging mechanics as averageMatAcrossSessions, but restricted to
% the supplied subset of column indices. Sessions are weighted equally
% (not by trial count), consistent with the rest of this project.
stack = [];
for k = 1:numel(sessionJ)
    entry = xcorrRezC{rowIdx, sessionJ(k)};
    if isempty(entry) || ~isstruct(entry)
        continue;
    end
    try
        M = extractFn(entry);
    catch
        continue;
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
