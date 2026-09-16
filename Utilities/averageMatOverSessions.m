function Manimal = averageMatOverSessions(xcorrRezC, rowIdx, sessionJ, extractFn)
% Mean matrix over the given sessions. Sessions are weighted equally, not
% by trial count, consistent with the rest of this project.
stack = [];
for k = 1:numel(sessionJ)
    entry = xcorrRezC{rowIdx, sessionJ(k)};
    if ~isUsableEntry(entry, extractFn), continue; end
    stack = cat(3, stack, extractFn(entry));
end

if isempty(stack)
    Manimal = [];
else
    Manimal = mean(stack, 3, 'omitnan');
end
end