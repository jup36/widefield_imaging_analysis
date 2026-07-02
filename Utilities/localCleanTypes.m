function T = localCleanTypes(T)
% Convert grouping variables to categorical.

if isempty(T)
    return
end

if ismember('mouseID', T.Properties.VariableNames)
    T.mouseID = categorical(T.mouseID);
end

if ismember('sessionID', T.Properties.VariableNames)
    T.sessionID = categorical(T.sessionID);
end

if ismember('sourceLabel', T.Properties.VariableNames)
    T.sourceLabel = categorical(T.sourceLabel);
end

end
