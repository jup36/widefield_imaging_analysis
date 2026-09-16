function binMap = extractBinIndexMap(X_names, prefix)
% EXTRACTBININDEXMAP
%   Finds columns in X_names starting with prefix (e.g. 'toneOnGo_rc'),
%   parses the trailing numeric lag index from each matching name, and
%   returns a containers.Map from that numeric bin index -> column
%   index in X_names. Used to pair Go/NoGo columns by their ACTUAL lag
%   number rather than by position, since the GLM's sparse-column drop
%   (finiteFrac<0.25) can remove a column from one side without removing
%   the matching column on the other side in a given session.
binMap = containers.Map('KeyType', 'double', 'ValueType', 'double');
colIdx = find(startsWith(X_names, prefix));
for ii = 1:numel(colIdx)
    tok = regexp(X_names{colIdx(ii)}, ['^' regexptranslate('escape', prefix) '(\d+)'], 'tokens', 'once');
    if ~isempty(tok)
        b = str2double(tok{1});
        binMap(b) = colIdx(ii);
    end
end
end