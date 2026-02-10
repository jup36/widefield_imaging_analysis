function idx = findHeaderInHeaderC(headerC, header)
%FINDHEADERINHEADERC  Find the [row, col] location of a session header in headerC.
%
% idx = findHeaderInHeaderC(headerC, header)
%
% INPUTS
%   headerC : cell array (J x S) of session headers (strings/chars or empty)
%   header  : char or string, e.g. 'm1045_122424'
%
% OUTPUT
%   idx     : [row, col] of the FIRST match
%             [] if no match is found
%
% NOTES
%   - Empty cells in headerC are ignored
%   - Exact string match is used
%   - If the header appears multiple times, the first occurrence
%     (row-major order) is returned

% -------------------- sanity --------------------
if nargin < 2 || isempty(headerC) || isempty(header)
    idx = [];
    return;
end

header = string(header);

% -------------------- flatten + filter empties --------------------
hdrFlat = headerC(:);

isNonEmpty = ~cellfun(@isempty, hdrFlat);
if ~any(isNonEmpty)
    idx = [];
    return;
end

hdrFlat = hdrFlat(isNonEmpty);
hdrFlatS = string(hdrFlat);

% -------------------- exact match --------------------
hit = find(hdrFlatS == header, 1, 'first');

if isempty(hit)
    idx = [];
    return;
end

% -------------------- map back to [row, col] --------------------
linIdxAll = find(isNonEmpty);      % linear indices into headerC
linIdx    = linIdxAll(hit);

[row, col] = ind2sub(size(headerC), linIdx);
idx = [row, col];

end
