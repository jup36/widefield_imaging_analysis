function [glmRez, found] = fetchByHeaderFromLookup_(hdr, uniqHdr, glmFlatFirst)
hdr = string(hdr);
j = find(uniqHdr == hdr, 1, 'first');
if isempty(j), glmRez = []; found = false; return; end
glmRez = glmFlatFirst{j};
found = ~isempty(glmRez) && isstruct(glmRez);
end