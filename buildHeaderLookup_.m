function [uniqHdr, idxFirst, glmFlatFirst] = buildHeaderLookup_(headerC, glmRezC)
% Flatten headerC/glmRezC and build unique header list + first index mapping.

hdrFlat = headerC(:);
glmFlat = glmRezC(:);

isHdr = ~cellfun(@isempty, hdrFlat);
hdrFlat = hdrFlat(isHdr);
glmFlat = glmFlat(isHdr);

hdrFlatS = string(hdrFlat);

[uniqHdr, ~, ic] = unique(hdrFlatS, 'stable');
if numel(uniqHdr) < numel(hdrFlatS)
    counts = accumarray(ic, 1);
    dup = uniqHdr(counts > 1);
    warning('collectPerMouseAndGlobalGlmPrjAxes:DuplicateHeaders', ...
        'Duplicate headers detected in headerC (matching uses first occurrence). Example(s): %s', ...
        strjoin(cellstr(dup(1:min(5,end))), ', '));
end

idxFirst = zeros(numel(uniqHdr),1);
glmFlatFirst = cell(numel(uniqHdr),1);
for i = 1:numel(uniqHdr)
    ii = find(hdrFlatS == uniqHdr(i), 1, 'first');
    idxFirst(i) = ii;
    glmFlatFirst{i} = glmFlat{ii};
end
end