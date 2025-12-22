function Xall = catSeq(seqC, layout)
% returns [K x Tall]
Xall = [];
for i = 1:numel(seqC)
    Xi = getSeq(seqC{i}, layout);
    Xall = [Xall, Xi]; %#ok<AGROW>
end
end
