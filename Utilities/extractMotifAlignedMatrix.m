
function pethMat = extractMotifAlignedMatrix(alignedCell, trialMask, motifK)
% alignedCell{1, tr}: K x time
% trialMask: logical trial vector
% pethMat: trials x time

trialIdx = find(trialMask(:)');

pethMat = [];

for ii = 1:numel(trialIdx)

    tr = trialIdx(ii);

    if tr > size(alignedCell, 2)
        continue;
    end

    if isempty(alignedCell{1, tr})
        continue;
    end

    A = alignedCell{1, tr};

    if motifK > size(A, 1)
        continue;
    end

    pethMat = [pethMat; A(motifK, :)]; %#ok<AGROW>
end

end