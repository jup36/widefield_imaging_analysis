function Diff_group = twoStageGroupAverage3D(sessDiffAll, nAnimals, K, nBins)
% 3D analog of twoStageGroupAverage, for the BINWISE DIFFERENCE
% definition: sessDiffAll{a} is [nSessionsPerRow x K x nBins]. Averages
% across sessions (dim 1) within each animal, then across animals (equal
% weighting), preserving the bin dimension throughout -- same two-stage
% philosophy, just carrying an extra dimension.
animalDiff = nan(nAnimals, K, nBins);
for a = 1:nAnimals
    Dmat = sessDiffAll{a};
    if isempty(Dmat) || all(isnan(Dmat(:)))
        continue;
    end
    animalDiff(a, :, :) = mean(Dmat, 1, 'omitnan');
end
Diff_group = squeeze(mean(animalDiff, 1, 'omitnan'));   % [K x nBins]
end