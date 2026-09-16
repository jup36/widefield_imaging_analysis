function [Go_group, NoGo_group] = twoStageGroupAverage(sessGoAll, sessNoGoAll, nAnimals, K)
% Two-stage average (session -> animal -> group, equal animal weighting),
% shared by both the observed computation and every permutation draw, so
% the null is built by the EXACT same aggregation procedure as the
% observed statistic (same principle as the xcorr pipeline's group-level
% null).
animalGo   = nan(nAnimals, K);
animalNoGo = nan(nAnimals, K);
for a = 1:nAnimals
    Gmat = sessGoAll{a};
    Nmat = sessNoGoAll{a};
    if isempty(Gmat) || all(isnan(Gmat(:)))
        continue;   % animal excluded (zero usable sessions)
    end
    animalGo(a, :)   = mean(Gmat, 1, 'omitnan');
    animalNoGo(a, :) = mean(Nmat, 1, 'omitnan');
end
Go_group   = mean(animalGo,   1, 'omitnan');
NoGo_group = mean(animalNoGo, 1, 'omitnan');
end

