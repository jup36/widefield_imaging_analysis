function Mgroup = averageXcorrAcrossAnimals(xcorrRezC, mIdC, animalList, fieldName)
% Two-stage average (matches the pipeline's established convention):
%   1) within each animal, average across that animal's valid sessions
%   2) across animals in animalList, average the resulting animal-level
%      matrices (equal animal weighting, regardless of session count)
mIdCol = mIdC(:,1);

animalMats = {};
for a = 1:numel(animalList)
    thisID = animalList{a};
    rowIdx = find(strcmp(mIdCol, thisID), 1, 'first');
    if isempty(rowIdx)
        warning('averageXcorrAcrossAnimals:animalNotFound', ...
            'Animal "%s" not found in mIdC -- skipping.', thisID);
        continue;
    end

    Manimal = averageXcorrAcrossSessions(xcorrRezC, rowIdx, fieldName);
    if isempty(Manimal)
        warning('averageXcorrAcrossAnimals:noValidSessions', ...
            'Animal "%s" has no valid sessions for field "%s" -- skipping.', thisID, fieldName);
        continue;
    end
    animalMats{end+1} = Manimal; %#ok<AGROW>
end

assert(~isempty(animalMats), ...
    'averageXcorrAcrossAnimals:noValidAnimals', ...
    'No animals in the requested list had valid data for field "%s".', fieldName);

Mgroup = mean(cat(3, animalMats{:}), 3, 'omitnan');
end