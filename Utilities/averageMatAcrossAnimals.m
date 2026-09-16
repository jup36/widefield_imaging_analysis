function Mgroup = averageMatAcrossAnimals(xcorrRezC, mIdC, animalList, extractFn)
% Generalized version of the earlier averageXcorrAcrossAnimals, using the
% same extractFn pattern as averageMatAcrossSessions above.
mIdCol = mIdC(:,1);

animalMats = {};
for a = 1:numel(animalList)
    thisID = animalList{a};
    rowIdx = find(strcmp(mIdCol, thisID), 1, 'first');
    if isempty(rowIdx)
        continue;
    end

    Manimal = averageMatAcrossSessions(xcorrRezC, rowIdx, extractFn);
    if isempty(Manimal)
        continue;
    end
    animalMats{end+1} = Manimal; %#ok<AGROW>
end

assert(~isempty(animalMats), 'averageMatAcrossAnimals: no animals had valid data.');
Mgroup = mean(cat(3, animalMats{:}), 3, 'omitnan');
end