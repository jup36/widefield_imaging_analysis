function reportExcluded(label, excludedMask, allAnimals)
if any(excludedMask)
    fprintf('WARNING: %s -- %d/%d animals excluded: %s\n', label, ...
        sum(excludedMask), numel(allAnimals), strjoin(allAnimals(excludedMask), ', '));
end
end