function [grandMean, rowMean] = meanNestedCellArray(C)

rowMean = nan(size(C,1), 1);

for r = 1:size(C,1)

    rowVals = [];

    for c = 1:size(C,2)

        thisEntry = C{r,c};

        if isempty(thisEntry)
            continue
        end

        if iscell(thisEntry)
            for k = 1:numel(thisEntry)
                if ~isempty(thisEntry{k}) && isnumeric(thisEntry{k})
                    rowVals = [rowVals; thisEntry{k}(:)];
                end
            end
        elseif isnumeric(thisEntry)
            rowVals = [rowVals; thisEntry(:)];
        end
    end

    rowVals = rowVals(isfinite(rowVals));

    if ~isempty(rowVals)
        rowMean(r) = mean(rowVals, 'omitnan');
    end
end

grandMean = mean(rowMean, 'omitnan');

end