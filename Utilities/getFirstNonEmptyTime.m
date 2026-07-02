function timeX = getFirstNonEmptyTime(alignedCell)

timeX = [];

for tr = 1:size(alignedCell, 2)

    if size(alignedCell, 1) >= 2 && ~isempty(alignedCell{2, tr})
        timeX = alignedCell{2, tr};
        return;
    end
end

end