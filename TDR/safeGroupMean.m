% ---- helper function (put at end of file or in its own file) ----
function m = safeGroupMean(B, idx)
    nRows = size(B, 1);
    validIdx = idx(idx <= nRows);      % keep only valid rows
    if isempty(validIdx)
        % no valid predictors for this group in this session
        m = nan(1, size(B, 2));        % 1 x nTime (or nMotifs) of NaNs
    else
        m = mean(B(validIdx, :), 1);   % average over existing rows only
    end
end