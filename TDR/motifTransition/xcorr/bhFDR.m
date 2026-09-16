function q = bhFDR(pvals)
% BHFDR
%   Benjamini-Hochberg FDR correction. pvals: vector of raw p-values.
%   Returns q: same-size vector of FDR-corrected q-values (monotonic,
%   capped at 1).
pvals = pvals(:);
n = numel(pvals);
[pSorted, sortIdx] = sort(pvals);
qSorted = pSorted .* n ./ (1:n)';
% enforce monotonicity (q must be non-decreasing when unsorted back by
% descending rank) -- standard BH step-up correction
qSorted = flipud(cummin(flipud(qSorted)));
qSorted = min(qSorted, 1);

q = nan(n, 1);
q(sortIdx) = qSorted;
end
