function [newIdx, mapOrigToFinal, mapsByRound] = relabel_cluster_idx_multiround(oldIdx, rounds)
%RELABEL_CLUSTER_IDX_MULTIRound  Relabel cluster assignments after multi-round motif curation.
%
% [newIdx, mapOrigToFinal, mapsByRound] = relabel_cluster_idx_multiround(oldIdx, rounds)
%
% PURPOSE
%   Given an original cluster label vector oldIdx (e.g., 9237×1) whose labels refer
%   to the ORIGINAL motif set (e.g., 1..38), apply multiple curation rounds where
%   motifs are excluded and the remaining motifs are compactly relabeled to 1..K.
%   Excluded motifs are assigned NaN.
%
% INPUTS
%   oldIdx : numeric vector (any shape). Labels in {1..nOrig}. Can contain NaN.
%   rounds : struct array (1xR or Rx1) with fields:
%              - nBefore        : number of motifs BEFORE this round (e.g., 38 then 31)
%              - exclude        : vector of motif indices excluded in THIS round’s indexing
%
% EXAMPLE (your case)
%   rounds(1) = struct('nBefore',38,'exclude',[3 8 11 20 22 23 36]);  % after round1 -> 31
%   rounds(2) = struct('nBefore',31,'exclude',[24]);                 % after round2 -> 30
%   [idxNew, map38to30] = relabel_cluster_idx_multiround(cluster_idxC{1,2}, rounds);
%
% OUTPUTS
%   newIdx         : same size as oldIdx; relabeled to 1..nFinal; excluded => NaN
%   mapOrigToFinal : [nOrig x 1] mapping from original motif id -> final motif id (NaN if excluded)
%   mapsByRound    : cell{R,1}, each map is [nBefore x 1] giving round-local mapping (NaN if excluded)
%
% NOTES
%   - This assumes each round relabels by "compacting" remaining labels (i.e., setdiff order).
%   - Round 2 exclusions are interpreted in the indexing AFTER round 1 (as you described).
%
% Junchol Park lab-style: explicit + deterministic + easy to audit.

% -------------------- sanity --------------------
assert(isnumeric(oldIdx), 'oldIdx must be numeric.');
assert(isstruct(rounds) && all(isfield(rounds, {'nBefore','exclude'})), ...
    'rounds must be a struct array with fields nBefore and exclude.');

oldIdxVec = oldIdx(:);
if all(isnan(oldIdxVec))
    newIdx = oldIdx;  % nothing to do
    mapOrigToFinal = [];
    mapsByRound = {};
    return;
end

% infer nOrig from round(1)
nOrig = rounds(1).nBefore;
assert(isscalar(nOrig) && nOrig >= 1, 'rounds(1).nBefore must be a positive scalar.');

% -------------------- build per-round maps --------------------
R = numel(rounds);
mapsByRound = cell(R,1);

for r = 1:R
    nBefore = rounds(r).nBefore;
    excl = rounds(r).exclude;

    assert(isscalar(nBefore) && nBefore >= 1, 'round %d: nBefore must be positive scalar.', r);
    if isempty(excl)
        excl = [];
    end

    % validate exclude indices
    excl = excl(:)';
    assert(all(isfinite(excl)) && all(excl==round(excl)), 'round %d: exclude must be integer indices.', r);
    assert(all(excl>=1 & excl<=nBefore), 'round %d: exclude indices out of range 1..%d.', r, nBefore);

    keep = setdiff(1:nBefore, excl, 'stable');
    map = nan(nBefore,1);
    map(keep) = 1:numel(keep);   % compact relabel
    mapsByRound{r} = map;
end

% -------------------- compose maps: original -> final --------------------
% Start with identity mapping for original space
mapOrigToFinal = (1:nOrig)';

for r = 1:R
    mapR = mapsByRound{r};  % maps from labels in that round's "before" space -> after space (or NaN)

    % At round r, the current mapping values should be within 1..nBefore (or NaN)
    nBefore = rounds(r).nBefore;

    % Safety check: any non-NaN entry must be <= nBefore
    cur = mapOrigToFinal;
    ok = ~isnan(cur);
    assert(all(cur(ok) >= 1 & cur(ok) <= nBefore), ...
        'Round %d composition mismatch: current map has labels outside 1..%d.', r, nBefore);

    % Compose: orig -> (label before round r) -> (label after round r)
    newMap = nan(size(mapOrigToFinal));
    newMap(~ok) = NaN;
    newMap(ok) = mapR(cur(ok));   % if mapR returns NaN -> excluded in this round
    mapOrigToFinal = newMap;
end

% -------------------- apply to oldIdx --------------------
newIdxVec = nan(size(oldIdxVec));

valid = ~isnan(oldIdxVec);
ix = oldIdxVec(valid);

% ensure original indices are in range (allow zeros? no)
assert(all(ix>=1 & ix<=nOrig & ix==round(ix)), ...
    'oldIdx contains labels outside 1..%d or non-integers.', nOrig);

newIdxVec(valid) = mapOrigToFinal(ix);

% reshape back
newIdx = reshape(newIdxVec, size(oldIdx));

end
