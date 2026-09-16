function [goVal, nogoVal] = extractGoNogoVals(Manimal, targetMotif, goIdx, nogoIdx, direction)
% Per-target extraction with PER-TARGET self-exclusion, then average across
% targets (two-stage average). Factored out so the pooled pass and the
% per-stage passes share one implementation.
%
% Direction convention: Manimal(row, col) -- COLUMN leads, ROW follows.
%   'into' : target = follower -> row trace Manimal(t, :)
%   'from' : target = leader   -> column trace Manimal(:, t)
%
% The diagonal never enters: t is removed from whichever preference set
% contains it, and the two sets are disjoint by assertion upstream.
nTargets = numel(targetMotif);
goVal_perTarget   = nan(1, nTargets);
nogoVal_perTarget = nan(1, nTargets);

for ti = 1:nTargets
    t = targetMotif(ti);

    goIdx_noSelf_t   = goIdx(goIdx ~= t);
    nogoIdx_noSelf_t = nogoIdx(nogoIdx ~= t);

    switch lower(direction)
        case 'into'
            traceVec = Manimal(t, :);    % target = follower
        case 'from'
            traceVec = Manimal(:, t)';   % target = leader
        otherwise
            error('Unrecognized direction "%s".', direction);
    end

    goVal_perTarget(ti)   = mean(traceVec(goIdx_noSelf_t),   'omitnan');
    nogoVal_perTarget(ti) = mean(traceVec(nogoIdx_noSelf_t), 'omitnan');
end

goVal   = mean(goVal_perTarget,   'omitnan');
nogoVal = mean(nogoVal_perTarget, 'omitnan');
end
