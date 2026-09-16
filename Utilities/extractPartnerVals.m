function vals = extractPartnerVals(Manimal, targetMotif, partnerIdx, direction)
% Mean coupling between the target set and each partner class. Computed per
% target motif -- with that target removed from whichever class contains
% it, so the diagonal never enters -- then averaged across targets. Two
% stages: per-target, then across-target, matching the per-session /
% per-animal convention used throughout this project.
%
% Manimal(row, col): COLUMN leads, ROW follows.
%   'into' -> row trace    Manimal(t, :)   (target follows)
%   'from' -> column trace Manimal(:, t)   (target leads)
%
% Returns a 1 x nClasses row vector, in partnerIdx order.
nTargets = numel(targetMotif);
nClasses = numel(partnerIdx);
perTarget = nan(nTargets, nClasses);

for ti = 1:nTargets
    t = targetMotif(ti);

    switch lower(direction)
        case 'into', traceVec = Manimal(t, :);
        case 'from', traceVec = Manimal(:, t)';
        otherwise,   error('Unrecognized direction "%s".', direction);
    end

    for c = 1:nClasses
        idx = partnerIdx{c};
        idx = idx(idx ~= t);            % self-exclusion, per target
        if isempty(idx), continue; end  % class had only this target in it
        perTarget(ti, c) = mean(traceVec(idx), 'omitnan');
    end
end

vals = mean(perTarget, 1, 'omitnan');
end