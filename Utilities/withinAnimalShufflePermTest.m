function [r_obs, pVal, nullDist, nAnimalsUsed, nSessUsed] = withinAnimalShufflePermTest(dprimeVec, proxyVec, mIdVec, nShuffle)
% WITHINANIMALSHUFFLEPERMTEST
%   Pooled Pearson correlation between dprimeVec and proxyVec (across ALL
%   sessions from ALL animals), tested against a null built by
%   independently permuting each animal's OWN d' values among its OWN
%   sessions (proxy values held fixed), repeated nShuffle times. See
%   SECTION 3's header comment above for the full rationale -- this
%   preserves between-animal differences identically in the observed
%   statistic and every null draw, isolating within-animal (session-to-
%   session) covariation specifically.
%
%   dprimeVec, proxyVec : column vectors, already filtered to finite pairs
%   mIdVec              : cell array of animal ID strings, same length
%   nShuffle             : number of permutation draws
%
%   r_obs        : observed pooled Pearson correlation coefficient
%   pVal         : two-sided permutation p-value
%   nullDist     : [nShuffle x 1] null distribution of pooled r
%   nAnimalsUsed : number of unique animals contributing at least one session
%   nSessUsed    : total number of sessions used

dprimeVec = dprimeVec(:);
proxyVec  = proxyVec(:);
mIdVec    = mIdVec(:);

nSessUsed = numel(dprimeVec);
uAnimals  = unique(mIdVec, 'stable');
nAnimalsUsed = numel(uAnimals);

Robs = corrcoef(dprimeVec, proxyVec);
r_obs = Robs(1, 2);

% Pre-index each animal's rows once (not per shuffle) for efficiency
animalRowIdx = cell(nAnimalsUsed, 1);
for a = 1:nAnimalsUsed
    animalRowIdx{a} = find(strcmp(mIdVec, uAnimals{a}));
end

nullDist = nan(nShuffle, 1);
for sh = 1:nShuffle
    dShuf = dprimeVec;
    for a = 1:nAnimalsUsed
        rows = animalRowIdx{a};
        if numel(rows) < 2, continue; end   % single-session animal: nothing to permute, stays fixed
        permOrder = rows(randperm(numel(rows)));
        dShuf(rows) = dprimeVec(permOrder);
    end
    Rshuf = corrcoef(dShuf, proxyVec);
    nullDist(sh) = Rshuf(1, 2);
end

muNull = mean(nullDist, 'omitnan');
devObs  = abs(r_obs - muNull);
devNull = abs(nullDist - muNull);
pVal = (1 + sum(devNull >= devObs, 'omitnan')) / (nShuffle + 1);
end