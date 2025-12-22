function [pi0, A0, mu0, var0] = initParams(seqC, S, opt, muGlob, varGlob)
K = numel(muGlob);

pi0 = normalizeDim(rand(S,1), 1);

% A init
if opt.sticky > 0
    if S==1
        A0 = 1;
    else
        A0 = ones(S) * (1-opt.sticky)/(S-1);
        A0(1:S+1:end) = opt.sticky;
    end
else
    A0 = normalizeDim(rand(S,S), 2);
end

% emission init around global stats
mu0  = repmat(muGlob,1,S) + 0.1*randn(K,S);
var0 = repmat(varGlob,1,S);

% optional kmeans on sampled observations
if isfield(opt,'init') && strcmpi(opt.init,'kmeans')
    Xall = catSeq(seqC, opt.dataLayout)'; % [Tall x K]
    nSamp = min(size(Xall,1), 5000);
    idx = randperm(size(Xall,1), nSamp);
    Xs = Xall(idx,:);

    % guard against all-NaN rows
    good = all(isfinite(Xs),2);
    Xs = Xs(good,:);
    if size(Xs,1) >= S
        lab = kmeans(Xs, S, 'Replicates', 3, 'MaxIter', 200);
        for s = 1:S
            mu0(:,s) = mean(Xs(lab==s,:), 1, 'omitnan')';
            vv = var(Xs(lab==s,:), 0, 1, 'omitnan')';
            vv(~isfinite(vv) | vv < opt.varFloor) = opt.varFloor;
            var0(:,s) = vv;
        end
    end
end
end