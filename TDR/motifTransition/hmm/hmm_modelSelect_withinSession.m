function out = hmm_modelSelect_withinSession(seqC, Slist, opt)
%HMM_MODELSELECT_WITHINSESSION  Model-selection wrapper for diagonal-Gaussian HMM.
%
% out = hmm_modelSelect_withinSession(seqC, Slist, opt)
%
% INPUT
%   seqC  : 1xN cell, each cell is KxT_i or TxK (we will coerce to KxT)
%   Slist : vector of candidate # states, e.g. [2 4 6 8 10]
%   opt   : struct with fields (recommended defaults below)
%
% OUTPUT
%   out.table : table with per-S metrics (LL, BIC, occupancy, dwell, CV LL, etc.)
%   out.models: 1xnumel(Slist) cell of best models
%   out.details: struct array with hist, stats, etc.

% ---------------- defaults ----------------
if nargin < 2 || isempty(Slist), Slist = [2 4 6 8 10]; end
if nargin < 3, opt = struct(); end

opt = setDefaultOpt(opt);

% coerce seqC to KxT and build global stats
seqC = coerceSeqC(seqC);           % each is KxT
K = size(seqC{1},1);

% global mean/var for init
[Xall, ~] = stackSeqC(seqC);       % K x (sumT)
muGlob  = mean(Xall, 2, 'omitnan'); % Kx1
varGlob = var(Xall, 0, 2, 'omitnan');
varGlob(~isfinite(varGlob) | varGlob < opt.varFloor) = opt.varFloor;

% total # observations (time bins across trials)
Nobs = sum(cellfun(@(x) size(x,2), seqC));

% allocate
nS = numel(Slist);
models  = cell(1,nS);
details = repmat(struct(), 1, nS);

% optional CV partition (over sequences/trials)
if opt.doCV
    cvFoldIdx = makeKfoldIdx(numel(seqC), opt.cvK, opt.rng);
else
    cvFoldIdx = [];
end

rows = struct('S',[],'LL',[],'BIC',[],'AIC',[],'nParams',[], ...
              'minOcc',[],'occEntropy',[],'meanDwell',[],'medianDwell',[], ...
              'cvLL_perObs',[],'cvLL_sem',[]);

for ii = 1:nS
    S = Slist(ii);

    % fit on full data (best restart)
    [modelBest, bestLL, histBest] = hmm_fit_restarts_gauss_diag(seqC, S, opt, muGlob, varGlob);

    % post-fit stats on full data (for occupancy/dwell)
    stats = estep_gauss_diag(seqC, modelBest, opt);  % uses forward-backward
    [occ, dwell] = occupancy_and_dwell(stats, Nobs);

    % parameter count for diag-Gauss HMM
    nParams = hmm_numParams_diag(S, K);

    % BIC/AIC (using bestLL on full data)
    BIC = -2*bestLL + nParams*log(max(1, Nobs));
    AIC = -2*bestLL + 2*nParams;

    % cross-validated LL (per observation)
    if opt.doCV
        [cvMean, cvSem] = hmm_cv_loglik(seqC, S, opt, muGlob, varGlob, cvFoldIdx);
    else
        cvMean = NaN; cvSem = NaN;
    end

    % save
    models{ii} = modelBest;
    details(ii).hist = histBest;
    details(ii).stats = stats;
    details(ii).occ = occ;
    details(ii).dwell = dwell;
    details(ii).bestLL = bestLL;

    rows(ii).S = S;
    rows(ii).LL = bestLL;
    rows(ii).BIC = BIC;
    rows(ii).AIC = AIC;
    rows(ii).nParams = nParams;
    rows(ii).minOcc = min(occ);
    rows(ii).occEntropy = -nansum(occ(:) .* log(occ(:)+eps));
    rows(ii).meanDwell = mean(dwell,'omitnan');
    rows(ii).medianDwell = median(dwell,'omitnan');
    rows(ii).cvLL_perObs = cvMean;
    rows(ii).cvLL_sem = cvSem;

    if opt.verbose
        fprintf('S=%d | LL=%.3f | BIC=%.3f | minOcc=%.3g | meanDwell=%.3f | CVLL=%.4f\n', ...
            S, bestLL, BIC, rows(ii).minOcc, rows(ii).meanDwell, cvMean);
    end
end

T = struct2table(rows);

% ---------------- plots ----------------
if opt.doPlots
    plot_modelSelection_summary(T, opt);
end

out = struct();
out.table  = T;
out.models = models;
out.details = details;
out.opt = opt;

end

%% ============================ helpers ===================================

function opt = setDefaultOpt(opt)
opt = setFieldDefault(opt,'maxIters',200);
opt = setFieldDefault(opt,'nRestarts',10);
opt = setFieldDefault(opt,'tolRel',1e-4);
opt = setFieldDefault(opt,'sticky',0.90);     % strong self-transition init helps
opt = setFieldDefault(opt,'varFloor',1e-3);
opt = setFieldDefault(opt,'verbose',true);
opt = setFieldDefault(opt,'rng',1);

% CV
opt = setFieldDefault(opt,'doCV',true);
opt = setFieldDefault(opt,'cvK',5);
opt = setFieldDefault(opt,'cvRestarts',5);   % fewer to save time
opt = setFieldDefault(opt,'cvMaxIters',150);

% plotting
opt = setFieldDefault(opt,'doPlots',true);
end

function opt = setFieldDefault(opt, name, val)
if ~isfield(opt,name) || isempty(opt.(name))
    opt.(name) = val;
end
end

function seqC = coerceSeqC(seqC)
% Ensure each sequence is KxT (motifs x time)
for n = 1:numel(seqC)
    X = seqC{n};
    if isempty(X), continue; end
    if size(X,1) < size(X,2)
        % assume already KxT
    else
        % could be TxK
        if size(X,2) < size(X,1)
            X = X';
        end
    end
    seqC{n} = X;
end
end

function [Xall, idxC] = stackSeqC(seqC)
% Stack into K x sumT and keep indices per trial
K = size(seqC{1},1);
Ttot = sum(cellfun(@(x) size(x,2), seqC));
Xall = nan(K, Ttot);
idxC = cell(1,numel(seqC));
c = 1;
for n = 1:numel(seqC)
    Tn = size(seqC{n},2);
    Xall(:, c:c+Tn-1) = seqC{n};
    idxC{n} = c:c+Tn-1;
    c = c + Tn;
end
end

function nParams = hmm_numParams_diag(S, K)
% pi: S-1, A: S*(S-1), mu: S*K, var: S*K
nParams = (S-1) + S*(S-1) + 2*S*K;
end

function [modelBest, bestLL, histBest] = hmm_fit_restarts_gauss_diag(seqC, S, opt, muGlob, varGlob)
% Fit with multiple restarts, keep best LL. Uses EM with monotonic safeguard.

bestLL = -inf;
modelBest = [];
histBest = struct('LL',[],'restart',NaN);

for r = 1:opt.nRestarts
    [pi0, A0, mu0, var0] = initParams(seqC, S, opt, muGlob, varGlob);
    model = struct('pi',pi0,'A',A0,'mu',mu0,'var',var0);

    LL_hist = nan(opt.maxIters,1);
    prevModel = model;
    prevLL    = -inf;

    for it = 1:opt.maxIters
        stats = estep_gauss_diag(seqC, model, opt);
        LL = stats.LL;
        LL_hist(it) = LL;

        if opt.verbose && (it<=5 || mod(it,10)==0)
            if it==1
                fprintf('Restart %d | EM %3d | LL=%.3f\n', r, it, LL);
            else
                relImpro = (LL - prevLL) / max(1, abs(prevLL));
                fprintf('Restart %d | EM %3d | LL=%.3f | relImpro=%.3e\n', r, it, LL, relImpro);
            end
        end

        % monotonicity safeguard
        if it > 1 && LL < prevLL - 1e-8
            if opt.verbose
                fprintf('Restart %d | EM %3d | LL decreased (Δ=%.6g). Rollback & stop.\n', ...
                    r, it, LL - prevLL);
            end
            model = prevModel;
            LL_hist(it) = prevLL;
            break;
        end

        % convergence
        if it > 1
            relImpro = (LL - prevLL) / max(1, abs(prevLL));
            if relImpro >= 0 && relImpro < opt.tolRel
                break;
            end
        end

        % M-step
        prevModel = model;
        prevLL = LL;
        model = mstep_gauss_diag(stats, opt);

        assert(all(isfinite(model.pi)) && all(isfinite(model.A(:))) && ...
               all(isfinite(model.mu(:))) && all(isfinite(model.var(:))), ...
               'Non-finite parameter produced in M-step.');
    end

    finalLL = LL_hist(find(isfinite(LL_hist),1,'last'));
    if finalLL > bestLL
        bestLL = finalLL;
        modelBest = model;
        histBest = struct('LL', LL_hist, 'restart', r);
    end
end

end

function [pi0, A0, mu0, var0] = initParams(seqC, S, opt, muGlob, varGlob)
% Simple, stable init:
%   pi: random simplex
%   A : sticky self-transition + uniform off-diagonal
%   mu: perturb global mean
%   var: global var (floored)

K = numel(muGlob);

pi0 = normalizeDim(rand(S,1), 1);

if S==1
    A0 = 1;
else
    if opt.sticky > 0
        A0 = ones(S) * (1-opt.sticky)/(S-1);
        A0(1:S+1:end) = opt.sticky;
    else
        A0 = normalizeDim(rand(S,S), 2);
    end
end

% mu0: SxK
mu0 = repmat(muGlob(:)', S, 1) + 0.1*randn(S,K);

% var0: SxK
var0 = repmat(varGlob(:)', S, 1);
var0(var0 < opt.varFloor) = opt.varFloor;

end

function X = normalizeDim(X, dim)
s = sum(X, dim);
s(s==0) = 1;
X = X ./ s;
end

function stats = estep_gauss_diag(seqC, model, opt)
% Forward-backward for diagonal Gaussian emissions.
% Returns:
%   stats.LL: total log-likelihood across sequences
%   stats.gammaSum: Sx1 sum_t p(z_t=s|x)
%   stats.gamma1Sum: Sx1 sum over p(z_1=s|x)
%   stats.xiSum: SxS sum_t p(z_t=i,z_{t+1}=j|x)
%   stats.xSum: SxK sum_t gamma_s(t)*x_t
%   stats.x2Sum: SxK sum_t gamma_s(t)*x_t.^2

S = numel(model.pi);
K = size(model.mu,2);

gammaSum  = zeros(S,1);
gamma1Sum = zeros(S,1);
xiSum     = zeros(S,S);
xSum      = zeros(S,K);
x2Sum     = zeros(S,K);
LLtot     = 0;

for n = 1:numel(seqC)
    X = seqC{n};  % KxT
    if isempty(X), continue; end
    T = size(X,2);

    logB = log_emission_gauss_diag(X, model.mu, model.var); % SxT
    logpi = log(model.pi(:) + eps);
    logA  = log(model.A + eps);

    % forward with scaling
    logalpha = zeros(S,T);
    c = zeros(1,T);

    logalpha(:,1) = logpi + logB(:,1);
    [logalpha(:,1), c(1)] = logsumexp_normalize(logalpha(:,1));

    for t = 2:T
        tmp = logA' + logalpha(:,t-1); % SxS: (j,i)?? we want log sum_i alpha_i + A_i->j
        logalpha(:,t) = logB(:,t) + logsumexp(tmp, 2);
        [logalpha(:,t), c(t)] = logsumexp_normalize(logalpha(:,t));
    end

    % total loglik for this sequence
    LL = sum(c);
    LLtot = LLtot + LL;

    % backward
    logbeta = zeros(S,T);
    logbeta(:,T) = 0; % normalized space

    for t = T-1:-1:1
        tmp = logA + (logB(:,t+1) + logbeta(:,t+1))'; % SxS
        logbeta(:,t) = logsumexp(tmp, 2);
        % keep in same normalized space as alpha (subtract c(t+1))
        logbeta(:,t) = logbeta(:,t) - c(t+1);
    end

    % gamma
    loggamma = logalpha + logbeta;
    loggamma = loggamma - logsumexp(loggamma,1);
    gamma = exp(loggamma); % SxT

    gammaSum = gammaSum + sum(gamma,2);
    gamma1Sum = gamma1Sum + gamma(:,1);

    % xi sums
    for t = 1:T-1
        logxi = logalpha(:,t) + logA + (logB(:,t+1) + logbeta(:,t+1))';
        logxi = logxi - logsumexp(logxi(:),1);
        xiSum = xiSum + exp(logxi);
    end

    % emission stats
    for s = 1:S
        gs = gamma(s,:)'; % Tx1
        xSum(s,:)  = xSum(s,:)  + (gs' * X');           % 1xK
        x2Sum(s,:) = x2Sum(s,:) + (gs' * (X'.^2));      % 1xK
    end
end

stats = struct();
stats.LL = LLtot;
stats.gammaSum = gammaSum;
stats.gamma1Sum = gamma1Sum;
stats.xiSum = xiSum;
stats.xSum = xSum;
stats.x2Sum = x2Sum;

end

function model = mstep_gauss_diag(stats, opt)
S = numel(stats.gammaSum);
K = size(stats.xSum,2);

% pi
pi = stats.gamma1Sum;
pi = pi / max(eps, sum(pi));

% A
A = stats.xiSum;
A = A ./ max(eps, sum(A,2));

% mu, var
mu = zeros(S,K);
varr = zeros(S,K);

for s = 1:S
    w = max(eps, stats.gammaSum(s));
    mu(s,:) = stats.xSum(s,:) / w;
    Ex2 = stats.x2Sum(s,:) / w;
    varr(s,:) = Ex2 - mu(s,:).^2;
end

varr(~isfinite(varr) | varr < opt.varFloor) = opt.varFloor;

model = struct('pi',pi,'A',A,'mu',mu,'var',varr);

end

function logB = log_emission_gauss_diag(X, mu, varr)
% X: KxT, mu/varr: SxK -> logB: SxT
S = size(mu,1);
K = size(mu,2);
T = size(X,2);

logB = zeros(S,T);
Xk = X; % KxT

for s = 1:S
    m = mu(s,:)';       % Kx1
    v = varr(s,:)';     % Kx1
    % log N(x; m, diag(v)) summed over K
    % = -0.5 * [sum log(2πv) + sum (x-m)^2/v]
    const = -0.5 * sum(log(2*pi*v + eps));
    dif2 = (Xk - m).^2;          % KxT
    quad = -0.5 * sum(dif2 ./ (v + eps), 1); % 1xT
    logB(s,:) = const + quad;
end
end

function [vNorm, c] = logsumexp_normalize(v)
% subtract logsumexp to normalize in log space; return scaling term c
c = logsumexp(v,1);
vNorm = v - c;
end

function y = logsumexp(A, dim)
if nargin < 2, dim = 1; end
amax = max(A, [], dim);
y = amax + log(sum(exp(A - amax), dim) + eps);
end

function [occ, dwell] = occupancy_and_dwell(stats, Nobs)
% occupancy: expected fraction of time in each state
occ = stats.gammaSum / max(eps, sum(stats.gammaSum));

% dwell time estimate:
% segments_i ≈ expected #entries into i
% entries = gamma1(i) + sum_{t} sum_{j!=i} xi(j,i,t)
% We only have xi summed over time, so:
S = numel(stats.gammaSum);
entries = stats.gamma1Sum(:);
xi = stats.xiSum;
for i = 1:S
    entries(i) = entries(i) + sum(xi(:,i)) - xi(i,i);
end

timeInState = stats.gammaSum(:);       % expected # time bins in state i
dwell = timeInState ./ max(eps, entries);

% (Optional) you can also report the Markov implied dwell: 1/(1-Aii)
% dwell_markov = 1 ./ max(eps, 1 - diag(A));
end

function idx = makeKfoldIdx(N, K, rngSeed)
rng(rngSeed);
perm = randperm(N);
idx = zeros(N,1);
for i = 1:N
    idx(perm(i)) = mod(i-1, K) + 1;
end
end

function [cvMean, cvSem] = hmm_cv_loglik(seqC, S, opt, muGlob, varGlob, foldIdx)
% Fit on train folds, evaluate LL on test folds (forward only).
Kfold = opt.cvK;
ll_perObs = nan(Kfold,1);

% use cheaper opts for CV
optCV = opt;
optCV.nRestarts = opt.cvRestarts;
optCV.maxIters  = opt.cvMaxIters;
optCV.verbose   = false;

for k = 1:Kfold
    tr = find(foldIdx ~= k);
    te = find(foldIdx == k);

    [model, ~] = hmm_fit_restarts_gauss_diag(seqC(tr), S, optCV, muGlob, varGlob);

    % evaluate test LL (no EM)
    LLtest = hmm_loglik_gauss_diag(seqC(te), model);
    Ntest  = sum(cellfun(@(x) size(x,2), seqC(te)));
    ll_perObs(k) = LLtest / max(1, Ntest);
end

cvMean = mean(ll_perObs,'omitnan');
cvSem  = std(ll_perObs,'omitnan') / sqrt(sum(isfinite(ll_perObs)));

end

function LLtot = hmm_loglik_gauss_diag(seqC, model)
% Forward algorithm only, returns total LL across sequences.
S = numel(model.pi);
logpi = log(model.pi(:) + eps);
logA  = log(model.A + eps);

LLtot = 0;
for n = 1:numel(seqC)
    X = seqC{n}; % KxT
    if isempty(X), continue; end
    T = size(X,2);

    logB = log_emission_gauss_diag(X, model.mu, model.var); % SxT

    logalpha = zeros(S,T);
    c = zeros(1,T);

    logalpha(:,1) = logpi + logB(:,1);
    [logalpha(:,1), c(1)] = logsumexp_normalize(logalpha(:,1));

    for t = 2:T
        tmp = logA' + logalpha(:,t-1);
        logalpha(:,t) = logB(:,t) + logsumexp(tmp, 2);
        [logalpha(:,t), c(t)] = logsumexp_normalize(logalpha(:,t));
    end

    LLtot = LLtot + sum(c);
end
end

function plot_modelSelection_summary(T, opt)
figure('Color','w'); 
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

% BIC
nexttile; 
plot(T.S, T.BIC, '-o','LineWidth',1.5);
xlabel('# states S'); ylabel('BIC'); title('Model selection (BIC)'); grid on;

% CV LL
nexttile;
if any(isfinite(T.cvLL_perObs))
    errorbar(T.S, T.cvLL_perObs, T.cvLL_sem, '-o','LineWidth',1.5);
    xlabel('# states S'); ylabel('CV loglik / obs'); title('Generalization (CV)'); grid on;
else
    text(0.1,0.5,'CV disabled','Units','normalized');
    axis off;
end

% min occupancy
nexttile;
plot(T.S, T.minOcc, '-o','LineWidth',1.5);
xlabel('# states S'); ylabel('min occupancy'); title('Smallest state occupancy'); grid on;

% dwell
nexttile;
plot(T.S, T.meanDwell, '-o','LineWidth',1.5); hold on;
plot(T.S, T.medianDwell, '-s','LineWidth',1.5);
xlabel('# states S'); ylabel('dwell (bins)'); title('State dwell time'); legend({'mean','median'},'Box','off');
grid on;
end
