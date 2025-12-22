function out = hmm_modelSelect_withinSession_cv(seqC, varargin)
% hmm_modelSelect_withinSession_cv
% -------------------------------------------------------------------------
% Drop-in, within-session diagonal-Gaussian HMM model selection with K-fold CV.
%
% What this does (synopsis)
%   1) Takes trial-by-trial data in a cell array seqC, where each trial is a
%      matrix of motif activity with shape K×T (K features/motifs, T time bins).
%      (If a trial is T×K, this code will automatically transpose it safely.)
%
%   2) For each candidate number of states S in opt.S_list:
%        - Runs K-fold cross-validation over trials:
%            * Fit HMM on training trials using EM (Baum–Welch), with multiple
%              random restarts (keeps best LL).
%            * Evaluate held-out (test) trials with the fitted model
%              (NO M-step on test; only computes log-likelihood).
%        - Returns per-fold and pooled per-observation test LL.
%
%   3) Optionally fits a final model on all trials at the selected S.
%
% Key design choices & safeguards (important)
%   - NO CV leakage: global mean/variance initializers are computed from
%     TRAINING folds only (not from all data).
%   - Stable log-space forward/backward with robust logsumexp handling -Inf.
%   - Diagonal Gaussian emissions; variance floor prevents collapse.
%   - Reproducibility: one rng seed controls folds + restarts deterministically.
%   - Safe orientation: each trial is forced to K×T, with K inferred robustly.
%
% Inputs
%   seqC : 1×N cell array. Each cell is K×T or T×K numeric.
%
% Name-Value options
%   'S_list'        : vector of candidate state counts (default 2:12)
%   'Kfold'         : number of folds (default 5)
%   'maxIter'       : EM iterations (default 200)
%   'tolLL'         : relative LL improvement tolerance (default 1e-6)
%   'nRestarts'     : random restarts per fold (default 5)
%   'sticky'        : sticky transition weight (default 0.0)
%   'varFloor'      : minimum variance per feature (default 1e-4)
%   'transPrior'    : Dirichlet pseudocount for A (default 1e-2)
%   'piPrior'       : Dirichlet pseudocount for pi (default 1e-2)
%   'rng'           : rng seed (default 0)
%   'verbose'       : true/false (default true)
%
% Output (struct) out
%   .opt                 : options used
%   .S_list              : candidate state counts
%   .cv                  : struct with fields per S:
%       .foldLL_perObs   : Kfold×1 per-fold test LL per observation
%       .foldNobs        : Kfold×1 number of test observations
%       .meanFoldLL      : mean over folds of per-fold perObs LL
%       .pooledLL        : pooled test LL / pooled Nobs
%   .bestS_by_pooledLL   : S with max pooledLL
%   .bestS_by_meanFoldLL : S with max meanFoldLL
%   .finalModel          : (optional) HMM fit on all data at bestS_by_pooledLL
%
% -------------------------------------------------------------------------
% NOTE: Observation model is diagonal Gaussian: p(x_t | z_t=s) = N(mu_s, diag(sig2_s)).
% -------------------------------------------------------------------------

% --------------------------- parse options ------------------------------
opt = struct();
opt.S_list     = 2:10;
opt.Kfold      = 5;
opt.maxIter    = 200;
opt.tolLL      = 1e-6;
opt.nRestarts  = 5;
opt.sticky     = 0.0;
opt.varFloor   = 1e-4;
opt.transPrior = 1e-2;
opt.piPrior    = 1e-2;
opt.rng        = 0;
opt.verbose    = true;
opt.fitFinalOnAll = true;

opt = parseNameValue(opt, varargin{:});

% Reproducibility for folds + restarts
rng(opt.rng);

% --------------------------- validate & coerce ---------------------------
seqC = coerceSeqC_KxT(seqC);

Ntr = numel(seqC);
K   = size(seqC{1}, 1);

% --------------------------- build folds --------------------------------
foldIdx = makeKfoldIdx(Ntr, opt.Kfold, opt.rng);

out = struct();
out.opt    = opt;
out.S_list = opt.S_list;
out.cv     = struct();

if opt.verbose
    fprintf('[HMM-CV] Trials=%d, Features(K)=%d, Kfold=%d, Restarts=%d\n', ...
        Ntr, K, opt.Kfold, opt.nRestarts);
end

% --------------------------- CV over S ----------------------------------
for sIdx = 1:numel(opt.S_list)
    S = opt.S_list(sIdx);

    foldLL_perObs = nan(opt.Kfold,1);
    foldNobs      = nan(opt.Kfold,1);

    if opt.verbose
        fprintf('\n[HMM-CV] S=%d\n', S);
    end

    for kf = 1:opt.Kfold
        testI  = (foldIdx == kf);
        trainI = ~testI;

        trainSeq = seqC(trainI);
        testSeq  = seqC(testI);

        % TRAIN-only global stats for initialization (prevents CV leakage)
        [Xall_train, ~] = stackSeqC(trainSeq);
        muGlob  = mean(Xall_train, 2, 'omitnan');
        varGlob = var(Xall_train, 0, 2, 'omitnan');
        varGlob = max(varGlob, opt.varFloor);

        % Fit best model over restarts on training set
        [modelBest, llTrainBest, histBest] = hmm_fit_restarts_gauss_diag( ...
            trainSeq, S, muGlob, varGlob, opt);

        % Evaluate held-out test LL (no parameter update)
        [llTest, nObsTest] = hmm_loglik_dataset(testSeq, modelBest);

        foldLL_perObs(kf) = llTest / max(1, nObsTest);
        foldNobs(kf)      = nObsTest;

        if opt.verbose
            fprintf('  Fold %d/%d: trainLL=%.3e, testLL/obs=%.6f (Nobs=%d)\n', ...
                kf, opt.Kfold, llTrainBest, foldLL_perObs(kf), nObsTest);
        end

        %#ok<NASGU> histBest; % keep if you want fold-by-fold EM curves
    end

    pooledLL = nansum(foldLL_perObs .* foldNobs) / max(1, nansum(foldNobs));
    meanFold = mean(foldLL_perObs, 'omitnan');

    out.cv(sIdx).S              = S;
    out.cv(sIdx).foldLL_perObs  = foldLL_perObs;
    out.cv(sIdx).foldNobs       = foldNobs;
    out.cv(sIdx).meanFoldLL     = meanFold;
    out.cv(sIdx).pooledLL       = pooledLL;

    if opt.verbose
        fprintf('  -> S=%d: meanFoldLL=%.6f, pooledLL=%.6f\n', S, meanFold, pooledLL);
    end
end

% --------------------------- choose best S -------------------------------
pooledLLs = arrayfun(@(c) c.pooledLL, out.cv);
meanLLs   = arrayfun(@(c) c.meanFoldLL, out.cv);

[~, iP] = max(pooledLLs);
[~, iM] = max(meanLLs);

out.bestS_by_pooledLL   = out.cv(iP).S;
out.bestS_by_meanFoldLL = out.cv(iM).S;

if opt.verbose
    fprintf('\n[HMM-CV] bestS_by_pooledLL=%d, bestS_by_meanFoldLL=%d\n', ...
        out.bestS_by_pooledLL, out.bestS_by_meanFoldLL);
end

% --------------------------- optional: final model on all data -----------
out.finalModel = [];
out.finalTrainLL = [];
out.finalHist = [];

if opt.fitFinalOnAll
    Sfinal = out.bestS_by_pooledLL;

    [Xall, ~] = stackSeqC(seqC);
    muGlob  = mean(Xall, 2, 'omitnan');
    varGlob = var(Xall, 0, 2, 'omitnan');
    varGlob = max(varGlob, opt.varFloor);

    if opt.verbose
        fprintf('[HMM] Fitting final model on all data at S=%d...\n', Sfinal);
    end

    [modelFinal, llFinal, histFinal] = hmm_fit_restarts_gauss_diag( ...
        seqC, Sfinal, muGlob, varGlob, opt);

    out.finalModel   = modelFinal;
    out.finalTrainLL = llFinal;
    out.finalHist    = histFinal;
end
end

% ========================================================================
%                               Helpers
% ========================================================================

function opt = parseNameValue(opt, varargin)
if mod(numel(varargin), 2) ~= 0
    error('Options must be name-value pairs.');
end
for i = 1:2:numel(varargin)
    name = varargin{i};
    val  = varargin{i+1};
    if ~isfield(opt, name)
        error('Unknown option: %s', name);
    end
    opt.(name) = val;
end
end

function seqC = coerceSeqC_KxT(seqC)
% Force each trial to K×T consistently.
% Strategy:
%   - If all trials share a common smaller dimension (K) and larger (T), we
%     infer K as the mode of min(size) across trials, then enforce K×T.
%   - If inference is ambiguous (rare), we fall back to heuristic: if rows>cols,
%     treat as T×K and transpose, else keep.
N = numel(seqC);
if N == 0
    error('seqC is empty.');
end

sz = zeros(N,2);
for i = 1:N
    X = seqC{i};
    if ~isnumeric(X) || isempty(X) || ndims(X) ~= 2
        error('seqC{%d} must be a non-empty 2D numeric matrix.', i);
    end
    sz(i,:) = size(X);
end

minDim = min(sz, [], 2);
Kmode  = mode(minDim);

for i = 1:N
    X = seqC{i};
    [r,c] = size(X);

    if r == Kmode
        % already K×T
        seqC{i} = X;
    elseif c == Kmode
        % transpose to K×T
        seqC{i} = X';
    else
        % ambiguous; fallback heuristic
        if r > c
            seqC{i} = X'; % likely T×K
        else
            seqC{i} = X;  % keep
        end
    end
end

% final sanity: enforce same K across all
K = size(seqC{1}, 1);
for i = 2:N
    if size(seqC{i}, 1) ~= K
        error('Inconsistent K after coercion: seqC{1} has K=%d but seqC{%d} has K=%d.', ...
            K, i, size(seqC{i},1));
    end
end
end

function foldIdx = makeKfoldIdx(N, Kfold, seed)
rng(seed);
perm = randperm(N);
foldIdx = zeros(N,1);
for i = 1:N
    foldIdx(perm(i)) = mod(i-1, Kfold) + 1;
end
end

function [Xall, lens] = stackSeqC(seqC)
% Stack K×T trials into K×sumT, return per-trial lengths.
N = numel(seqC);
K = size(seqC{1}, 1);
lens = zeros(N,1);
totT = 0;
for i = 1:N
    Xi = seqC{i};
    if size(Xi,1) ~= K
        error('stackSeqC: inconsistent feature dimension at trial %d.', i);
    end
    lens(i) = size(Xi,2);
    totT = totT + lens(i);
end
Xall = zeros(K, totT);
idx = 1;
for i = 1:N
    T = lens(i);
    Xall(:, idx:(idx+T-1)) = seqC{i};
    idx = idx + T;
end
end

% ========================================================================
%                         HMM fit + CV scoring
% ========================================================================

function [modelBest, bestLL, histBest] = hmm_fit_restarts_gauss_diag(seqC, S, muGlob, varGlob, opt)
bestLL   = -Inf;
modelBest = [];
histBest = [];

for r = 1:opt.nRestarts
    % deterministic restart stream (depends on S and restart index)
    rng(opt.rng + 1000*S + r);

    model0 = hmm_init_gauss_diag(S, muGlob, varGlob, opt);
    [model, LLhist] = hmm_em_gauss_diag(seqC, model0, opt);

    LLend = LLhist(find(~isnan(LLhist), 1, 'last'));
    if isempty(LLend), LLend = -Inf; end

    if LLend > bestLL
        bestLL   = LLend;
        modelBest = model;
        histBest.LL = LLhist;
        histBest.restart = r;
    end
end
end

function model = hmm_init_gauss_diag(S, muGlob, varGlob, opt)
K = numel(muGlob);

% init pi: near-uniform
pi0 = rand(S,1) + opt.piPrior;
pi0 = pi0 / sum(pi0);

% init A: sticky + random + pseudocount
A0 = rand(S,S) + opt.transPrior;
if opt.sticky > 0
    A0 = A0 + opt.sticky * eye(S);
end
A0 = A0 ./ sum(A0, 2);

% init emissions: spread means around global mean
mu0  = repmat(muGlob(:)', S, 1) + 0.1*randn(S,K).*sqrt(varGlob(:)');
sig20 = repmat(varGlob(:)', S, 1);

% enforce var floor
sig20 = max(sig20, opt.varFloor);

model = struct();
model.S    = S;
model.K    = K;
model.pi   = pi0;
model.A    = A0;
model.mu   = mu0;    % S×K
model.sig2 = sig20;  % S×K (diag variances)
end

function [model, LL_hist] = hmm_em_gauss_diag(seqC, model0, opt)
model = model0;

maxIter = opt.maxIter;
LL_hist = nan(maxIter,1);

prevLL = -Inf;

for it = 1:maxIter
    % E-step: expected sufficient stats across sequences
    [ll, stats] = estep_dataset_gauss_diag(seqC, model);

    LL_hist(it) = ll;

    % convergence check (relative improvement)
    if it > 1
        denom = max(1, abs(prevLL));
        relImpro = (ll - prevLL) / denom;
        if relImpro < opt.tolLL
            break;
        end
        % monotonic safeguard (rare with floors/prior, but keep)
        if ll < prevLL - 1e-9
            % rollback one step by stopping (model is still last updated from M-step,
            % so we just stop here; if you want strict rollback, store modelPrev)
            LL_hist(it) = prevLL;
            break;
        end
    end

    % M-step
    model = mstep_gauss_diag(stats, opt);

    prevLL = ll;
end

% truncate for cleanliness
last = find(~isnan(LL_hist), 1, 'last');
LL_hist = LL_hist(1:last);
end

function model = mstep_gauss_diag(stats, opt)
% stats fields:
%   gamma1Sum : S×1
%   xiSum     : S×S
%   gammaSum  : S×1  (total expected occupancy)
%   xSum      : S×K
%   x2Sum     : S×K
S = size(stats.xiSum,1);
K = size(stats.xSum,2);

% pi with pseudocount
pi = stats.gamma1Sum + opt.piPrior;
pi = pi / sum(pi);

% A with pseudocount + optional sticky encouragement handled in init only
A  = stats.xiSum + opt.transPrior;
A  = A ./ sum(A, 2);

% emissions
gamma = max(stats.gammaSum, eps); % S×1
mu  = stats.xSum ./ gamma;        % S×K (broadcast)
Ex2 = stats.x2Sum ./ gamma;       % S×K
sig2 = Ex2 - mu.^2;
sig2 = max(sig2, opt.varFloor);

model = struct();
model.S = S;
model.K = K;
model.pi = pi;
model.A  = A;
model.mu = mu;
model.sig2 = sig2;
end

function [ll, stats] = estep_dataset_gauss_diag(seqC, model)
S = model.S;
K = model.K;

gamma1Sum = zeros(S,1);
xiSum     = zeros(S,S);
gammaSum  = zeros(S,1);
xSum      = zeros(S,K);
x2Sum     = zeros(S,K);

ll = 0;

for n = 1:numel(seqC)
    X = seqC{n}; % K×T
    [lln, gamma, xi, ~] = estep_gauss_diag(X, model);

    ll = ll + lln;

    gamma1Sum = gamma1Sum + gamma(:,1);
    xiSum     = xiSum + xi;

    gsum = sum(gamma, 2); % S×1
    gammaSum = gammaSum + gsum;

    % accumulate emission stats
    % gamma: S×T, X: K×T
    % xSum(s,k) += sum_t gamma(s,t) * X(k,t)
    xSum  = xSum  + gamma * X';         % (S×T)*(T×K)=S×K
    x2Sum = x2Sum + gamma * (X'.^2);    % S×K

end

stats = struct();
stats.gamma1Sum = gamma1Sum;
stats.xiSum     = xiSum;
stats.gammaSum  = gammaSum;
stats.xSum      = xSum;
stats.x2Sum     = x2Sum;
end

function [ll, gamma, xiSum, logalpha] = estep_gauss_diag(X, model)
% E-step for a single K×T sequence under diagonal Gaussian HMM.
% Outputs:
%   ll      : log-likelihood of X
%   gamma   : S×T posterior p(z_t=s | X)
%   xiSum   : S×S expected transitions sum_{t} p(z_t=i, z_{t+1}=j | X)
%   logalpha: S×T forward messages (log)

S = model.S;

logpi = log(model.pi + eps);
logA  = log(model.A  + eps);

logB = log_emission_gauss_diag(X, model); % S×T

T = size(X,2);
logalpha = -Inf(S,T);
c = zeros(1,T); % scaling (log)

% forward init
logalpha(:,1) = logpi + logB(:,1);
c(1) = logsumexp(logalpha(:,1), 1);
logalpha(:,1) = logalpha(:,1) - c(1);

% forward recursion
for t = 2:T
    % tmp(j,i) = logA(i->j) + logalpha(i,t-1)
    tmp = logA' + logalpha(:,t-1); % S×S (j,i)
    logalpha(:,t) = logB(:,t) + logsumexp(tmp, 2);
    c(t) = logsumexp(logalpha(:,t), 1);
    logalpha(:,t) = logalpha(:,t) - c(t);
end

ll = sum(c);

% backward
logbeta = -Inf(S,T);
logbeta(:,T) = 0; % in scaled space

for t = T-1:-1:1
    % logbeta(i,t) = logsumexp_j [ logA(i->j) + logB(j,t+1) + logbeta(j,t+1) ]
    tmp = logA + (logB(:,t+1) + logbeta(:,t+1))'; % S×S (i,j)
    logbeta(:,t) = logsumexp(tmp, 2);
    logbeta(:,t) = logbeta(:,t) - c(t+1);
end

% gamma
loggamma = logalpha + logbeta;
loggamma = loggamma - logsumexp(loggamma, 1); % normalize per t
gamma = exp(loggamma);

% xi sum
xiSum = zeros(S,S);
for t = 1:T-1
    % log xi(i,j,t) ∝ logalpha(i,t) + logA(i,j) + logB(j,t+1) + logbeta(j,t+1)
    logxi = logalpha(:,t) + logA + (logB(:,t+1) + logbeta(:,t+1))';
    logxi = logxi - logsumexp(logxi(:), 1);
    xiSum = xiSum + exp(logxi);
end
end

function logB = log_emission_gauss_diag(X, model)
% Compute log p(x_t | z_t=s) for diagonal Gaussian.
% X: K×T
% mu: S×K, sig2: S×K
mu   = model.mu;     % S×K
sig2 = model.sig2;   % S×K
sig2 = max(sig2, 1e-12);

K = size(X,1);
T = size(X,2);
S = model.S;

logB = zeros(S,T);

% log N(x | mu, diag(sig2)) =
%  -0.5 * [ sum_k log(2*pi*sig2) + sum_k (x-mu)^2/sig2 ]
const = -0.5 * sum(log(2*pi*sig2), 2); % S×1

for s = 1:S
    % (K×T) - (K×1) -> K×T
    Xm = X - mu(s,:)';
    quad = -0.5 * sum((Xm.^2) ./ sig2(s,:)', 1); % 1×T
    logB(s,:) = const(s) + quad;
end
end

function [ll, nObs] = hmm_loglik_dataset(seqC, model)
ll = 0;
nObs = 0;
for n = 1:numel(seqC)
    X = seqC{n};
    ll = ll + hmm_loglik_single(X, model);
    nObs = nObs + size(X,2);
end
end

function ll = hmm_loglik_single(X, model)
% log-likelihood via scaled forward pass only
S = model.S;
logpi = log(model.pi + eps);
logA  = log(model.A  + eps);

logB = log_emission_gauss_diag(X, model); % S×T
T = size(X,2);

logalpha = -Inf(S,T);
c = zeros(1,T);

logalpha(:,1) = logpi + logB(:,1);
c(1) = logsumexp(logalpha(:,1), 1);
logalpha(:,1) = logalpha(:,1) - c(1);

for t = 2:T
    tmp = logA' + logalpha(:,t-1); % S×S (j,i)
    logalpha(:,t) = logB(:,t) + logsumexp(tmp, 2);
    c(t) = logsumexp(logalpha(:,t), 1);
    logalpha(:,t) = logalpha(:,t) - c(t);
end

ll = sum(c);
end

% ========================================================================
%                               Math utils
% ========================================================================

function y = logsumexp(A, dim)
% logsumexp with robust handling of all -Inf slices.
if nargin < 2, dim = 1; end
amax = max(A, [], dim);

% If amax is -Inf, the result should be -Inf (since all entries are -Inf)
isNegInf = ~isfinite(amax); % covers -Inf and NaN; NaN shouldn't happen, but guard

% subtract amax safely
Ashift = bsxfun(@minus, A, amax);
Ashift(~isfinite(Ashift)) = -Inf;

s = sum(exp(Ashift), dim);

y = amax + log(s + eps);

% force -Inf where appropriate
y(isNegInf) = -Inf;
end
