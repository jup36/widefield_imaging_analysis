function out = hmm_modelSelect_withinSession_cv(seqC, varargin)
% hmm_modelSelect_withinSession_cv
% -------------------------------------------------------------------------
% Drop-in, within-session diagonal-Gaussian HMM model selection with K-fold CV.
%
% UPDATED (per request):
%   - Instead of saving only out.finalModel for the single "best" S,
%     this version saves *all* fitted full-data models as:
%         out.models{sIdx}  (cell array; each is an HMM struct)
%     aligned with out.S_list(sIdx).
%
% What this does (synopsis)
%   1) Takes trial-by-trial data in seqC, each trial is K×T (or T×K; will coerce).
%   2) For each candidate S:
%        - K-fold CV across trials:
%            * Fit on train trials (EM with restarts), score test trials.
%        - Store per-fold and pooled test log-likelihood per bin.
%   3) Fit a *full-data* model for EACH S (best over restarts) and store it in out.models.
%      This supports post-hoc elbow selection (choose S by curve shape, then pull model).
%
% Inputs
%   seqC : 1×N cell array, each cell is K×T or T×K numeric
%
% Name-Value options
%   'S_list'        : candidate state counts (default 2:10)
%   'Kfold'         : folds (default 5)
%   'maxIter'       : EM iterations (default 200)
%   'tolLL'         : relative LL tolerance (default 1e-6)
%   'nRestarts'     : random restarts (default 5)
%   'sticky'        : sticky init boost on diag(A) (default 0.0)
%   'varFloor'      : variance floor (default 1e-4)
%   'transPrior'    : Dirichlet pseudocount for A (default 1e-2)
%   'piPrior'       : Dirichlet pseudocount for pi (default 1e-2)
%   'rng'           : rng seed (default 0)
%   'verbose'       : true/false (default true)
%   'fitModelsOnAll': fit and store full-data model for each S (default true)
%
% Output (struct) out
%   .opt                 : options used
%   .S_list              : candidate state counts
%   .cv(sIdx)            : CV results for S_list(sIdx)
%      .S
%      .foldLL_perObs
%      .foldNobs
%      .meanFoldLL
%      .pooledLL
%   .bestS_by_pooledLL
%   .bestS_by_meanFoldLL
%   .models{sIdx}        : full-data fitted model for S_list(sIdx) (if enabled)
%   .modelTrainLL(sIdx)  : full-data train LL for that model (sum over all trials)
%   .modelHist{sIdx}     : EM LL history for that model (best restart)
%
% -------------------------------------------------------------------------

% --------------------------- parse options ------------------------------
opt = struct();
opt.S_list        = 2:10;
opt.Kfold         = 5;
opt.maxIter       = 200;
opt.tolLL         = 1e-6;
opt.nRestarts     = 5;
opt.sticky        = 0.0;
opt.varFloor      = 1e-4;
opt.transPrior    = 1e-2;
opt.piPrior       = 1e-2;
opt.rng           = 0;
opt.verbose       = true;
opt.fitModelsOnAll = true;

opt = parseNameValue(opt, varargin{:});

% Reproducibility for folds + restarts
rng(opt.rng);

% --------------------------- validate & coerce ---------------------------
seqC = coerceSeqC_KxT(seqC);

Ntr = numel(seqC);
K   = size(seqC{1}, 1);

% --------------------------- build folds --------------------------------
foldIdx = makeKfoldIdx(Ntr, opt.Kfold, opt.rng);

% --------------------------- init output --------------------------------
out = struct();
out.opt    = opt;
out.S_list = opt.S_list;
out.cv     = struct();

% NEW: store all full-data models
out.models      = cell(numel(opt.S_list), 1);
out.modelTrainLL = nan(numel(opt.S_list), 1);
out.modelHist   = cell(numel(opt.S_list), 1);

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
        [modelBest, llTrainBest, ~] = hmm_fit_restarts_gauss_diag( ...
            trainSeq, S, muGlob, varGlob, opt);

        % Evaluate held-out test LL (no parameter update)
        [llTest, nObsTest] = hmm_loglik_dataset(testSeq, modelBest);

        foldLL_perObs(kf) = llTest / max(1, nObsTest);
        foldNobs(kf)      = nObsTest;

        if opt.verbose
            fprintf('  Fold %d/%d: trainLL=%.3e, testLL/obs=%.6f (Nobs=%d)\n', ...
                kf, opt.Kfold, llTrainBest, foldLL_perObs(kf), nObsTest);
        end
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

    % ------------------- Fit full-data model for this S ------------------
    if opt.fitModelsOnAll
        [Xall, ~] = stackSeqC(seqC);
        muGlobAll  = mean(Xall, 2, 'omitnan');
        varGlobAll = var(Xall, 0, 2, 'omitnan');
        varGlobAll = max(varGlobAll, opt.varFloor);

        if opt.verbose
            fprintf('  [HMM] Fitting full-data model for S=%d...\n', S);
        end

        [modelAll, llAll, histAll] = hmm_fit_restarts_gauss_diag( ...
            seqC, S, muGlobAll, varGlobAll, opt);

        out.models{sIdx}       = modelAll;
        out.modelTrainLL(sIdx) = llAll;
        out.modelHist{sIdx}    = histAll;
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

% Convenience: index lookup for a chosen S
out.getModelIdx = @(S) find(out.S_list == S, 1, 'first');
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
%   - infer K as mode(min(size(trial))) across trials
%   - if a trial is T×K, transpose to K×T
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
        seqC{i} = X;
    elseif c == Kmode
        seqC{i} = X';
    else
        % ambiguous; fallback heuristic
        if r > c
            seqC{i} = X';
        else
            seqC{i} = X;
        end
    end
end

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
histBest = struct('LL', [], 'restart', []);

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

pi0 = rand(S,1) + opt.piPrior;
pi0 = pi0 / sum(pi0);

A0 = rand(S,S) + opt.transPrior;
if opt.sticky > 0
    A0 = A0 + opt.sticky * eye(S);
end
A0 = A0 ./ sum(A0, 2);

mu0   = repmat(muGlob(:)', S, 1) + 0.1*randn(S,K).*sqrt(varGlob(:)');
sig20 = repmat(varGlob(:)', S, 1);
sig20 = max(sig20, opt.varFloor);

model = struct();
model.S    = S;
model.K    = K;
model.pi   = pi0;
model.A    = A0;
model.mu   = mu0;    % S×K
model.sig2 = sig20;  % S×K
end

function [model, LL_hist] = hmm_em_gauss_diag(seqC, model0, opt)
model = model0;

LL_hist = nan(opt.maxIter,1);
prevLL  = -Inf;

for it = 1:opt.maxIter
    [ll, stats] = estep_dataset_gauss_diag(seqC, model);
    LL_hist(it) = ll;

    if it > 1
        denom = max(1, abs(prevLL));
        relImpro = (ll - prevLL) / denom;
        if relImpro < opt.tolLL
            break;
        end
        if ll < prevLL - 1e-9
            LL_hist(it) = prevLL;
            break;
        end
    end

    model = mstep_gauss_diag(stats, opt);
    prevLL = ll;
end

last = find(~isnan(LL_hist), 1, 'last');
LL_hist = LL_hist(1:last);
end

function model = mstep_gauss_diag(stats, opt)
S = size(stats.xiSum,1);
K = size(stats.xSum,2);

pi = stats.gamma1Sum + opt.piPrior;
pi = pi / sum(pi);

A  = stats.xiSum + opt.transPrior;
A  = A ./ sum(A, 2);

gamma = max(stats.gammaSum, eps); % S×1
mu  = stats.xSum ./ gamma;        % S×K
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
    X = seqC{n};
    [lln, gamma, xi, ~] = estep_gauss_diag(X, model);

    ll = ll + lln;

    gamma1Sum = gamma1Sum + gamma(:,1);
    xiSum     = xiSum + xi;

    gammaSum = gammaSum + sum(gamma, 2);

    xSum  = xSum  + gamma * X';
    x2Sum = x2Sum + gamma * (X'.^2);
end

stats = struct();
stats.gamma1Sum = gamma1Sum;
stats.xiSum     = xiSum;
stats.gammaSum  = gammaSum;
stats.xSum      = xSum;
stats.x2Sum     = x2Sum;
end

function [ll, gamma, xiSum, logalpha] = estep_gauss_diag(X, model)
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
    tmp = logA' + logalpha(:,t-1); % tmp(j,i) = logA(i->j)+logalpha(i)
    logalpha(:,t) = logB(:,t) + logsumexp(tmp, 2);
    c(t) = logsumexp(logalpha(:,t), 1);
    logalpha(:,t) = logalpha(:,t) - c(t);
end

ll = sum(c);

logbeta = -Inf(S,T);
logbeta(:,T) = 0;

for t = T-1:-1:1
    tmp = logA + (logB(:,t+1) + logbeta(:,t+1))';
    logbeta(:,t) = logsumexp(tmp, 2);
    logbeta(:,t) = logbeta(:,t) - c(t+1);
end

loggamma = logalpha + logbeta;
loggamma = loggamma - logsumexp(loggamma, 1);
gamma = exp(loggamma);

xiSum = zeros(S,S);
for t = 1:T-1
    logxi = logalpha(:,t) + logA + (logB(:,t+1) + logbeta(:,t+1))';
    logxi = logxi - logsumexp(logxi(:), 1);
    xiSum = xiSum + exp(logxi);
end
end

function logB = log_emission_gauss_diag(X, model)
mu   = model.mu;     % S×K
sig2 = max(model.sig2, 1e-12);

S = model.S;
T = size(X,2);

logB = zeros(S,T);

const = -0.5 * sum(log(2*pi*sig2), 2); % S×1
for s = 1:S
    Xm = X - mu(s,:)';
    quad = -0.5 * sum((Xm.^2) ./ sig2(s,:)', 1);
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
S = model.S;
logpi = log(model.pi + eps);
logA  = log(model.A  + eps);

logB = log_emission_gauss_diag(X, model);
T = size(X,2);

logalpha = -Inf(S,T);
c = zeros(1,T);

logalpha(:,1) = logpi + logB(:,1);
c(1) = logsumexp(logalpha(:,1), 1);
logalpha(:,1) = logalpha(:,1) - c(1);

for t = 2:T
    tmp = logA' + logalpha(:,t-1);
    logalpha(:,t) = logB(:,t) + logsumexp(tmp, 2);
    c(t) = logsumexp(logalpha(:,t), 1);
    logalpha(:,t) = logalpha(:,t) - c(t);
end

ll = sum(c);
end

function y = logsumexp(A, dim)
if nargin < 2, dim = 1; end
amax = max(A, [], dim);
isNegInf = ~isfinite(amax);

Ashift = bsxfun(@minus, A, amax);
Ashift(~isfinite(Ashift)) = -Inf;

s = sum(exp(Ashift), dim);
y = amax + log(s + eps);

y(isNegInf) = -Inf;
end
