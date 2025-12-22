function [modelBest, histBest] = hmm_em_gauss_diag_restart(seqC, S, opt)
%HMM_EM_GAUSS_DIAG_RESTART  Minimal Gaussian HMM (diag cov) with EM + restarts.
%
% DROP-IN USAGE
%   opt = struct('nRestarts',20,'maxIters',300,'tolRel',1e-6,'varFloor',1e-3, ...
%                'init','kmeans','sticky',0.85,'dataLayout','KxT','verbose',true);
%   [model, hist] = hmm_em_gauss_diag_restart(seqC, 4, opt);
%
% INPUT
%   seqC : 1xN cell, each sequence is either [K x T] or [T x K]
%   S    : # states
%   opt  : struct with optional fields
%          .maxIters    (default 200)
%          .tolRel      (default 1e-6)
%          .varFloor    (default 1e-3)
%          .nRestarts   (default 10)
%          .init        ('random'|'kmeans', default 'random')
%          .sticky      (0..1, default 0)  % A init with large diagonal if >0
%          .dataLayout  ('KxT'|'TxK', default 'KxT')
%          .verbose     (default true)
%
% OUTPUT
%   modelBest : struct with fields
%       .pi   : [S x 1]
%       .A    : [S x S] row-stochastic
%       .mu   : [K x S]
%       .var  : [K x S]  (diagonal variances per state)
%   histBest : struct
%       .LL        : [maxIters x 1] LL trace (NaN padded)
%       .restart   : best restart index
%
% NOTES
%   - E-step uses scaled forward-backward (stable) and returns monotone LL
%     for most practical cases. If LL drops, we rollback and stop that restart.
%   - Emissions are independent across K dims with diagonal variance.

% ---------------- defaults ----------------
if nargin < 3, opt = struct(); end
opt = setDefaultOpt(opt, 'maxIters',    200);
opt = setDefaultOpt(opt, 'tolRel',      1e-6);
opt = setDefaultOpt(opt, 'varFloor',    1e-3);
opt = setDefaultOpt(opt, 'nRestarts',   5);
opt = setDefaultOpt(opt, 'init',        'random');
opt = setDefaultOpt(opt, 'sticky',      0);
opt = setDefaultOpt(opt, 'dataLayout',  'KxT');
opt = setDefaultOpt(opt, 'verbose',     true);

% ---------------- global stats ----------------
Xall = catSeq(seqC, opt.dataLayout); % [K x Tall]
muGlob  = mean(Xall, 2, 'omitnan');
varGlob = var(Xall, 0, 2, 'omitnan');
varGlob(~isfinite(varGlob) | varGlob < opt.varFloor) = opt.varFloor;

bestLL   = -inf;
modelBest = [];
histBest  = [];

for r = 1:opt.nRestarts

    [pi0, A0, mu0, var0] = initParams(seqC, S, opt, muGlob, varGlob);
    model = struct('pi',pi0,'A',A0,'mu',mu0,'var',var0);

    LL_hist = nan(opt.maxIters,1);

    prevModel = model;
    prevLL    = -inf;

    for it = 1:opt.maxIters
        % ---------- E-step ----------
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

        % ---------- monotonicity safeguard ----------
        if it > 1 && LL < prevLL - 1e-8
            if opt.verbose
                fprintf('Restart %d | EM %3d | LL decreased (Δ=%.6g). Rollback & stop.\n', ...
                    r, it, LL - prevLL);
            end
            model = prevModel;
            LL_hist(it) = prevLL;
            break;
        end

        % ---------- convergence ----------
        if it > 1
            relImpro = (LL - prevLL) / max(1, abs(prevLL));
            if relImpro >= 0 && relImpro < opt.tolRel
                break;
            end
        end

        % ---------- M-step ----------
        prevModel = model;
        prevLL = LL;

        model = mstep_gauss_diag(stats, opt);

        % ---------- sanity ----------
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

if opt.verbose
    fprintf('Best restart: %d | best LL = %.3f\n', histBest.restart, bestLL);
end

end % ===== end main =====


%% ========================================================================
%                            Helper functions
% ========================================================================

function opt = setDefaultOpt(opt, field, val)
if ~isfield(opt, field) || isempty(opt.(field))
    opt.(field) = val;
end
end

function X = getSeq(Xin, layout)
if strcmpi(layout,'KxT')
    X = Xin;
elseif strcmpi(layout,'TxK')
    X = Xin';
else
    error('opt.dataLayout must be ''KxT'' or ''TxK''.');
end
end

function Xall = catSeq(seqC, layout)
% returns [K x Tall]
Xall = [];
for i = 1:numel(seqC)
    Xi = getSeq(seqC{i}, layout);
    Xall = [Xall, Xi]; %#ok<AGROW>
end
end

function Y = normalizeDim(X, dim)
Y = X ./ (sum(X, dim) + eps);
end

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

function stats = estep_gauss_diag(seqC, model, opt)
% Accumulate expected sufficient stats over sequences.

S = size(model.A,1);
K = size(model.mu,1);

sumGamma1 = zeros(S,1);
sumXi     = zeros(S,S);
sumGamma  = zeros(S,1);
sumX      = zeros(K,S);
sumXX     = zeros(K,S);

LL_total  = 0;

for n = 1:numel(seqC)
    X = getSeq(seqC{n}, opt.dataLayout); % [K x T]
    T = size(X,2);

    logB = log_emission_gauss_diag(X, model.mu, model.var); % [S x T]

    [~, ~, gamma, xi, LL] = fb_scaled(logB, model.pi, model.A);

    LL_total = LL_total + LL;

    sumGamma1 = sumGamma1 + gamma(:,1);
    sumXi     = sumXi     + sum(xi, 3);
    gsum      = sum(gamma,2);
    sumGamma  = sumGamma  + gsum;

    sumX  = sumX  + X * gamma';      % KxS
    sumXX = sumXX + (X.^2) * gamma'; % KxS
end

stats = struct();
stats.LL = LL_total;
stats.sumGamma1 = sumGamma1;
stats.sumXi = sumXi;
stats.sumGamma = sumGamma;
stats.sumX = sumX;
stats.sumXX = sumXX;
end

function model = mstep_gauss_diag(stats, opt)
S = numel(stats.sumGamma1);

piNew = stats.sumGamma1 / (sum(stats.sumGamma1) + eps);

Anew = stats.sumXi ./ (sum(stats.sumXi, 2) + eps); % row-normalize

muNew = stats.sumX ./ (stats.sumGamma' + eps);

Ex2   = stats.sumXX ./ (stats.sumGamma' + eps);
varNew = Ex2 - muNew.^2;
varNew(~isfinite(varNew) | varNew < opt.varFloor) = opt.varFloor;

model = struct('pi',piNew,'A',Anew,'mu',muNew,'var',varNew);
end

function logB = log_emission_gauss_diag(X, mu, varr)
% log p(x_t | z=s) with diagonal covariance (independent dims).
% X    : [K x T]
% mu   : [K x S]
% varr : [K x S]
% logB : [S x T]

[K, T] = size(X);
S = size(mu,2);

varr(~isfinite(varr) | varr <= 0) = eps;

logB = zeros(S, T);

% Compute per state: sum_k -0.5*(log(2πσ^2) + (x-μ)^2/σ^2)
const = -0.5 * log(2*pi);

for s = 1:S
    mu_s  = mu(:,s);
    var_s = varr(:,s);

    % broadcast across time
    d = X - mu_s;
    term = const - 0.5*log(var_s) - 0.5*(d.^2)./var_s;  % [K x T]
    logB(s,:) = sum(term, 1);                            % [1 x T]
end
end

function [alpha, beta, gamma, xi, LL] = fb_scaled(logB, pi0, A)
% Forward-backward with scaling (numerically stable).
%
% logB : [S x T] log emission likelihoods
% pi0  : [S x 1]
% A    : [S x S] row-stochastic
%
% Outputs:
%   alpha, beta : [S x T]
%   gamma       : [S x T]
%   xi          : [S x S x (T-1)]
%   LL          : scalar log-likelihood (scaled-domain; suitable for EM)

S = size(logB,1);
T = size(logB,2);

% Convert to stabilized emission probs
mx = max(logB, [], 1);
B  = exp(logB - mx);  % [S x T], each column scaled by exp(-max)

alpha = zeros(S,T);
c = zeros(1,T);

alpha(:,1) = pi0(:) .* B(:,1);
c(1) = sum(alpha(:,1)) + eps;
alpha(:,1) = alpha(:,1) / c(1);

for t = 2:T
    alpha(:,t) = (A' * alpha(:,t-1)) .* B(:,t);
    c(t) = sum(alpha(:,t)) + eps;
    alpha(:,t) = alpha(:,t) / c(t);
end

beta = zeros(S,T);
beta(:,T) = 1 / c(T);
for t = T-1:-1:1
    beta(:,t) = A * (B(:,t+1) .* beta(:,t+1));
    beta(:,t) = beta(:,t) / (c(t) + eps);
end

gamma = alpha .* beta;
gamma = gamma ./ (sum(gamma,1) + eps);

xi = zeros(S,S,T-1);
for t = 1:T-1
    % unnormalized xi(i,j,t) ∝ alpha(i,t)*A(i,j)*B(j,t+1)*beta(j,t+1)
    tmp = (alpha(:,t) * ( (B(:,t+1) .* beta(:,t+1))' )); % SxS
    tmp = tmp .* A;
    xi(:,:,t) = tmp / (sum(tmp,'all') + eps);
end

% LL: sum log c(t) plus the removed max(logB) per time bin
LL = sum(log(c + eps)) + sum(mx);  % add back mx to approximate true log-likelihood
end
