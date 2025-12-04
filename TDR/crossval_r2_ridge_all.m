function [cvR2_full, cvR2_global] = crossval_r2_ridge_all(X, Y, lambdas, cv)
% Use this 'cv_global_r2_ridge.m' instead for 'trial-aware' cross validation
% CROSSVAL_R2_RIDGE_ALL  Cross-validated R^2 for a ridge GLM (per motif and global).
%
% SYNTAX
%   [cvR2_full, cvR2_global] = crossval_r2_ridge_all(X, Y, lambdas, cv)
%
% INPUTS
%   X        : (M x P) design matrix, already standardized if desired.
%              M = #samples (e.g., trials×time-bins), P = #predictors.
%   Y        : (M x K) response matrix (K motifs/units). Each column is a target.
%   lambdas  : (1 x K) ridge penalties for each motif (can be scalar -> applied to all).
%   cv       : cvpartition object over M samples; F = cv.NumTestSets folds.
%
% OUTPUTS
%   cvR2_full   : (K x 1) CV R^2 for each motif (variance explained on held-out data).
%   cvR2_global : scalar, single “global” CV R^2 pooling SSE/SST across all motifs.
%
% DEFINITIONS / NOTES
%   - For each fold f = 1..F:
%       * Fit ridge per motif k on training data:  β_k = (X_tr'X_tr + λ_k I) \ (X_tr' y_tr,k).
%       * Predict on held-out data: ŷ_te,k = X_te β_k.
%       * Compute SSE_te,k = sum((y_te,k − ŷ_te,k).^2).
%       * Compute SST_te,k = sum((y_te,k − μ_tr,k).^2), where μ_tr,k = mean(y_tr,k).
%         (Using the training mean within each fold to define chance/baseline.)
%   - cvR2_full(k) = 1 − (Σ_f SSE_te,k) / (Σ_f SST_te,k).
%   - cvR2_global  = 1 − (Σ_f,k SSE_te,k) / (Σ_f,k SST_te,k). This is a variance-weighted
%     aggregate across motifs (not a simple average of per-motif R^2).
%   - X and Y rows must be aligned. Any standardization of X should be done before calling.
%   - Y can be z-scored per motif beforehand if you want comparable scales; this does not
%     change R^2 but can stabilize λ selection elsewhere.
%
% Junchol Park / Buschman Lab — 2025

[M, P] = size(X);
K = size(Y, 2);
F = cv.NumTestSets;

% Allow scalar lambda → broadcast to all motifs
if isscalar(lambdas), lambdas = repmat(lambdas, 1, K); end
assert(isvector(lambdas) && numel(lambdas)==K, 'lambdas must be 1xK or scalar.');

I = speye(P);

% Accumulators across folds
SSE_k = zeros(K,1);
SST_k = zeros(K,1);
SSE_global = 0;
SST_global = 0;

for f = 1:F
    tr = training(cv, f);
    te = test(cv, f);

    Xtr = X(tr,:);  Xte = X(te,:);
    Ytr = Y(tr,:);  Yte = Y(te,:);

    % Baseline per motif uses TRAINING mean (fold-safe)
    mu_tr = mean(Ytr, 1);

    for k = 1:K
        lam  = lambdas(k);
        beta = (Xtr.'*Xtr + lam*I) \ (Xtr.'*Ytr(:,k));   % ridge closed form
        yhat = Xte * beta;

        err  = Yte(:,k) - yhat;
        SSE  = sum(err.^2);
        SST  = sum( (Yte(:,k) - mu_tr(k)).^2 );

        SSE_k(k) = SSE_k(k) + SSE;
        SST_k(k) = SST_k(k) + SST;

        SSE_global = SSE_global + SSE;
        SST_global = SST_global + SST;
    end
end

% Per-motif CV R^2
cvR2_full = 1 - SSE_k ./ max(SST_k, eps);

% Global CV R^2 pooled across motifs
cvR2_global = 1 - SSE_global / max(SST_global, eps);
end
