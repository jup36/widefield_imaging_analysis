function [R2_per, R2_combined, details] = cv_global_r2_ridge(X, Y, rowTrialIdx, varargin)
%CV_GLOBAL_R2_RIDGE  Trial-wise CV R^2 for multi-target ridge GLM (NaN-robust).
%
% [R2_per, R2_combined, details] = cv_global_r2_ridge(X, Y, rowTrialIdx, ...
%   'outerFolds',5, 'lambda',[], 'lambdaGrid',logspace(-3,3,15), ...
%   'lambdaMode','perTarget', 'innerFolds',5, 'standardize',true, 'rng',1)
%
% NaN handling:
%   • X: standardize on train; set NaNs -> 0 in z-space (column mean).
%   • Y: for each target, fit using only rows with finite Y; evaluate on
%        test rows with finite Y.

% ----- parse -----
p = inputParser;
p.addParameter('outerFolds', 5, @(x)isscalar(x)&&x>=2);
p.addParameter('lambda', [], @(x) isempty(x) || isscalar(x) || isvector(x));
p.addParameter('lambdaGrid', logspace(-3,3,15), @(x)isvector(x)&&all(x>=0));
p.addParameter('lambdaMode','perTarget', @(s) any(strcmpi(s,{'perTarget','shared'})));
p.addParameter('innerFolds', 5, @(x)isscalar(x)&&x>=2);
p.addParameter('standardize', true, @(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('rng', 1, @(x)isscalar(x));
p.parse(varargin{:});
opt = p.Results;

% ----- sanity -----
X = double(X);
Y = double(Y);
[nRows, pX] = size(X);
[ny, q] = size(Y);
if ny ~= nRows, error('Y must have same #rows as X.'); end
if numel(rowTrialIdx) ~= nRows, error('rowTrialIdx must be nRows-by-1.'); end

% ----- outer folds by trials -----
uTrials = unique(rowTrialIdx(:));
nTrials = numel(uTrials);
rng(opt.rng);
perm = uTrials(randperm(nTrials));
cuts = round(linspace(0,nTrials,opt.outerFolds+1));
foldTrials = cell(opt.outerFolds,1);
for f = 1:opt.outerFolds
    foldTrials{f} = perm(cuts(f)+1 : cuts(f+1));
end

% ----- accumulators -----
yhat_all = nan(nRows, q);
testMask = false(nRows, opt.outerFolds);

SSE_tot  = zeros(1,q);
SST_tot  = zeros(1,q);
foldR2   = nan(opt.outerFolds, q);
lam_used = nan(opt.outerFolds, q);
nTestFinite = zeros(1,q);   % count #finite test samples per target

% store scalers (optional)
muX_cell = cell(opt.outerFolds,1);
sdX_cell = cell(opt.outerFolds,1);
muy_cell = cell(opt.outerFolds,1);     % [1 x q]

I = speye(pX);

for f = 1:opt.outerFolds
    teTr = foldTrials{f};
    te   = ismember(rowTrialIdx, teTr);
    tr   = ~te;

    Xtr = X(tr,:);    Xte = X(te,:);
    Ytr = Y(tr,:);    Yte = Y(te,:);

    % ---- standardize X on TRAIN only ----
    if opt.standardize
        muX = mean(Xtr,1,'omitnan');
        sdX = std(Xtr,0,1,'omitnan');
        sdX(~isfinite(sdX)|sdX<1e-6) = 1e-6;
        XtrZ = (Xtr - muX) ./ sdX;
        XteZ = (Xte - muX) ./ sdX;
    else
        muX = zeros(1,pX); sdX = ones(1,pX);
        XtrZ = Xtr; XteZ = Xte;
    end
    % Impute NaNs in X (train/test) to 0 in z-space
    XtrZ(~isfinite(XtrZ)) = 0;
    XteZ(~isfinite(XteZ)) = 0;

    % Center Y per target on TRAIN (baseline)
    muy = mean(Ytr,1,'omitnan');   % [1 x q]
    YtrC = Ytr - muy;              % OK to contain NaNs

    % choose lambda(s)
    if isempty(opt.lambda)
        switch lower(opt.lambdaMode)
            case 'pertarget'
                lam_f = pick_lambda_per_target_nanrobust(XtrZ, YtrC, rowTrialIdx(tr), opt.lambdaGrid(:)', opt.innerFolds, opt.rng+f);
            case 'shared'
                lam_shared = pick_lambda_shared_nanrobust(XtrZ, YtrC, rowTrialIdx(tr), opt.lambdaGrid(:)', opt.innerFolds, opt.rng+f);
                lam_f = lam_shared * ones(1,q);
        end
    else
        if isscalar(opt.lambda)
            lam_f = opt.lambda * ones(1,q);
        else
            if numel(opt.lambda) ~= q
                error('lambda vector must have length equal to #targets (columns of Y).');
            end
            lam_f = opt.lambda(:)';   % 1 x q
        end
    end

    % fit/predict per target with NaN-safe masks
    XtX_full = XtrZ' * XtrZ;   % can reuse when masks are dense, but we’ll refit with masks
    for t = 1:q
        lam = lam_f(t);

        mask_tr = isfinite(YtrC(:,t));
        if ~any(mask_tr)
            % no training data for this target in this fold
            yhat_te = nan(sum(te),1);
        else
            % Solve with masked rows
            Xtr_t = XtrZ(mask_tr,:);
            ytr_t = YtrC(mask_tr,t);
            beta  = (Xtr_t' * Xtr_t + lam * I) \ (Xtr_t' * ytr_t);
            yhat_te = XteZ * beta + muy(t);
        end

        % accumulate metrics on TEST using finite Y only
        y    = Yte(:,t);
        mask_te = isfinite(y);
        nTestFinite(t) = nTestFinite(t) + sum(mask_te);
        err  = y(mask_te) - yhat_te(mask_te);
        sse  = sum(err.^2);
        sst  = sum( (y(mask_te) - muy(t)).^2 );

        SSE_tot(t) = SSE_tot(t) + sse;
        SST_tot(t) = SST_tot(t) + sst;

        foldR2(f,t) = 1 - sse / max(sst, eps);

        % stash predictions
        yhat_all(te,t) = yhat_te;
        lam_used(f,t)  = lam;
    end

    testMask(te,f) = true;
    muX_cell{f} = muX; sdX_cell{f} = sdX; muy_cell{f} = muy;
end

% ----- final metrics -----
R2_per = 1 - SSE_tot ./ max(SST_tot, eps);                 % [1 x q]
R2_per(nTestFinite==0) = NaN;   % mark targets with no finite test data
R2_combined = 1 - sum(SSE_tot) / max(sum(SST_tot), eps);   % scalar

% ----- outputs -----
details = struct();
details.foldR2     = foldR2;
details.lambda     = lam_used;
details.yhat       = yhat_all;
details.testMask   = testMask;
details.muX        = muX_cell;
details.sdX        = sdX_cell;
details.muy        = muy_cell;
details.outerFolds = opt.outerFolds;
details.innerFolds = opt.innerFolds;
details.lambdaMode = opt.lambdaMode;

end % main


% ===== helpers (NaN-robust) =====

function lam_vec = pick_lambda_per_target_nanrobust(XtrZ, YtrC, trialIdx_tr, lamGrid, K, seed)
% Independently pick λ per target by inner CV (trial-wise), NaN-robust.
    rng(seed);
    uTr = unique(trialIdx_tr(:)); nTr = numel(uTr);
    perm = uTr(randperm(nTr));
    cuts = round(linspace(0,nTr,K+1));

    q = size(YtrC,2); p = size(XtrZ,2);
    I = speye(p);
    lam_vec = nan(1,q);

    % Pre-impute NaNs in X to 0 (already z-space)
    XtrZ(~isfinite(XtrZ)) = 0;

    for t = 1:q
        mses = zeros(numel(lamGrid), K);
        for k = 1:K
            teTr = perm(cuts(k)+1 : cuts(k+1));
            te   = ismember(trialIdx_tr, teTr); tr = ~te;

            Xtr_k = XtrZ(tr,:); Xte_k = XtrZ(te,:); % already NaN->0
            ytr   = YtrC(tr,t);
            yte   = YtrC(te,t) + mean(YtrC(tr,t),'omitnan'); % not used directly for mse baseline here

            % training mask: finite y only
            mtr = isfinite(ytr);
            Xtr_t = Xtr_k(mtr,:);  ytr_t = ytr(mtr);

            % test mask: finite y only (for mse)
            mte = isfinite(YtrC(te,t) + mean(YtrC(tr,t),'omitnan'));
            yte_raw = YtrC(te,t) + mean(YtrC(tr,t),'omitnan'); % back to original scale

            for i = 1:numel(lamGrid)
                lam = lamGrid(i);
                if any(mtr)
                    beta = (Xtr_t' * Xtr_t + lam*I) \ (Xtr_t' * ytr_t);
                    yhat = Xte_k * beta + mean(YtrC(tr,t),'omitnan');
                else
                    yhat = nan(sum(te),1);
                end
                % MSE on finite test rows
                err = yte_raw(mte) - yhat(mte);
                mses(i,k) = mean(err.^2, 'omitnan');
            end
        end
        [~,ix] = min(mean(mses,2,'omitnan'));
        lam_vec(t) = lamGrid(ix);
    end
end

function lam = pick_lambda_shared_nanrobust(XtrZ, YtrC, trialIdx_tr, lamGrid, K, seed)
% One λ minimizing sum of MSE across all targets (trial-wise), NaN-robust.
    rng(seed);
    uTr = unique(trialIdx_tr(:)); nTr = numel(uTr);
    perm = uTr(randperm(nTr));
    cuts = round(linspace(0,nTr,K+1));

    p = size(XtrZ,2);
    I = speye(p);
    mse_sum = zeros(numel(lamGrid), K);

    % Impute NaNs in X to 0 (z-space)
    XtrZ(~isfinite(XtrZ)) = 0;

    for k = 1:K
        teTr = perm(cuts(k)+1 : cuts(k+1));
        te   = ismember(trialIdx_tr, teTr); tr = ~te;

        Xtr_k = XtrZ(tr,:); Xte_k = XtrZ(te,:);
        Ytr_k = YtrC(tr,:); Yte_k = YtrC(te,:);

        mu_tr = mean(Ytr_k,1,'omitnan');  % baseline per target

        for i = 1:numel(lamGrid)
            lam_i = lamGrid(i);
            % Fit per target with its own finite mask, but same lambda
            Yhat_te = nan(sum(te), size(Ytr_k,2));
            for t = 1:size(Ytr_k,2)
                mtr = isfinite(Ytr_k(:,t));
                if any(mtr)
                    beta_t = (Xtr_k(mtr,:)'*Xtr_k(mtr,:) + lam_i*I) \ (Xtr_k(mtr,:)'*Ytr_k(mtr,t));
                    Yhat_te(:,t) = Xte_k * beta_t + mu_tr(t);
                end
            end
            % MSE summed across targets on finite test rows per target
            err2_sum = 0; ncount = 0;
            for t = 1:size(Ytr_k,2)
                yte_t = Yte_k(:,t) + mu_tr(t);
                mte   = isfinite(yte_t);
                if any(mte)
                    err2_sum = err2_sum + mean((yte_t(mte) - Yhat_te(mte,t)).^2, 'omitnan');
                    ncount   = ncount + 1;
                end
            end
            mse_sum(i,k) = err2_sum / max(ncount,1);
        end
    end
    [~,ix] = min(mean(mse_sum,2,'omitnan'));
    lam = lamGrid(ix);
end
