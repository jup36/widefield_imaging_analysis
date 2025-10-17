function [lambdaBest, lossPerLambda] = pick_lambda_ridge(Xz, y, lambdaGrid, K)
% Standardizes X, does K-fold CV ridge, returns best lambda by MSE.
% X: (n x p), y: (n x 1)

if nargin < 4, K = 5; end
cvp = cvpartition(numel(y), 'KFold', K);

% standardize once (like Mussall et al.)
% muX = mean(X,1,'omitnan');
% sdX = std (X,0,1,'omitnan'); 
% sdX(~isfinite(sdX) | sdX < 1e-3) = 1e-3;
% Xz = (X - muX) ./ sdX;

lossPerLambda = zeros(numel(lambdaGrid),1);

for li = 1:numel(lambdaGrid)
    lam = lambdaGrid(li);
    mse_k = zeros(K,1);

    for k = 1:K
        tr = training(cvp,k); te = test(cvp,k);
        Xtr = Xz(tr,:);  ytr = y(tr);
        Xte = Xz(te,:);  yte = y(te);

        % closed-form ridge (bias handled by pre-centering)
        B = (Xtr.'*Xtr + lam*eye(size(Xtr,2))) \ (Xtr.'*ytr);
        yhat = Xte*B;

        mse_k(k) = mean((yte - yhat).^2, 'omitnan');
    end
    lossPerLambda(li) = mean(mse_k);
end

[~, idx] = min(lossPerLambda);
lambdaBest = lambdaGrid(idx);
end
