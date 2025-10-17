function r2 = crossval_r2_ridge(X, y, lambda, cv)
% Returns mean R^2 across folds (held-out).
r2fold = zeros(cv.NumTestSets,1);
p = size(X,2);
I = speye(p);

for f = 1:cv.NumTestSets
    tr = training(cv,f);
    te = test(cv,f);

    Xtr = X(tr,:); ytr = y(tr);
    Xte = X(te,:); yte = y(te);

    % Ridge solution via Cholesky (stable & fast)
    A = Xtr' * Xtr + lambda * I;
    b = Xtr' * ytr;
    beta = A \ b;

    yhat = Xte * beta;
    ss_res = sum((yte - yhat).^2);
    ss_tot = sum((yte - mean(yte)).^2);
    r2fold(f) = 1 - ss_res / max(ss_tot, eps);
end
r2 = mean(r2fold);
end