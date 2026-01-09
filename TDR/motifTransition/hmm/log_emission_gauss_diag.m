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