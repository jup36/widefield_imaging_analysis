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