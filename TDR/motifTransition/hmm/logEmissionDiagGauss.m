function logB = logEmissionDiagGauss(X, mu, v)
% X: D x T, mu: D x S, v: D x S
% returns logB: S x T
[D,T] = size(X);
S = size(mu,2);

% compute log N(x | mu_s, diag(v_s)) for all s,t
logB = zeros(S,T);
const = -0.5 * D * log(2*pi);

for s = 1:S
    vs = v(:,s);
    invv = 1 ./ vs;
    logdet = -0.5 * sum(log(vs));
    % quadratic term per time
    Xm = X - mu(:,s);
    quad = -0.5 * sum((Xm.^2) .* invv, 1);   % 1 x T
    logB(s,:) = const + logdet + quad;
end
end