function [gamma, xi, ll] = fbScaled(logB, pi, A)
% Forward-backward with scaling for numerical stability.
% logB: S x T  (log emission probs)
% pi : S x 1
% A  : S x S
% gamma: S x T
% xi: S x S x (T-1)
% ll: scalar log-likelihood

[S,T] = size(logB);

B = exp(logB - max(logB,[],1));   % stabilize exp; columnwise
% (rescale columns later via alpha scaling anyway)

alpha = zeros(S,T);
c = zeros(1,T);

alpha(:,1) = pi .* B(:,1);
c(1) = sum(alpha(:,1)) + eps;
alpha(:,1) = alpha(:,1) / c(1);

for t = 2:T
    alpha(:,t) = (A.' * alpha(:,t-1)) .* B(:,t);
    c(t) = sum(alpha(:,t)) + eps;
    alpha(:,t) = alpha(:,t) / c(t);
end

beta = zeros(S,T);
beta(:,T) = ones(S,1) / c(T);

for t = T-1:-1:1
    beta(:,t) = A * (B(:,t+1) .* beta(:,t+1));
    beta(:,t) = beta(:,t) / (c(t) + eps);
end

gamma = alpha .* beta;
gamma = gamma ./ (sum(gamma,1) + eps);

xi = zeros(S,S,max(T-1,1));
if T > 1
    for t = 1:T-1
        tmp = (alpha(:,t) * (B(:,t+1) .* beta(:,t+1)).'); % S x S (outer)
        tmp = tmp .* A;
        xi(:,:,t) = tmp ./ (sum(tmp,'all') + eps);
    end
end

ll = sum(log(c + eps)); % log-likelihood
end