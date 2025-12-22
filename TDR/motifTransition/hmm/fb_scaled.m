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
