function [ll, gamma, xiSum, logalpha] = estep_gauss_diag(X, model)
S = model.S;

logpi = log(model.pi + eps);
logA  = log(model.A  + eps);

logB = log_emission_gauss_diag(X, model); % S×T
T = size(X,2);

logalpha = -Inf(S,T);
c = zeros(1,T);

logalpha(:,1) = logpi + logB(:,1);
c(1) = logsumexp(logalpha(:,1), 1);
logalpha(:,1) = logalpha(:,1) - c(1);

for t = 2:T
    tmp = logA' + logalpha(:,t-1); % tmp(j,i) = logA(i->j)+logalpha(i)
    logalpha(:,t) = logB(:,t) + logsumexp(tmp, 2);
    c(t) = logsumexp(logalpha(:,t), 1);
    logalpha(:,t) = logalpha(:,t) - c(t);
end

ll = sum(c);

logbeta = -Inf(S,T);
logbeta(:,T) = 0;

for t = T-1:-1:1
    tmp = logA + (logB(:,t+1) + logbeta(:,t+1))';
    logbeta(:,t) = logsumexp(tmp, 2);
    logbeta(:,t) = logbeta(:,t) - c(t+1);
end

loggamma = logalpha + logbeta;
loggamma = loggamma - logsumexp(loggamma, 1);
gamma = exp(loggamma);

xiSum = zeros(S,S);
for t = 1:T-1
    logxi = logalpha(:,t) + logA + (logB(:,t+1) + logbeta(:,t+1))';
    logxi = logxi - logsumexp(logxi(:), 1);
    xiSum = xiSum + exp(logxi);
end
end



function y = logsumexp(A, dim)
if nargin < 2, dim = 1; end
amax = max(A, [], dim);
isNegInf = ~isfinite(amax);

Ashift = bsxfun(@minus, A, amax);
Ashift(~isfinite(Ashift)) = -Inf;

s = sum(exp(Ashift), dim);
y = amax + log(s + eps);

y(isNegInf) = -Inf;
end
