function stats = estep_gauss_diag(seqC, model, opt)
% Accumulate expected sufficient stats over sequences.

S = size(model.A,1);
K = size(model.mu,1);

sumGamma1 = zeros(S,1);
sumXi     = zeros(S,S);
sumGamma  = zeros(S,1);
sumX      = zeros(K,S);
sumXX     = zeros(K,S);

LL_total  = 0;

for n = 1:numel(seqC)
    X = getSeq(seqC{n}, opt.dataLayout); % [K x T]
    T = size(X,2);

    logB = log_emission_gauss_diag(X, model.mu, model.var); % [S x T]

    [~, ~, gamma, xi, LL] = fb_scaled(logB, model.pi, model.A);

    LL_total = LL_total + LL;

    sumGamma1 = sumGamma1 + gamma(:,1);
    sumXi     = sumXi     + sum(xi, 3);
    gsum      = sum(gamma,2);
    sumGamma  = sumGamma  + gsum;

    sumX  = sumX  + X * gamma';      % KxS
    sumXX = sumXX + (X.^2) * gamma'; % KxS
end

stats = struct();
stats.LL = LL_total;
stats.sumGamma1 = sumGamma1;
stats.sumXi = sumXi;
stats.sumGamma = sumGamma;
stats.sumX = sumX;
stats.sumXX = sumXX;
end