function dec = hmm_decode_trial(X, model)
% Return time-by-time emission log-likelihoods and state posteriors for one trial.
% X: K×T

% Emission log-likelihood (S×T)
logB = log_emission_gauss_diag(X, model);

% Run forward-backward to get gamma and xiSum
[~, gamma, xiSum, logalpha] = estep_gauss_diag(X, model);

dec = struct();
dec.logB     = logB;      % S×T (log P(x_t | z_t=s))
dec.B        = exp(logB); % S×T (often extremely small; log is safer)
dec.gamma    = gamma;     % S×T (P(z_t=s | x_1:T))
dec.xiSum    = xiSum;     % S×S (expected transitions in this trial)
dec.logalpha = logalpha;  % S×T (optional, filtered info up to t)
end
