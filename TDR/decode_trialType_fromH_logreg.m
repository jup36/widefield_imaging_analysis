function decRez = decode_trialType_fromH_logreg(HsY3d, Htime, trI, varargin)
%DECODE_TRIALTYPE_FROMH_LOGREG  Go vs NoGo decoding from motif activity H.
%
% decRez = decode_trialType_fromH_logreg(HsY3d, Htime, trI, ...
%             'timeWindow', [0 2], ...
%             'KFold', 5, ...
%             'Standardize', true, ...
%             'rng', 1);
%
% INPUTS
%   HsY3d : [N x K x T]  motif activity (e.g., HsZ.Y3 from stack_trials_H)
%           N = # trials, K = # motifs, T = # time bins
%
%   Htime : [1 x T] time vector for the 3rd dimension of HsY3d (e.g., Hs.winCtrs)
%
%   trI   : struct with logical trial-type masks
%           .goI   : [N x 1] true for Go trials
%           .nogoI : [N x 1] true for NoGo trials
%
% Name–Value pairs
%   'timeWindow' : [tStart tEnd] (default = full range of Htime)
%                  Time window (in same units as Htime) over which to
%                  summarize H per trial (e.g., [0 2] for peri-tone window).
%
%   'KFold'      : # of outer CV folds (default 5)
%
%   'Standardize': whether to z-score features (motif activity) across
%                  trials within each fold (default true)
%
%   'rng'        : RNG seed (default 1)
%
% OUTPUT
%   decRez : struct with fields
%       .acc_mean      : scalar mean accuracy across folds
%       .acc_sem       : SEM of accuracy across folds
%       .acc_folds     : [1 x KFold] accuracy per fold
%       .yhat_trial    : [N x 1] predicted labels (0/1; 1 = Go) per trial
%       .pHat_trial    : [N x 1] predicted P(Go) per trial
%       .beta_folds    : [K x KFold] motif weights per fold
%       .bias_folds    : [1 x KFold] bias per fold
%       .beta_mean     : [K x 1] mean motif weight across folds
%       .motifIdx      : [K x 1] (1:K)
%       .timeWindow    : [1 x 2] time window used
%       .KFold         : scalar
%
% NOTES
%   • This is a **trial-level** decoder: each trial → one K-dim feature
%     vector summarizing H over the specified time window.
%   • Decoder is logistic regression with L2 (ridge) using fitclinear.
%   • Class labels: Go = 1, NoGo = 0.

% ------------------- parse inputs -------------------
p = inputParser;
p.addParameter('timeWindow', [], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('KFold', 5, @(x)isscalar(x)&&x>=2);
p.addParameter('Standardize', true, @(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('rng', 1, @(x)isscalar(x));
p.parse(varargin{:});
opt = p.Results;

[N, K, T] = size(HsY3d); %#ok<NASGU>
if ~isvector(Htime) || numel(Htime) ~= T
    error('Htime must be [1 x T] matching the 3rd dimension of HsY3d.');
end

goI   = trI.goI(:);
nogoI = trI.nogoI(:);
if numel(goI) ~= N || numel(nogoI) ~= N
    error('trI.goI / trI.nogoI must be [N x 1].');
end

% trial labels (0/1)
yTrial = goI;          % Go = 1, NoGo = 0
yTrial = double(yTrial);

% ------------------- choose time window -------------------
if isempty(opt.timeWindow)
    tStart = min(Htime);
    tEnd   = max(Htime);
else
    tStart = opt.timeWindow(1);
    tEnd   = opt.timeWindow(2);
end

tMask = (Htime >= tStart) & (Htime <= tEnd);
if ~any(tMask)
    error('decode_trialType_fromH_logreg: no time bins within specified timeWindow.');
end

% ------------------- summarize H per trial -------------------
% HsY3d: [N x K x T]; take mean over selected time bins → [N x K]
H_feat = nan(N, K);
for i = 1:N
    Hi = squeeze(HsY3d(i, :, tMask));   % [K x Tsel]
    % average across time bins (dim 2)
    H_feat(i,:) = mean(Hi, 2, 'omitnan')';
end

% Optional: you could choose other summary stats (max, etc.) later.

% ------------------- CV folds over trials -------------------
uTrials = (1:N)';   % each trial is one sample
rng(opt.rng);
perm = uTrials(randperm(N));
cuts = round(linspace(0, N, opt.KFold+1));

foldTrials = cell(opt.KFold,1);
for f = 1:opt.KFold
    foldTrials{f} = perm(cuts(f)+1 : cuts(f+1));
end

% ------------------- accumulators -------------------
acc_folds    = nan(1, opt.KFold);
yhat_trial   = nan(N, 1);
pHat_trial   = nan(N, 1);
beta_folds   = nan(K, opt.KFold);
bias_folds   = nan(1, opt.KFold);

for f = 1:opt.KFold
    teTr = foldTrials{f};
    te   = ismember(uTrials, teTr);
    tr   = ~te;

    Ftr = H_feat(tr,:);   % [nTr x K]
    Fte = H_feat(te,:);   % [nTe x K]
    ytr = yTrial(tr);     % [nTr x 1]
    yte = yTrial(te);     % [nTe x 1]

    % standardize across trials on TRAIN only
    if opt.Standardize
        muF = mean(Ftr, 1, 'omitnan');
        sdF = std(Ftr, 0, 1, 'omitnan');
        sdF(~isfinite(sdF) | sdF < 1e-6) = 1e-6;

        FtrZ = (Ftr - muF) ./ sdF;
        FteZ = (Fte - muF) ./ sdF;
    else
        FtrZ = Ftr;
        FteZ = Fte;
    end

    % NaNs in features → 0
    FtrZ(~isfinite(FtrZ)) = 0;
    FteZ(~isfinite(FteZ)) = 0;

    % guard: need both classes in this fold
    if numel(unique(ytr)) < 2
        warning('decode_trialType_fromH_logreg: fold %d has only one class; skipping.', f);
        continue;
    end

    % ----- fit logistic regression -----
    mdl = fitclinear(FtrZ, ytr, ...
        'Learner', 'logistic', ...
        'Regularization', 'ridge', ...
        'Solver', 'lbfgs', ...
        'ClassNames', [0 1]);

    beta = mdl.Beta;   % [K x 1]
    bias = mdl.Bias;   % scalar

    % ----- predict on test trials -----
    scores = FteZ * beta + bias;
    pHat   = 1 ./ (1 + exp(-scores));    % P(Go)
    yhat   = double(pHat >= 0.5);

    acc_folds(f) = mean(yhat == yte);

    % stash
    yhat_trial(te) = yhat;
    pHat_trial(te) = pHat;

    beta_folds(:,f) = beta;
    bias_folds(1,f) = bias;
end

% ------------------- summary stats -------------------
acc_mean = mean(acc_folds, 'omitnan');
acc_sem  = std(acc_folds, 'omitnan') ./ sqrt(sum(isfinite(acc_folds)));

decRez = struct();
decRez.acc_mean    = acc_mean;
decRez.acc_sem     = acc_sem;
decRez.acc_folds   = acc_folds;
decRez.yhat_trial  = yhat_trial;
decRez.pHat_trial  = pHat_trial;
decRez.beta_folds  = beta_folds;
decRez.bias_folds  = bias_folds;
decRez.beta_mean   = mean(beta_folds, 2, 'omitnan');   % [K x 1]
decRez.motifIdx    = (1:K)';                           % motif indices
decRez.timeWindow  = [tStart tEnd];
decRez.KFold       = opt.KFold;

end
