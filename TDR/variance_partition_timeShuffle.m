function out = variance_partition_timeShuffle(Xz, YbigVal, lambdaBest, varargin)
% VARIANCE_PARTITION_TIMESHUFFLE
% Cross-validated explained variance and time-shuffle variance partitioning
% in the style of Musall/Kaufman/Churchland (2019 Nat Neurosci).
%
% ROW / COL NOTATION
%   M  = N * nW  (rows). N = #trials, nW = #time bins (bin-major row order).
%   P  = #predictors (columns of Xz).
%   K  = #motifs / outputs (columns of YbigVal).
%   F  = #CV folds (cv.NumTestSets).
%
% INPUTS
%   Xz         : (M x P) design matrix, standardized (z-scored), bin-major row order
%   YbigVal    : (M x K) responses (motif activities; aligned with Xz rows)
%   lambdaBest : (1 x K) ridge penalty per motif (scalar allowed -> broadcast)
%
% NAME–VALUE OPTIONS
%   'KFold'        : integer, default 5
%   'RowIndex'     : (M x 1) trial index (1..N) per row for *within_trial* shuffles.
%                    For bin-major rows, build as: RowIndex = repmat((1:N)', nW, 1).
%   'Groups'       : cell array of structs with fields:
%                       .name  (string)
%                       .cols  (vector of predictor columns in Xz)
%   'ShuffleMode'  : 'within_trial' | 'global' | 'circular' | 'block'  (default 'within_trial')
%   'BlockSize'    : integer block length for 'block' mode (default 10 rows)
%   'RngSeed'      : scalar RNG seed (default 0)
%   'Verbose'      : true/false (default true)
%
% OUTPUT (struct)
%   .cvR2_full       : (K x 1) full-model CV R^2 per motif
%   .cvR2_global     : scalar, global pooled CV R^2 across all motifs
%   .cvR2_unique     : (P x K) unique ΔR^2 per predictor (full − shuffle-only-this)
%   .cvR2_upper      : (P x K) upper-bound CV R^2 per predictor (keep-only-this; shuffle others)
%   .cvR2_unique_grp : (G x K) unique ΔR^2 per group (if Groups provided)
%   .cvR2_upper_grp  : (G x K) upper-bound CV R^2 per group (if Groups provided)
%   .folds           : cvpartition used
%   .opts            : resolved options
%
% NOTES
%   • “Unique” uses time-shuffling ONLY the variable (or group) of interest.
%   • “Upper-bound” keeps ONLY the variable (or group) time-locked; all others are time-shuffled.
%   • CV baseline uses the *training* mean per fold to compute SST (fold-safe).
%   • Ensure Xz and YbigVal have rows filtered to remove NaN/Inf before calling.

% ---------------- options ----------------
p = inputParser;
p.addParameter('KFold',        5, @(x)isnumeric(x)&&isscalar(x)&&x>=2);
p.addParameter('RowIndex',     [], @(x)isnumeric(x)&&isvector(x));
p.addParameter('Groups',       {}, @(x)iscell(x) || isempty(x));
p.addParameter('ShuffleMode',  'within_trial', @(s)ischar(s) || isstring(s));
p.addParameter('BlockSize',    10, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('RngSeed',      0, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('Verbose',      true, @(x)islogical(x)||ismember(x,[0 1]));
p.parse(varargin{:});
opt = p.Results;

rng(opt.RngSeed);

[M, P] = size(Xz);
K = size(YbigVal, 2);
if isscalar(lambdaBest), lambdaBest = repmat(lambdaBest, 1, K); end
assert(numel(lambdaBest)==K, 'lambdaBest must be length K (one per motif).');
assert(size(YbigVal,1)==M, 'Row count of Xz and YbigVal must match.');

% ----- trial indices for within-trial shuffling -----
if strcmpi(opt.ShuffleMode,'within_trial')
    if isempty(opt.RowIndex)
        warning('RowIndex missing for within_trial shuffle; falling back to global shuffle.');
        opt.ShuffleMode = 'global';
        rowTrialIdx = [];
    else
        rowTrialIdx = opt.RowIndex(:);
        assert(numel(rowTrialIdx)==M, 'RowIndex length must equal #rows of Xz.');
    end
else
    rowTrialIdx = [];
end

% -------------- CV folds --------------
cv = cvpartition(M, 'KFold', opt.KFold);
F = cv.NumTestSets;

% -------------- FULL MODEL cvR2 --------------
if opt.Verbose, fprintf('Fitting FULL model (%d folds)...\n', F); end
[cvR2_full, cvR2_global] = crossval_r2_ridge_all(Xz, YbigVal, lambdaBest, cv);

% -------------- per-predictor UNIQUE & UPPER --------------
if opt.Verbose, fprintf('Computing per-predictor unique ΔR^2 and upper bounds...\n'); end
cvR2_unique = zeros(P,K);
cvR2_upper  = zeros(P,K);

for pcol = 1:P
    % (1) shuffle ONLY this column => reduced model R^2
    X_reduced = Xz;
    X_reduced(:,pcol) = shuffle_column(Xz(:,pcol), opt.ShuffleMode, rowTrialIdx, opt.BlockSize);
    for k = 1:K
        r2_reduced = crossval_r2_ridge_foldsafe(X_reduced, YbigVal(:,k), lambdaBest(k), cv);
        cvR2_unique(pcol,k) = cvR2_full(k) - r2_reduced;
    end

    % (2) shuffle ALL OTHER columns => keep-only-this upper bound
    X_upper = Xz;
    others = setdiff(1:P, pcol);
    X_upper(:,others) = shuffle_matrix(X_upper(:,others), opt.ShuffleMode, rowTrialIdx, opt.BlockSize);
    for k = 1:K
        cvR2_upper(pcol,k) = crossval_r2_ridge_foldsafe(X_upper, YbigVal(:,k), lambdaBest(k), cv);
    end
end

% -------------- per-group UNIQUE & UPPER --------------
cvR2_unique_grp = [];
cvR2_upper_grp  = [];
if ~isempty(opt.Groups)
    G = numel(opt.Groups);
    cvR2_unique_grp = zeros(G,K);
    cvR2_upper_grp  = zeros(G,K);
    if opt.Verbose, fprintf('Computing per-GROUP unique ΔR^2 and upper bounds...\n'); end

    for g = 1:G
        cols = opt.Groups{g}.cols(:)';     % indices into Xz
        cols = cols(cols>=1 & cols<=P);
        if isempty(cols), continue; end

        % Shuffle ONLY this group (reduced)
        X_reduced = Xz;
        X_reduced(:,cols) = shuffle_matrix(X_reduced(:,cols), opt.ShuffleMode, rowTrialIdx, opt.BlockSize);
        for k = 1:K
            r2_reduced = crossval_r2_ridge_foldsafe(X_reduced, YbigVal(:,k), lambdaBest(k), cv);
            cvR2_unique_grp(g,k) = cvR2_full(k) - r2_reduced;
        end

        % Keep ONLY this group's timing (upper)
        X_keep = Xz;
        keepMask = false(1,P); keepMask(cols) = true;
        X_keep(:,~keepMask) = shuffle_matrix(X_keep(:,~keepMask), opt.ShuffleMode, rowTrialIdx, opt.BlockSize);
        for k = 1:K
            cvR2_upper_grp(g,k) = crossval_r2_ridge_foldsafe(X_keep, YbigVal(:,k), lambdaBest(k), cv);
        end
    end
end

% -------------- pack output --------------
out = struct();
out.cvR2_full       = cvR2_full;            % K x 1
out.cvR2_global     = cvR2_global;          % scalar
out.cvR2_unique     = cvR2_unique;          % P x K
out.cvR2_upper      = cvR2_upper;           % P x K
out.cvR2_unique_grp = cvR2_unique_grp;      % G x K (if groups provided)
out.cvR2_upper_grp  = cvR2_upper_grp;       % G x K (if groups provided)
out.folds           = cv;
out.opts            = opt;
end

% ======== helpers ========

function r2 = crossval_r2_ridge_foldsafe(X, y, lambda, cv)
% Mean CV R^2 across folds using TRAINING-mean baseline for SST (fold-safe).
F = cv.NumTestSets;
r2fold = zeros(F,1);
p = size(X,2);
I = speye(p);

for f = 1:F
    tr = training(cv,f);
    te = test(cv,f);

    Xtr = X(tr,:); ytr = y(tr);
    Xte = X(te,:); yte = y(te);

    % Ridge fit (closed form)
    beta = (Xtr' * Xtr + lambda * I) \ (Xtr' * ytr);

    % Predict & compute fold-safe R^2
    yhat = Xte * beta;
    ss_res = sum((yte - yhat).^2);
    mu_tr  = mean(ytr);
    ss_tot = sum((yte - mu_tr).^2);
    r2fold(f) = 1 - ss_res / max(ss_tot, eps);
end
r2 = mean(r2fold);
end

function xsh = shuffle_column(x, mode, trialIdx, blockSize)
% Shuffle one column according to the chosen policy.
switch lower(string(mode))
    case "global"
        xsh = x(randperm(numel(x)));

    case "within_trial"
        assert(~isempty(trialIdx), 'RowIndex required for within_trial mode.');
        xsh = x;
        uTrials = unique(trialIdx(:))';
        for t = uTrials
            rows = find(trialIdx==t);
            xsh(rows) = x(rows(randperm(numel(rows))));
        end

    case "circular"
        % single random circular shift
        sh = randi([0, numel(x)-1], 1, 1);
        xsh = circshift(x, sh);

    case "block"
        % shuffle chunks of length blockSize (approximate)
        n = numel(x);
        idx = 1:blockSize:n;
        blocks = arrayfun(@(i) i:min(i+blockSize-1,n), idx,'uni',0);
        order = randperm(numel(blocks));
        xsh = zeros(n,1);
        pos = 1;
        for k = 1:numel(order)
            b = blocks{order(k)};
            m = numel(b);
            xsh(pos:pos+m-1) = x(b);
            pos = pos + m;
        end

    otherwise
        error('Unknown ShuffleMode: %s', mode);
end
end

function Xsh = shuffle_matrix(X, mode, trialIdx, blockSize)
% Apply the same shuffling policy column-wise.
Xsh = X;
for j = 1:size(X,2)
    Xsh(:,j) = shuffle_column(X(:,j), mode, trialIdx, blockSize);
end
end
