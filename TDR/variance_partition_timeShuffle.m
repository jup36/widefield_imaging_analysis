function out = variance_partition_timeShuffle(X, Y, lambda, rowTrialIdx, varargin)
% VARIANCE_PARTITION_TIMESHUFFLE
% Cross-validated explained variance and time-shuffle variance partitioning
% in the style of Musall/Kaufman/Churchland (2019 Nat Neurosci).
%
% ROW / COL NOTATION
%   M  = N * nW  (rows). N = #trials, nW = #time bins (bin-major row order).
%   P  = #predictors (columns of X).
%   K  = #motifs / outputs (columns of Y).
%
% INPUTS
%   X         : (M x P) design matrix, standardized (z-scored), bin-major row order
%   Y    : (M x K) responses (motif activities; aligned with X rows)
%   lambda : (1 x K) ridge penalty per motif (scalar allowed -> broadcast)
%
% NAME–VALUE OPTIONS
%   'KFold'        : integer, default 5
%   'Groups'       : cell array of structs with fields:
%                       .name  (string)
%                       .cols  (vector of predictor columns in X)
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

%   .opts            : resolved options
%
% NOTES
%   • “Unique” uses time-shuffling ONLY the variable (or group) of interest.
%   • “Upper-bound” keeps ONLY the variable (or group) time-locked; all others are time-shuffled.
%   • CV baseline uses the *training* mean per fold to compute SST (fold-safe).
%   • Ensure X and Y have rows filtered to remove NaN/Inf before calling.

% ---------------- options ----------------
p = inputParser;
p.addParameter('KFold',        5, @(x)isnumeric(x)&&isscalar(x)&&x>=2);
p.addParameter('Groups',       {}, @(x)iscell(x) || isempty(x));
p.addParameter('ShuffleMode',  'within_trial', @(s)ischar(s) || isstring(s));
p.addParameter('BlockSize',    10, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('RngSeed',      0, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('Verbose',      true, @(x)islogical(x)||ismember(x,[0 1]));
p.parse(varargin{:});
opt = p.Results;

rng(opt.RngSeed);

[M, P] = size(X);
K = size(Y, 2);
if isscalar(lambda), lambda = repmat(lambda, 1, K); end
assert(numel(lambda)==K, 'lambda must be length K (one per motif).');
assert(size(Y,1)==M, 'Row count of X and Y must match.');

% -------------- FULL MODEL cvR2 --------------
if opt.Verbose, fprintf('Fitting FULL model...\n'); end
[cvR2_full, cvR2_global] = cv_global_r2_ridge(X, Y, rowTrialIdx, 'lambda',lambda, 'outerFolds', opt.KFold);

% -------------- per-predictor UNIQUE & UPPER --------------
if opt.Verbose, fprintf('Computing per-predictor unique ΔR^2 and upper bounds...\n'); end
cvR2_unique = zeros(P,K);
cvR2_upper  = zeros(P,K);

for pcol = 1:P
    % (1) shuffle ONLY this column => reduced model R^2
    X_reduced = X;
    X_reduced(:,pcol) = shuffle_column(X(:,pcol), opt.ShuffleMode, rowTrialIdx, opt.BlockSize);
    cvR2_reduced = cv_global_r2_ridge(X_reduced, Y, rowTrialIdx, 'lambda',lambda, 'outerFolds', opt.KFold);
    cvR2_unique(pcol, :) = max(cvR2_full-cvR2_reduced, 0); 

    % (2) shuffle ALL OTHER columns => keep-only-this upper bound
    X_upper = X;
    others = setdiff(1:P, pcol);
    X_upper(:,others) = shuffle_matrix(X_upper(:,others), opt.ShuffleMode, rowTrialIdx, opt.BlockSize); 
    cvR2_keep = cv_global_r2_ridge(X_upper, Y,  rowTrialIdx, 'lambda',lambda, 'outerFolds', opt.KFold);
    cvR2_upper(pcol,:) = max(cvR2_keep, 0); 
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
        cols = opt.Groups{g}.cols(:)';     % indices into X
        cols = cols(cols>=1 & cols<=P);
        if isempty(cols), continue; end

        % Shuffle ONLY this group (reduced)
        X_reduced = X;
        X_reduced(:,cols) = shuffle_matrix(X_reduced(:,cols), opt.ShuffleMode, rowTrialIdx, opt.BlockSize);
        cvR2_reduced_group = cv_global_r2_ridge(X_reduced, Y,  rowTrialIdx, 'lambda',lambda, 'outerFolds', opt.KFold);
        cvR2_unique_grp(g,:) = max(cvR2_full - cvR2_reduced_group, 0);
  
        % Keep ONLY this group's timing (upper)
        X_keep = X;
        keepMask = false(1,P); keepMask(cols) = true;
        X_keep(:,~keepMask) = shuffle_matrix(X_keep(:,~keepMask), opt.ShuffleMode, rowTrialIdx, opt.BlockSize);
        cvR2_upper_group = cv_global_r2_ridge(X_keep, Y, rowTrialIdx, 'lambda',lambda, 'outerFolds', opt.KFold);
        cvR2_upper_grp(g,:) = max(cvR2_upper_group, 0); 
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
out.opts            = opt;
end



% ======== helpers ========
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
