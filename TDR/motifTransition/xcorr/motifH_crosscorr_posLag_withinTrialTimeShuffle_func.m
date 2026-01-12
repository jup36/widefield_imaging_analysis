function S = motifH_crosscorr_posLag_withinTrialTimeShuffle_func(filePath, fileKeyword, varargin)
% MOTIFH_CROSSCORR_POSLAG_WITHINTRIALTIMESHUFFLE_FUNC
%   Observed + within-trial time-shuffle null for motif–motif cross-correlation.
%
%   S = motifH_crosscorr_posLag_withinTrialTimeShuffle_func(filePath, fileKeyword, ...)
%
% INPUTS
%   filePath    : full path to session "task" folder
%   fileKeyword : suffix for H-file
%
% NAME–VALUE PAIRS
%   'postLagPeriod'   : seconds, default 0.5
%   'nShuffle'        : integer, default 1000
%   'shuffleMethod'   : 'circshift' (default) or 'permute'
%   'rngSeed'         : [] (default) or scalar seed
%   'zscoreHs'        : true (default)
%   'useSymmetry'     : true (default) ~2x speedup for xcorr matrix
%
% OUTPUT
%   S is a struct containing observed matrices, shuffle distributions, and stats.

%% -------------------- Parse inputs --------------------
p = inputParser;
p.addParameter('postLagPeriod', 0.5, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('nShuffle', 1000, @(x) isnumeric(x) && isscalar(x) && x>=1 && x==round(x));
p.addParameter('shuffleMethod', 'circshift', @(s) ischar(s) || isstring(s));
p.addParameter('rngSeed', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('zscoreHs', true, @(x) islogical(x) && isscalar(x));
p.addParameter('useSymmetry', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});

postLagPeriod_sec = p.Results.postLagPeriod;
nShuffle          = p.Results.nShuffle;
shuffleMethod     = lower(string(p.Results.shuffleMethod));
rngSeed           = p.Results.rngSeed;
doZscore          = p.Results.zscoreHs;
useSymmetry       = p.Results.useSymmetry;

if ~ismember(shuffleMethod, ["circshift","permute"])
    error('shuffleMethod must be ''circshift'' or ''permute''.');
end
if ~isempty(rngSeed)
    rng(rngSeed);
end

%% -------------------- 0) Grab files --------------------
header      = extract_date_animalID_header(filePath);
keyword_beh = '_alignedPupilOrofacial.mat';

filePath_matfiles = cell2mat(GrabFiles_sort_trials('Matfiles', 0, {filePath}));
filePath_H        = cell2mat(GrabFiles_sort_trials([header '*' fileKeyword], 0, {filePath_matfiles}));
filePath_B        = cell2mat(GrabFiles_sort_trials([header '*' keyword_beh], 0, {filePath_matfiles}));

load(filePath_H, 'tbytDat_hAligned');
load(filePath_B, 'tbytDat');

trI = trialTypeInfoAuditoryGngTbytDat(tbytDat);
assert(size(tbytDat_hAligned, 2) == numel(tbytDat), 'Mismatch in trial counts.');

%% -------------------- 1) Stack trials --------------------
Hs = stack_trials_H(tbytDat_hAligned, 'zscore', doZscore);
[~, K, T] = size(Hs.Y3);

stepSec = Hs.params.Step;
nBin = round(postLagPeriod_sec / stepSec);   % seconds → bins

%% -------------------- 2) Observed cross-correlograms --------------------
[XcorrMat_sess, lags] = computeMotifXcorr_fast(Hs.Y3, nBin, useSymmetry);

if any(trI.goI)
    X_go = Hs.Y3(trI.goI, :, :);
    XcorrMat_sess_go = computeMotifXcorr_fast(X_go, nBin, useSymmetry);
else
    XcorrMat_sess_go = nan(K, K, 2*nBin + 1);
end

if any(trI.nogoI)
    X_ng = Hs.Y3(trI.nogoI, :, :);
    XcorrMat_sess_nogo = computeMotifXcorr_fast(X_ng, nBin, useSymmetry);
else
    XcorrMat_sess_nogo = nan(K, K, 2*nBin + 1);
end

posLagMask = (lags > 0 & lags <= nBin);

xcorrPosLagMat      = squeeze(mean(XcorrMat_sess(:, :,      posLagMask), 3));
xcorrPosLagMat_go   = squeeze(mean(XcorrMat_sess_go(:, :,   posLagMask), 3));
xcorrPosLagMat_nogo = squeeze(mean(XcorrMat_sess_nogo(:, :, posLagMask), 3));

%% -------------------- 3) Shuffle null distributions --------------------
shuf_all  = nan(K, K, nShuffle);
shuf_go   = nan(K, K, nShuffle);
shuf_nogo = nan(K, K, nShuffle);

X_all = Hs.Y3;

X_go = [];
X_ng = [];
if any(trI.goI),   X_go = Hs.Y3(trI.goI, :, :); end
if any(trI.nogoI), X_ng = Hs.Y3(trI.nogoI, :, :); end

for s = 1:nShuffle
    % ALL
    Xs_all = withinTrialShuffle_independent(X_all, shuffleMethod);
    Xc_all = computeMotifXcorr_fast(Xs_all, nBin, useSymmetry);
    shuf_all(:, :, s) = squeeze(mean(Xc_all(:, :, posLagMask), 3));

    % GO
    if ~isempty(X_go)
        Xs_go = withinTrialShuffle_independent(X_go, shuffleMethod);
        Xc_go = computeMotifXcorr_fast(Xs_go, nBin, useSymmetry);
        shuf_go(:, :, s) = squeeze(mean(Xc_go(:, :, posLagMask), 3));
    end

    % NOGO
    if ~isempty(X_ng)
        Xs_ng = withinTrialShuffle_independent(X_ng, shuffleMethod);
        Xc_ng = computeMotifXcorr_fast(Xs_ng, nBin, useSymmetry);
        shuf_nogo(:, :, s) = squeeze(mean(Xc_ng(:, :, posLagMask), 3));
    end
end

%% -------------------- 4) Stats vs null --------------------
S = struct();

S.meta = struct('header', header);

S.params = struct();
S.params.postLagPeriod_sec = postLagPeriod_sec;
S.params.stepSec           = stepSec;
S.params.nBin              = nBin;
S.params.lags              = lags;
S.params.posLagMask        = posLagMask;
S.params.nShuffle          = nShuffle;
S.params.shuffleMethod     = shuffleMethod;
S.params.rngSeed           = rngSeed;
S.params.zscoreHs          = doZscore;
S.params.useSymmetry       = useSymmetry;

S.obs = struct();
S.obs.XcorrMat_sess      = XcorrMat_sess;
S.obs.XcorrMat_sess_go   = XcorrMat_sess_go;
S.obs.XcorrMat_sess_nogo = XcorrMat_sess_nogo;

S.obs.xcorrPosLagMat      = xcorrPosLagMat;
S.obs.xcorrPosLagMat_go   = xcorrPosLagMat_go;
S.obs.xcorrPosLagMat_nogo = xcorrPosLagMat_nogo;

S.shuf = struct();
S.shuf.xcorrPosLagMat_all  = shuf_all;
S.shuf.xcorrPosLagMat_go   = shuf_go;
S.shuf.xcorrPosLagMat_nogo = shuf_nogo;

[S.shuf.mean_all,  S.shuf.std_all,  S.shuf.z_all,  S.shuf.p_all]  = shufStats(xcorrPosLagMat,      shuf_all);
[S.shuf.mean_go,   S.shuf.std_go,   S.shuf.z_go,   S.shuf.p_go]   = shufStats(xcorrPosLagMat_go,   shuf_go);
[S.shuf.mean_nogo, S.shuf.std_nogo, S.shuf.z_nogo, S.shuf.p_nogo] = shufStats(xcorrPosLagMat_nogo, shuf_nogo);

end

%% ========================================================================
function Xs = withinTrialShuffle_independent(X, method)
% Independent within-trial shuffle per trial×motif (NOT a shared shift).
% X: [N x K x T]
[N, K, T] = size(X);
Xs = X;

switch method
    case "circshift"
        % sh(n,k) ~ Uniform{0,...,T-1}
        shMat = randi(T, N, K) - 1;
        for n = 1:N
            for k = 1:K
                Xs(n, k, :) = circshift(squeeze(X(n, k, :)), shMat(n, k));
            end
        end

    case "permute"
        for n = 1:N
            for k = 1:K
                permIdx = randperm(T);
                Xs(n, k, :) = X(n, k, permIdx);
            end
        end

    otherwise
        error('Unknown shuffle method.');
end
end

%% ------------------------------------------------------------------------
function [mu, sig, z, p] = shufStats(obsMat, shufMat)
% obsMat  : [K x K]
% shufMat : [K x K x S]
mu  = mean(shufMat, 3, 'omitnan');
sig = std(shufMat, 0, 3, 'omitnan');

sig(sig == 0) = eps;   % guard
z = (obsMat - mu) ./ sig;

% Empirical two-sided p:
% p = (1 + #{|shuf - mu| >= |obs - mu|}) / (S + 1)
S = size(shufMat, 3);
devObs  = abs(obsMat - mu);
devShuf = abs(shufMat - mu);  % implicit expansion
countExtreme = sum(devShuf >= devObs, 3, 'omitnan');
p = (1 + countExtreme) ./ (S + 1);
end

%% ------------------------------------------------------------------------
function [XcorrMat, lags] = computeMotifXcorr_fast(HsY3, maxLag, useSymmetry)
% COMPUTEMOTIFXCORR_FAST
%   XcorrMat(i,j,:) = xcorr(Xi, Xj, maxLag, 'coeff')
%   where Xi, Xj are flattened across trials/time.
%
%   If useSymmetry=true, compute only upper triangle and fill lower using:
%     xcorr(Xj, Xi) = flip(xcorr(Xi, Xj))  (for all lags)
%
% INPUT
%   HsY3   : [N x K x T]
%   maxLag : nonnegative integer
%
% OUTPUT
%   XcorrMat : [K x K x (2*maxLag+1)]
%   lags     : [-maxLag: maxLag]

if nargin < 2 || isempty(maxLag), maxLag = 10; end
if nargin < 3 || isempty(useSymmetry), useSymmetry = true; end
if maxLag < 0 || maxLag ~= round(maxLag)
    error('maxLag must be a nonnegative integer.');
end

[~, K, ~] = size(HsY3);
L = 2*maxLag + 1;
lags = -maxLag:maxLag;

XcorrMat = zeros(K, K, L);

% Flatten each motif once
Xflat = cell(1, K);
for k = 1:K
    xk = reshape(HsY3(:, k, :), [], 1);
    xk(~isfinite(xk)) = 0;
    Xflat{k} = xk;
end

if useSymmetry
    for i = 1:K
        Xi = Xflat{i};
        for j = i:K
            Xj = Xflat{j};
            cij = xcorr(Xi, Xj, maxLag, 'coeff');  % [L x 1]
            XcorrMat(i, j, :) = cij;
            if j ~= i
                XcorrMat(j, i, :) = flipud(cij);    % symmetry fill
            end
        end
    end
else
    for i = 1:K
        Xi = Xflat{i};
        for j = 1:K
            Xj = Xflat{j};
            XcorrMat(i, j, :) = xcorr(Xi, Xj, maxLag, 'coeff');
        end
    end
end

end
