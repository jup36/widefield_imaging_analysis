function S = motifH_perTrial_crosscorr_posLag_withinTrialTimeShuffle_func(filePath, fileKeyword, varargin)
% MOTIFH_PERTRIAL_CROSSCORR_POSLAG_WITHINTRIALTIMESHUFFLE_FUNC
%   Per-trial motif–motif xcorr (full lag window) averaged across trials using Fisher-z,
%   with optional within-trial time-shuffle null. Also computes a pooled positive-lag
%   summary over a user-defined positive-lag epoch (e.g., 0–0.5 s) that can be shorter
%   than the full lag window (e.g., ±1 s).
%
%   S = motifH_perTrial_crosscorr_posLag_withinTrialTimeShuffle_func(filePath, fileKeyword, ...)
%
% NAME–VALUE PAIRS (explicit)
%   'xcorrLagWindowSec'     : total lag window for xcorr curve in seconds (default 1.0)
%   'posLagPoolWindowSec'   : positive-lag pooling window in seconds (default 0.5)
%   'doTimeShuffle'         : true (default). If false, skip shuffles and return observed only.
%   'nShuffle'              : integer, default 1000
%   'shuffleMethod'         : 'circshift' (default) or 'permute'
%   'rngSeed'               : [] (default) or scalar seed
%   'zscoreHs'              : true (default)
%   'useSymmetry'           : true (default)
%   'doFisherZ'             : true (default)
%   'clipR'                 : scalar in (0,1), default 0.999
%   'showProgress'          : true (default). Prints progress during shuffles (works w/ parfor).
%   'progressEvery'         : positive integer, default 25. Print every N shuffles.
%
% OUTPUT
%   See prior version + progress printing if enabled.

%% -------------------- Parse inputs --------------------
p = inputParser;

% Explicit lag windows
p.addParameter('xcorrLagWindowSec',   1.0, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('posLagPoolWindowSec', 0.5, @(x) isnumeric(x) && isscalar(x) && x>0);

% Shuffle + preprocessing options
p.addParameter('doTimeShuffle', true, @(x) islogical(x) && isscalar(x));
p.addParameter('nShuffle', 1000, @(x) isnumeric(x) && isscalar(x) && x>=1 && x==round(x));
p.addParameter('shuffleMethod', 'circshift', @(s) ischar(s) || isstring(s));
p.addParameter('rngSeed', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('zscoreHs', true, @(x) islogical(x) && isscalar(x));
p.addParameter('useSymmetry', true, @(x) islogical(x) && isscalar(x));
p.addParameter('doFisherZ', true, @(x) islogical(x) && isscalar(x));
p.addParameter('clipR', 0.999, @(x) isnumeric(x) && isscalar(x) && x>0 && x<1);

% NEW: progress tracking
p.addParameter('showProgress', true, @(x) islogical(x) && isscalar(x));
p.addParameter('progressEvery', 25, @(x) isnumeric(x) && isscalar(x) && x>=1 && x==round(x));

p.parse(varargin{:});

xcorrLagWindowSec   = p.Results.xcorrLagWindowSec;
posLagPoolWindowSec = p.Results.posLagPoolWindowSec;

doTimeShuffle       = p.Results.doTimeShuffle;
nShuffle            = p.Results.nShuffle;
shuffleMethod       = lower(string(p.Results.shuffleMethod));
rngSeed             = p.Results.rngSeed;
doZscore            = p.Results.zscoreHs;
useSymmetry         = p.Results.useSymmetry;
doFisherZ           = p.Results.doFisherZ;
clipR               = p.Results.clipR;

showProgress        = p.Results.showProgress;
progressEvery       = p.Results.progressEvery;

if posLagPoolWindowSec > xcorrLagWindowSec
    error('posLagPoolWindowSec (%.3f) must be <= xcorrLagWindowSec (%.3f).', ...
        posLagPoolWindowSec, xcorrLagWindowSec);
end

if doTimeShuffle
    if ~ismember(shuffleMethod, ["circshift","permute"])
        error('shuffleMethod must be ''circshift'' or ''permute''.');
    end
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
[~, K, ~] = size(Hs.Y3);

stepSec = Hs.params.Step;

% Full xcorr lag window (curve)
nBin_total = round(xcorrLagWindowSec / stepSec);     % seconds → bins
lags = -nBin_total:nBin_total;
L = numel(lags);

% Positive-lag pooling window (summary)
nBin_pool = round(posLagPoolWindowSec / stepSec);    % seconds → bins
posLagMask_pool = (lags > 0 & lags <= nBin_pool);

%% -------------------- 2) Observed per-trial averaged xcorr --------------------
X_all = Hs.Y3;

XcorrMat_sess = computeMotifXcorr_perTrial_FisherZ(X_all, nBin_total, useSymmetry, doFisherZ, clipR);

if any(trI.goI)
    X_go = Hs.Y3(trI.goI, :, :);
    XcorrMat_sess_go = computeMotifXcorr_perTrial_FisherZ(X_go, nBin_total, useSymmetry, doFisherZ, clipR);
else
    XcorrMat_sess_go = nan(K, K, L);
end

if any(trI.nogoI)
    X_ng = Hs.Y3(trI.nogoI, :, :);
    XcorrMat_sess_nogo = computeMotifXcorr_perTrial_FisherZ(X_ng, nBin_total, useSymmetry, doFisherZ, clipR);
else
    XcorrMat_sess_nogo = nan(K, K, L);
end

% Pooled pos-lag summary uses ONLY 0–posLagPoolWindowSec
xcorrPosLagMat      = squeeze(mean(XcorrMat_sess(:, :,      posLagMask_pool), 3, 'omitnan'));
xcorrPosLagMat_go   = squeeze(mean(XcorrMat_sess_go(:, :,   posLagMask_pool), 3, 'omitnan'));
xcorrPosLagMat_nogo = squeeze(mean(XcorrMat_sess_nogo(:, :, posLagMask_pool), 3, 'omitnan'));

%% -------------------- 3) Package observed outputs --------------------
S = struct();
S.meta = struct('header', header);

S.params = struct();
S.params.xcorrLagWindowSec     = xcorrLagWindowSec;
S.params.posLagPoolWindowSec   = posLagPoolWindowSec;
S.params.stepSec               = stepSec;
S.params.nBin_total            = nBin_total;
S.params.nBin_pool             = nBin_pool;
S.params.lags                  = lags;
S.params.posLagMask_pool       = posLagMask_pool;

S.params.doTimeShuffle         = doTimeShuffle;
S.params.shuffleMethod         = shuffleMethod;
S.params.nShuffle              = nShuffle;
S.params.rngSeed               = rngSeed;
S.params.zscoreHs              = doZscore;
S.params.useSymmetry           = useSymmetry;
S.params.doFisherZ             = doFisherZ;
S.params.clipR                 = clipR;

S.params.showProgress          = showProgress;
S.params.progressEvery         = progressEvery;

S.obs = struct();
S.obs.XcorrMat_sess      = XcorrMat_sess;
S.obs.XcorrMat_sess_go   = XcorrMat_sess_go;
S.obs.XcorrMat_sess_nogo = XcorrMat_sess_nogo;

S.obs.xcorrPosLagMat      = xcorrPosLagMat;
S.obs.xcorrPosLagMat_go   = xcorrPosLagMat_go;
S.obs.xcorrPosLagMat_nogo = xcorrPosLagMat_nogo;

%% -------------------- 4) Optional shuffle null + stats --------------------
if ~doTimeShuffle
    return;
end

% pooled samples for p-values
shuf_all  = nan(K, K, nShuffle);
shuf_go   = nan(K, K, nShuffle);
shuf_nogo = nan(K, K, nShuffle);

% per-lag descriptives (mean/std) via sums
sum_all    = zeros(K, K, L);   sumsq_all   = zeros(K, K, L);   cnt_all   = zeros(K, K, L);
sum_go     = zeros(K, K, L);   sumsq_go    = zeros(K, K, L);   cnt_go    = zeros(K, K, L);
sum_nogo   = zeros(K, K, L);   sumsq_nogo  = zeros(K, K, L);   cnt_nogo  = zeros(K, K, L);

% Trial subsets
X_go = [];
X_ng = [];
if any(trI.goI),   X_go = Hs.Y3(trI.goI, :, :); end
if any(trI.nogoI), X_ng = Hs.Y3(trI.nogoI, :, :); end

% Initialize RNG on each worker (reduces chance of identical shuffles)
try
    pctRunOnAll rng('shuffle');
catch
end

% ---------------- Progress (parfor-safe) ----------------
dq = [];
tStart = tic;
if showProgress
    dq = parallel.pool.DataQueue;
    nDone = 0; % on client
    afterEach(dq, @updateProgress);
end

parfor s = 1:nShuffle
    % ---------------- ALL ----------------
    Xs_all = withinTrialShuffle_independent(X_all, shuffleMethod);
    Xc_all = computeMotifXcorr_perTrial_FisherZ(Xs_all, nBin_total, useSymmetry, doFisherZ, clipR);

    shuf_all(:, :, s) = squeeze(mean(Xc_all(:, :, posLagMask_pool), 3, 'omitnan'));

    m = isfinite(Xc_all);
    Xc0 = Xc_all; Xc0(~m) = 0;
    sum_all   = sum_all   + Xc0;
    sumsq_all = sumsq_all + Xc0.^2;
    cnt_all   = cnt_all   + double(m);

    % ---------------- GO ----------------
    if ~isempty(X_go)
        Xs_go = withinTrialShuffle_independent(X_go, shuffleMethod);
        Xc_go = computeMotifXcorr_perTrial_FisherZ(Xs_go, nBin_total, useSymmetry, doFisherZ, clipR);

        shuf_go(:, :, s) = squeeze(mean(Xc_go(:, :, posLagMask_pool), 3, 'omitnan'));

        m = isfinite(Xc_go);
        Xc0 = Xc_go; Xc0(~m) = 0;
        sum_go   = sum_go   + Xc0;
        sumsq_go = sumsq_go + Xc0.^2;
        cnt_go   = cnt_go   + double(m);
    end

    % ---------------- NOGO ----------------
    if ~isempty(X_ng)
        Xs_ng = withinTrialShuffle_independent(X_ng, shuffleMethod);
        Xc_ng = computeMotifXcorr_perTrial_FisherZ(Xs_ng, nBin_total, useSymmetry, doFisherZ, clipR);

        shuf_nogo(:, :, s) = squeeze(mean(Xc_ng(:, :, posLagMask_pool), 3, 'omitnan'));

        m = isfinite(Xc_ng);
        Xc0 = Xc_ng; Xc0(~m) = 0;
        sum_nogo   = sum_nogo   + Xc0;
        sumsq_nogo = sumsq_nogo + Xc0.^2;
        cnt_nogo   = cnt_nogo   + double(m);
    end

    % Progress ping
    if showProgress
        send(dq, s);
    end
end

% Finalize per-lag null mean/std
[mu_all,  sd_all]  = finalizeMeanStd(sum_all,  sumsq_all,  cnt_all);
[mu_go,   sd_go]   = finalizeMeanStd(sum_go,   sumsq_go,   cnt_go);
[mu_nogo, sd_nogo] = finalizeMeanStd(sum_nogo, sumsq_nogo, cnt_nogo);

S.shuf = struct();

S.shuf.xcorrPosLagMat_all  = shuf_all;
S.shuf.xcorrPosLagMat_go   = shuf_go;
S.shuf.xcorrPosLagMat_nogo = shuf_nogo;

[S.shuf.mean_all,  S.shuf.std_all,  S.shuf.z_all,  S.shuf.p_all]  = shufStats(xcorrPosLagMat,      shuf_all);
[S.shuf.mean_go,   S.shuf.std_go,   S.shuf.z_go,   S.shuf.p_go]   = shufStats(xcorrPosLagMat_go,   shuf_go);
[S.shuf.mean_nogo, S.shuf.std_nogo, S.shuf.z_nogo, S.shuf.p_nogo] = shufStats(xcorrPosLagMat_nogo, shuf_nogo);

S.shuf.curve = struct();
S.shuf.curve.mean_all  = mu_all;
S.shuf.curve.std_all   = sd_all;
S.shuf.curve.mean_go   = mu_go;
S.shuf.curve.std_go    = sd_go;
S.shuf.curve.mean_nogo = mu_nogo;
S.shuf.curve.std_nogo  = sd_nogo;

if showProgress
    fprintf('[%s] Completed %d/%d shuffles in %.1f sec\n', header, nShuffle, nShuffle, toc(tStart));
end

%% ---------------- nested: progress callback (client) ----------------
    function updateProgress(~)
        nDone = nDone + 1;
        if mod(nDone, progressEvery) == 0 || nDone == 1 || nDone == nShuffle
            fprintf('[%s] Shuffle progress: %d/%d (%.1f%%) elapsed %.1fs\n', ...
                header, nDone, nShuffle, 100*nDone/nShuffle, toc(tStart));
        end
    end

end

%% ========================================================================
function Xs = withinTrialShuffle_independent(X, method)
% Independent within-trial shuffle per trial×motif (NOT a shared shift).
% X: [N x K x T]
[N, K, T] = size(X);
Xs = X;

switch method
    case "circshift"
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
function [mu, sd] = finalizeMeanStd(sumX, sumsqX, cntX)
mu = nan(size(sumX));
sd = nan(size(sumX));

valid = cntX > 0;
mu(valid) = sumX(valid) ./ cntX(valid);

ex2 = nan(size(sumX));
ex2(valid) = sumsqX(valid) ./ cntX(valid);

v = ex2 - mu.^2;
v(v < 0 & v > -1e-12) = 0; % numeric guard
sd(valid) = sqrt(v(valid));
end

%% ------------------------------------------------------------------------
function [mu, sig, z, p] = shufStats(obsMat, shufMat)
% obsMat  : [K x K]
% shufMat : [K x K x S]
mu  = mean(shufMat, 3, 'omitnan');
sig = std(shufMat, 0, 3, 'omitnan');

sig(sig == 0) = eps;
z = (obsMat - mu) ./ sig;

S = size(shufMat, 3);
devObs  = abs(obsMat - mu);
devShuf = abs(shufMat - mu);
countExtreme = sum(devShuf >= devObs, 3, 'omitnan');
p = (1 + countExtreme) ./ (S + 1);
end

%% ------------------------------------------------------------------------
function XcorrMat = computeMotifXcorr_perTrial_FisherZ(HsY3, maxLag, useSymmetry, doFisherZ, clipR)
% COMPUTEMOTIFXCORR_PERTRIAL_FISHERZ
%   Compute per-trial xcorr for each motif pair, then average across trials.
%   Optionally uses Fisher-z averaging (atanh -> mean -> tanh).

if nargin < 2 || isempty(maxLag), maxLag = 10; end
if nargin < 3 || isempty(useSymmetry), useSymmetry = true; end
if nargin < 4 || isempty(doFisherZ), doFisherZ = true; end
if nargin < 5 || isempty(clipR), clipR = 0.999; end

if maxLag < 0 || maxLag ~= round(maxLag)
    error('maxLag must be a nonnegative integer.');
end

[~, K, ~] = size(HsY3);
L = 2*maxLag + 1;

Acc  = zeros(K, K, L);
nEff = zeros(K, K, L);

for n = 1:size(HsY3,1)
    Xn = squeeze(HsY3(n, :, :)); % [K x T]
    if any(~isfinite(Xn(:)))
        Xn(~isfinite(Xn)) = 0;
    end

    if useSymmetry
        for i = 1:K
            xi = Xn(i, :).';
            xi = xi - mean(xi, 'omitnan'); % mean-subtract per trial
            for j = i:K
                xj = Xn(j, :).';
                xj = xj - mean(xj, 'omitnan');

                r = xcorr(xi, xj, maxLag, 'coeff');  % [L x 1]
                r = max(min(r, clipR), -clipR);

                if doFisherZ
                    v = atanh(r);
                else
                    v = r;
                end

                finiteMask = isfinite(v);
                Acc(i,j,finiteMask)  = Acc(i,j,finiteMask)  + reshape(v(finiteMask), 1,1,[]);
                nEff(i,j,finiteMask) = nEff(i,j,finiteMask) + 1;

                if j ~= i
                    vflip = flipud(v);
                    finiteMask2 = isfinite(vflip);
                    Acc(j,i,finiteMask2)  = Acc(j,i,finiteMask2)  + reshape(vflip(finiteMask2), 1,1,[]);
                    nEff(j,i,finiteMask2) = nEff(j,i,finiteMask2) + 1;
                end
            end
        end
    else
        for i = 1:K
            xi = Xn(i, :).';
            xi = xi - mean(xi, 'omitnan');
            for j = 1:K
                xj = Xn(j, :).';
                xj = xj - mean(xj, 'omitnan');

                r = xcorr(xi, xj, maxLag, 'coeff');
                r = max(min(r, clipR), -clipR);

                if doFisherZ
                    v = atanh(r);
                else
                    v = r;
                end

                finiteMask = isfinite(v);
                Acc(i,j,finiteMask)  = Acc(i,j,finiteMask)  + reshape(v(finiteMask), 1,1,[]);
                nEff(i,j,finiteMask) = nEff(i,j,finiteMask) + 1;
            end
        end
    end
end

XcorrMat = nan(K, K, L);
valid = nEff > 0;
XcorrMat(valid) = Acc(valid) ./ nEff(valid);

if doFisherZ
    XcorrMat = tanh(XcorrMat);
end
end
