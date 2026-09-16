function S = globalDA_motifH_perTrial_xcorr_func(tbytDat_hAligned, tbytDat_DAglobalAligned, trI, varargin)
%GLOBALDA_MOTIFH_PERTRIAL_XCORR_FUNC
%   Per-trial cross-correlogram between the GLOBAL cortical DA transient
%   and each of K calcium motifs' H, built to match the motif-motif
%   pipeline (motifH_perTrial_xcorr_residual_posLag_withinTrialTimeShuff_func)
%   step for step, so the two results are directly comparable:
%
%     1. Both signals are windowed onto one common time grid with
%        stack_trials_H (100 ms bins, 50 ms step by default) and z-scored
%        per signal using global stats across all trials and bins. Running
%        DA through the SAME function as H is what guarantees identical
%        binning and normalization.
%     2. Trials are split by type (Hit, CR). Optionally, each signal's own
%        trial-type PSTH is subtracted from every trial of that type
%        (doPSTHSubtraction = true) so the correlogram reflects residual
%        within-trial temporal covariation rather than shared task locking.
%     3. Per trial, per motif: xcorr(da, h_k, maxLag, 'coeff'), clipped,
%        Fisher-z, averaged across trials, tanh back. Result: [K x L].
%     4. Two independent nulls, as in the motif-motif function:
%          within-trial circshift  -> destroys lag-specific alignment
%          across-trial permutation -> destroys trial pairing
%
%   DIRECTION: positive lag = H LEADS DA.
%     MATLAB's xcorr(x, y) at lag +k correlates x(t+k) with y(t), i.e. y
%     precedes x. With x = DA and y = H, a peak at positive lag means motif
%     activity precedes the DA transient. Negative lag = DA leads H.
%     (This matches the convention in the original chunk-level script.)
%
%   NO POSITIVE-LAG POOLING: the full correlogram is returned at 50 ms lag
%   resolution over +/- maxLagSec. Pool downstream if you need a scalar.
%
%   S = globalDA_motifH_perTrial_xcorr_func(tbytDat_hAligned, tbytDat_DAglobalAligned, trI, ...)
%
% INPUTS
%   tbytDat_hAligned        : {2 x N} cell, row 1 = [K x Tn] H, row 2 = [1 x Tn] time
%   tbytDat_DAglobalAligned : {2 x N} cell, row 1 = [1 x Tn] DA (dF/F), row 2 = time
%   trI                     : struct with logical hitI, crI (shared across both)
%
% NAME-VALUE
%   'doPSTHSubtraction' : true (default) | false
%   'maxLagSec'         : 2.0   lag window, each side
%   'Epoch'             : []    -> stack_trials_H default (full available range)
%   'Win', 'Step'       : 0.10, 0.05  stack_trials_H binning (Step sets the lag resolution)
%   'zscoreHs'          : true  global per-signal z-score inside stack_trials_H
%   'minCorrectTrials'  : 10    a stream with fewer trials is skipped (NaN outputs)
%   'doFisherZ'         : true
%   'clipR'             : 0.999
%   'doTimeShuffle'     : true  within-trial circshift null
%   'doTrialShuffle'    : true  across-trial permutation null
%   'nShuffle'          : 1000
%   'shuffleMethod'     : 'circshift' | 'permute'  (within-trial null only)
%   'rngSeed'           : []
%   'useParfor'         : false
%   'verbose'           : true
%
% OUTPUT (S)
%   S.meta   : header info, trial counts, skipped flags, K, L, lag axis
%   S.params : every option above, plus stack_trials_H params, for provenance
%   S.obs    : per stream (hit, cr):
%       psthH_<stream>        [K x nW]  motif PSTH (always computed)
%       psthDA_<stream>       [1 x nW]  DA PSTH   (always computed)
%       XcorrMat_<stream>     [K x L]   observed correlogram
%   S.shuf   : per stream and per null (suffix '' = within-trial, '_trialShuffle'):
%       XcorrMat_<stream><null>        [K x L x nShuffle]  full null draws
%       mean_/std_/z_/p_<stream><null> [K x L]             per-lag stats
%   Field names do NOT change with doPSTHSubtraction; S.params records it.

%% -------------------- parse --------------------
p = inputParser;
p.addParameter('doPSTHSubtraction', true, @(x) islogical(x) && isscalar(x));
p.addParameter('maxLagSec', 2.0, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('Epoch', [], @(v) isempty(v) || (isnumeric(v) && numel(v) == 2));
p.addParameter('Win', 0.10, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('Step', 0.05, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('zscoreHs', true, @(x) islogical(x) && isscalar(x));
p.addParameter('minCorrectTrials', 10, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('doFisherZ', true, @(x) islogical(x) && isscalar(x));
p.addParameter('clipR', 0.999, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
p.addParameter('doTimeShuffle', true, @(x) islogical(x) && isscalar(x));
p.addParameter('doTrialShuffle', true, @(x) islogical(x) && isscalar(x));
p.addParameter('nShuffle', 1000, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('shuffleMethod', 'circshift', @(s) any(strcmpi(string(s), ["circshift","permute"])));
p.addParameter('rngSeed', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('useParfor', false, @(x) islogical(x) && isscalar(x));
p.addParameter('verbose', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;
shuffleMethod = lower(string(opt.shuffleMethod));

if ~isempty(opt.rngSeed), rng(opt.rngSeed); end

assert(isfield(trI, 'hitI') && isfield(trI, 'crI'), 'trI must contain hitI and crI.');
N = size(tbytDat_hAligned, 2);
assert(size(tbytDat_DAglobalAligned, 2) == N, ...
    'H has %d trials but DA has %d -- trial indexing must match.', N, size(tbytDat_DAglobalAligned, 2));

%% -------------------- common grid: run BOTH through stack_trials_H --------------------
% Epoch must be identical for both calls or the grids won't line up. If
% the user didn't pin one, derive it from the INTERSECTION of the two
% signals' time ranges, so neither is extrapolated.
if isempty(opt.Epoch)
    [tH0, tH1]   = local_timeRange(tbytDat_hAligned);
    [tDA0, tDA1] = local_timeRange(tbytDat_DAglobalAligned);
    Epoch = [max(tH0, tDA0), min(tH1, tDA1)];
    assert(Epoch(2) > Epoch(1), 'H and DA time ranges do not overlap.');
else
    Epoch = opt.Epoch;
end

Hs  = stack_trials_H(tbytDat_hAligned,        'Epoch', Epoch, 'Win', opt.Win, 'Step', opt.Step, 'zscore', opt.zscoreHs);
DAs = stack_trials_H(tbytDat_DAglobalAligned, 'Epoch', Epoch, 'Win', opt.Win, 'Step', opt.Step, 'zscore', opt.zscoreHs, 'motifIdx', 1);

assert(isequal(Hs.winCtrs, DAs.winCtrs), 'stack_trials_H produced different grids for H and DA -- check Epoch/Win/Step.');

[~, K, nW] = size(Hs.Y3);            % Hs.Y3 : N x K x nW
DA3 = DAs.Y3;                         % N x 1 x nW
stepSec = opt.Step;

maxLag = round(opt.maxLagSec / stepSec);
lags   = -maxLag:maxLag;
lagSec = lags * stepSec;
L      = numel(lags);
assert(maxLag < nW, 'maxLagSec (%.2f s) exceeds the trial window (%d bins x %.3f s).', opt.maxLagSec, nW, stepSec);

if opt.verbose
    fprintf('Grid: %d bins x %.3f s (Epoch [%.2f %.2f]); K = %d motifs; lags +/-%d bins (%.2f s), L = %d.\n', ...
        nW, stepSec, Epoch(1), Epoch(2), K, maxLag, opt.maxLagSec, L);
end

%% -------------------- streams --------------------
streams = {'hit', 'cr'};
masks   = {logical(trI.hitI(:)'), logical(trI.crI(:)')};

S = struct();
S.meta = struct('K', K, 'nW', nW, 'L', L, 'lags', lags, 'lagSec', lagSec, ...
    'stepSec', stepSec, 'Epoch', Epoch, 'winCtrs', Hs.winCtrs, 'nTrialsTotal', N);
S.params = opt;
S.params.stackParams = Hs.params;
S.params.direction   = 'positive lag = H leads DA';
S.obs  = struct();
S.shuf = struct();

for si = 1:numel(streams)
    st = streams{si};
    mask = masks{si};
    assert(numel(mask) == N, 'trI.%sI has %d entries, expected %d.', st, numel(mask), N);

    nTr = sum(mask);
    S.meta.(['nTrials_' st]) = nTr;
    skipped = nTr < opt.minCorrectTrials;
    S.meta.(['skipped_' st]) = skipped;

    if skipped
        if opt.verbose
            warning('[%s] %d trials < minCorrectTrials=%d -- skipping stream.', st, nTr, opt.minCorrectTrials);
        end
        S.obs.(['psthH_' st])    = nan(K, nW);
        S.obs.(['psthDA_' st])   = nan(1, nW);
        S.obs.(['XcorrMat_' st]) = nan(K, L);
        continue;
    end

    Xh  = Hs.Y3(mask, :, :);      % nTr x K x nW
    Xda = DA3(mask, :, :);        % nTr x 1 x nW

    % ---- PSTH: always computed (for reference); subtracted only if asked ----
    psthH  = squeeze(mean(Xh,  1, 'omitnan'));   % K x nW
    psthDA = squeeze(mean(Xda, 1, 'omitnan'))';  % 1 x nW
    if K == 1, psthH = psthH(:)'; end

    if opt.doPSTHSubtraction
        Xh_used  = Xh  - reshape(psthH,  1, K, nW);
        Xda_used = Xda - reshape(psthDA, 1, 1, nW);
    else
        Xh_used  = Xh;
        Xda_used = Xda;
    end

    S.obs.(['psthH_' st])  = psthH;
    S.obs.(['psthDA_' st]) = psthDA;

    % ---- observed ----
    S.obs.(['XcorrMat_' st]) = local_xcorr_DA_vs_H(Xda_used, Xh_used, maxLag, opt.doFisherZ, opt.clipR);

    if opt.verbose
        fprintf('[%s] n=%d trials, PSTH subtraction %s -> [%d x %d] correlogram.\n', ...
            st, nTr, char(string(opt.doPSTHSubtraction)), K, L);
    end

    % ---- nulls ----
    if opt.doTimeShuffle
        nullT = local_nullLoop(Xda_used, Xh_used, maxLag, opt, @(da, h) local_withinTrialShuffle(da, h, shuffleMethod));
        S.shuf.(['XcorrMat_' st]) = nullT;
        [S.shuf.(['mean_' st]), S.shuf.(['std_' st]), S.shuf.(['z_' st]), S.shuf.(['p_' st])] = ...
            local_shufStats(S.obs.(['XcorrMat_' st]), nullT);
    end
    if opt.doTrialShuffle
        nullS = local_nullLoop(Xda_used, Xh_used, maxLag, opt, @(da, h) local_trialShuffle(da, h));
        S.shuf.(['XcorrMat_' st '_trialShuffle']) = nullS;
        [S.shuf.(['mean_' st '_trialShuffle']), S.shuf.(['std_' st '_trialShuffle']), ...
         S.shuf.(['z_' st '_trialShuffle']),    S.shuf.(['p_' st '_trialShuffle'])] = ...
            local_shufStats(S.obs.(['XcorrMat_' st]), nullS);
    end
end
end

%% ========================================================================
%  local functions
%% ========================================================================
function [t0, t1] = local_timeRange(alignedCell)
t0 = inf; t1 = -inf;
for n = 1:size(alignedCell, 2)
    t = alignedCell{2, n};
    if isempty(t), continue; end
    t0 = min(t0, min(t)); t1 = max(t1, max(t));
end
end


function XcorrMat = local_xcorr_DA_vs_H(Xda, Xh, maxLag, doFisherZ, clipR)
% Per-trial xcorr(da, h_k), Fisher-z averaged across trials.
%   Xda : nTr x 1 x nW ;  Xh : nTr x K x nW
% Mirrors computeMotifXcorr_perTrial_FisherZ: NaN bins -> 0 after
% mean-centering, clip, atanh, accumulate over finite lags, tanh back.
[nTr, K, ~] = size(Xh);
L = 2*maxLag + 1;
Acc  = zeros(K, L);
nEff = zeros(K, L);

for n = 1:nTr
    da = squeeze(Xda(n, 1, :));  da(~isfinite(da)) = 0;  da = da - mean(da);
    for k = 1:K
        h = squeeze(Xh(n, k, :)); h(~isfinite(h)) = 0;   h = h - mean(h);
        r = xcorr(da, h, maxLag, 'coeff');       % +lag: h precedes da
        r = max(min(r, clipR), -clipR);
        if doFisherZ, v = atanh(r); else, v = r; end
        ok = isfinite(v);
        Acc(k, ok)  = Acc(k, ok)  + v(ok)';
        nEff(k, ok) = nEff(k, ok) + 1;
    end
end

XcorrMat = nan(K, L);
valid = nEff > 0;
XcorrMat(valid) = Acc(valid) ./ nEff(valid);
if doFisherZ, XcorrMat = tanh(XcorrMat); end
end


function nullStack = local_nullLoop(Xda, Xh, maxLag, opt, shuffleFn)
% Runs the shuffle nShuffle times. shuffleFn(da, h) returns shuffled copies.
K = size(Xh, 2); L = 2*maxLag + 1;
nullStack = nan(K, L, opt.nShuffle);
doF = opt.doFisherZ; cR = opt.clipR; nS = opt.nShuffle;

if opt.useParfor
    parfor s = 1:nS
        [daS, hS] = shuffleFn(Xda, Xh);
        nullStack(:, :, s) = local_xcorr_DA_vs_H(daS, hS, maxLag, doF, cR);
    end
else
    for s = 1:nS
        [daS, hS] = shuffleFn(Xda, Xh);
        nullStack(:, :, s) = local_xcorr_DA_vs_H(daS, hS, maxLag, doF, cR);
        if opt.verbose && mod(s, 200) == 0, fprintf('    shuffle %d/%d\n', s, nS); end
    end
end
end


function [daS, hS] = local_withinTrialShuffle(Xda, Xh, method)
% Within-trial shift/permute of EACH signal independently, per trial --
% same as withinTrialShuffle_independent in the motif-motif function.
% Destroys lag-specific alignment; preserves trial pairing.
[nTr, K, nW] = size(Xh);
daS = Xda; hS = Xh;
switch method
    case "circshift"
        for n = 1:nTr
            daS(n, 1, :) = circshift(Xda(n, 1, :), randi(nW) - 1, 3);
            for k = 1:K
                hS(n, k, :) = circshift(Xh(n, k, :), randi(nW) - 1, 3);
            end
        end
    case "permute"
        for n = 1:nTr
            daS(n, 1, :) = Xda(n, 1, randperm(nW));
            for k = 1:K
                hS(n, k, :) = Xh(n, k, randperm(nW));
            end
        end
end
end


function [daS, hS] = local_trialShuffle(Xda, Xh)
% Permute WHICH TRIAL the DA trace comes from, leaving each trace intact.
% Destroys trial pairing; preserves within-trial structure. Only one
% signal needs permuting to break the pairing.
nTr = size(Xda, 1);
daS = Xda(randperm(nTr), :, :);
hS  = Xh;
end


function [mu, sig, z, pv] = local_shufStats(obsMat, shufMat)
mu  = mean(shufMat, 3, 'omitnan');
sig = std(shufMat, 0, 3, 'omitnan');
sig(sig == 0) = eps;
z = (obsMat - mu) ./ sig;
nS = size(shufMat, 3);
devObs  = abs(obsMat - mu);
devShuf = abs(shufMat - mu);
pv = (1 + sum(devShuf >= devObs, 3, 'omitnan')) ./ (nS + 1);
end