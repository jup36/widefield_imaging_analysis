function S = motifH_perTrial_xcorr_residual_posLag_withinTrialTimeShuff_func(filePath, fileKeyword, varargin)
% MOTIFH_PERTRIAL_XCORR_RESIDUAL_POSLAG_WITHINTRIALTIMESHUFFLE_FUNC
%   Residual (PSTH-subtracted) version of
%   motifH_perTrial_crosscorr_posLag_withinTrialTimeShuffle_func.
%
%   Purpose: the original function computes per-trial motif-motif xcorr on
%   RAW single-trial traces. Because every trial of a given type (Go/NoGo)
%   shares the same task-locked template (tone-locked rise, response-locked
%   rise, etc.), two motifs with different fixed latencies will show a
%   lagged correlation on every trial purely from that shared template --
%   a "signal correlation" -- even with zero real trial-to-trial coupling.
%
%   This function instead:
%     1) Computes each motif's PSTH (across-trial mean waveform) SEPARATELY
%        for Hit, CR, FA, and Miss trials (each trial type's own PSTH is
%        computed and subtracted independently -- never pooling trial
%        types before this step, since that would blend outcome-driven
%        mean-level differences into the template).
%     2) Subtracts that PSTH from each individual trial's motif trace,
%        leaving a residual (trial-to-trial fluctuation around the
%        canonical shape).
%     3) Runs the same per-trial xcorr + Fisher-z averaging + within-trial
%        circshift/permute shuffle null on the RESIDUALS instead of raw
%        traces.
%
%   Correlation that survives in the residuals is much harder to explain as
%   "two independently-templated pathways with different latencies" and is
%   better evidence of real dynamic coupling.
%
%   Hit, CR, FA, and Miss trial types are ALWAYS PSTH-subtracted
%   separately (never pooled before residualization) -- this is what makes
%   the two pooled outputs below ("combined_correct" and "combinedAll")
%   methodologically sound rather than just an average-of-raw-traces: each
%   trial keeps its own condition-specific PSTH removed, and only the
%   RESULTING residuals get pooled, so a genuine mean-level difference
%   between trial types (real signal, not noise) never leaks into the
%   pooled xcorr as spurious "coupling."
%
%   Sessions/trial-types with fewer than 'minCorrectTrials' trials are
%   SKIPPED (not silently downgraded) -- flagged in S.meta.skipped_* and
%   logged to the console, consistent with the correct-trials pipeline
%   design principle already in use. This threshold is applied identically
%   to Hit, CR, FA, and Miss (same parameter, same style of check for all
%   four trial types).
%
%   OUTPUTS -- FOUR pooled/per-type xcorr streams are computed:
%     - Hit-only, CR-only          : as before, unchanged.
%     - combined_correct           : Hit+CR residuals pooled. Skipped
%                                     entirely (not partially) if EITHER
%                                     Hit or CR was individually skipped --
%                                     "combined correct" without one of its
%                                     two constituents isn't a meaningful
%                                     combined-correct estimate.
%     - combinedAll                : Hit+CR+FA+Miss residuals pooled, using
%                                     WHICHEVER of the four individually
%                                     cleared minCorrectTrials this session
%                                     (session-level PARTIAL inclusion is
%                                     allowed here, unlike combined_correct
%                                     -- e.g. a session with too few FA
%                                     trials still gets a combinedAll built
%                                     from Hit+CR+Miss). The exact set of
%                                     trial types actually included is
%                                     recorded in
%                                     S.meta.included_trialTypes_combinedAll,
%                                     since this composition can legitimately
%                                     vary session to session.
%   FA-only and Miss-only xcorr are deliberately NOT computed/saved as
%   standalone outputs -- their residuals exist only internally, as inputs
%   to combinedAll.
%
%   S = motifH_perTrial_xcorr_residual_posLag_withinTrialTimeShuffle_func(filePath, fileKeyword, ...)
%
% NAME-VALUE PAIRS (explicit)
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
%   'minCorrectTrials'      : minimum trial count required, per trial type
%                             (Hit, CR, FA, AND Miss -- same threshold,
%                             same check, applied to all four), to compute
%                             a PSTH + residual for that type in this
%                             session (default 10). Below this, that trial
%                             type is SKIPPED (NaN'd out / excluded from
%                             pooling, and logged), not computed on a thin
%                             PSTH.
%   'showProgress'          : true (default). Prints progress during shuffles (works w/ parfor).
%   'progressEvery'         : positive integer, default 25. Print every N shuffles.
%
% OUTPUT (S)
%   S.meta   : header, minCorrectTrials, nTrials_hit, nTrials_cr,
%              nTrials_fa, nTrials_miss, skipped_hit, skipped_cr,
%              skipped_fa, skipped_miss, nTrials_combinedCorrect,
%              skipped_combinedCorrect, nTrials_combinedAll,
%              skipped_combinedAll, included_trialTypes_combinedAll
%   S.params : same param bookkeeping as the raw-trace version
%   S.obs    : psth_hit, psth_cr  [K x T] canonical waveforms that were
%              subtracted (FA/Miss PSTHs are computed internally but not
%              saved -- see note above);
%              XcorrMat_sess_hit_residual, XcorrMat_sess_cr_residual,
%              XcorrMat_sess_combinedCorrect_residual,
%              XcorrMat_sess_combinedAll_residual [K x K x L] full-lag-curve
%              xcorr on residuals;
%              xcorrPosLagMat_hit_residual, xcorrPosLagMat_cr_residual,
%              xcorrPosLagMat_combinedCorrect_residual,
%              xcorrPosLagMat_combinedAll_residual [K x K] pooled
%              positive-lag summary on residuals
%   S.shuf   : shuffle null distributions + stats (mean/std/z/p), mirroring
%              the raw-trace version's fields but suffixed _hit_residual /
%              _cr_residual / _combinedCorrect_residual /
%              _combinedAll_residual, plus per-lag null curves in S.shuf.curve

%% -------------------- Parse inputs --------------------
p = inputParser;

p.addParameter('xcorrLagWindowSec',   1.0, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('posLagPoolWindowSec', 0.5, @(x) isnumeric(x) && isscalar(x) && x>0);

p.addParameter('doTimeShuffle', true, @(x) islogical(x) && isscalar(x));
p.addParameter('nShuffle', 1000, @(x) isnumeric(x) && isscalar(x) && x>=1 && x==round(x));
p.addParameter('shuffleMethod', 'circshift', @(s) ischar(s) || isstring(s));
p.addParameter('rngSeed', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('zscoreHs', true, @(x) islogical(x) && isscalar(x));
p.addParameter('useSymmetry', true, @(x) islogical(x) && isscalar(x));
p.addParameter('doFisherZ', true, @(x) islogical(x) && isscalar(x));
p.addParameter('clipR', 0.999, @(x) isnumeric(x) && isscalar(x) && x>0 && x<1);

% Minimum trial count required to compute a PSTH/residual for a given
% trial type. Below this, SKIP (do not compute on a thin/unstable PSTH --
% an under-sampled PSTH under-subtracts the template and leaves signal-
% correlation residue behind, which would silently look like the
% "coupling survived" result you're trying not to fool yourself with).
% Applied identically to Hit, CR, FA, and Miss.
p.addParameter('minCorrectTrials', 10, @(x) isnumeric(x) && isscalar(x) && x>=1 && x==round(x));

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

minCorrectTrials    = p.Results.minCorrectTrials;

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
filePath_H        = cell2mat(find_keyword_containing_files(filePath_matfiles, fileKeyword, 'recursive', true));
filePath_B        = cell2mat(GrabFiles_sort_trials([header '*' keyword_beh], 0, {filePath_matfiles}));

load(filePath_H, 'tbytDat_hAligned');
load(filePath_B, 'tbytDat');

trI = trialTypeInfoAuditoryGngTbytDat(tbytDat);
assert(size(tbytDat_hAligned, 2) == numel(tbytDat), 'Mismatch in trial counts.');

%% ASSUMPTION: trI exposes logical fields 'hitI', 'crI', 'faI', and 'missI'
%% (matching the trIdC convention -- hitI/missI/crI/faI/goI/nogoI -- used
%% elsewhere in this pipeline). If trialTypeInfoAuditoryGngTbytDat names
%% these differently, update the lines below accordingly -- everything
%% downstream is agnostic to the field name.
if ~isfield(trI, 'hitI') || ~isfield(trI, 'crI')
    error(['trI from trialTypeInfoAuditoryGngTbytDat is missing ''hitI'' and/or ''crI''. ' ...
           'This function requires correct-trial-only masks for Go (Hit) and NoGo (CR) trials.']);
end
if ~isfield(trI, 'faI') || ~isfield(trI, 'missI')
    error(['trI from trialTypeInfoAuditoryGngTbytDat is missing ''faI'' and/or ''missI''. ' ...
           'These are required (internally) to build the combinedAll (correct+incorrect) pool.']);
end

%% -------------------- 1) Stack trials --------------------
Hs = stack_trials_H(tbytDat_hAligned, 'zscore', doZscore);
[~, K, T] = size(Hs.Y3);

stepSec = Hs.params.Step;

% Full xcorr lag window (curve)
nBin_total = round(xcorrLagWindowSec / stepSec);     % seconds -> bins
lags = -nBin_total:nBin_total;
L = numel(lags);

% Positive-lag pooling window (summary)
nBin_pool = round(posLagPoolWindowSec / stepSec);    % seconds -> bins
posLagMask_pool = (lags > 0 & lags <= nBin_pool);

%% -------------------- 2) Build per-trial-type residuals (Hit, CR) --------------------
nTrials_hit = sum(trI.hitI);
nTrials_cr  = sum(trI.crI);

skipped_hit = nTrials_hit < minCorrectTrials;
skipped_cr  = nTrials_cr  < minCorrectTrials;

if skipped_hit
    warning('[%s] Hit trial count (%d) < minCorrectTrials (%d). Skipping Hit residual xcorr for this session.', ...
        header, nTrials_hit, minCorrectTrials);
    psth_hit    = nan(K, T);
    X_hit_resid = [];
else
    X_hit = Hs.Y3(trI.hitI, :, :);
    [psth_hit, X_hit_resid] = computeResidual_perMotif(X_hit);
end

if skipped_cr
    warning('[%s] CR trial count (%d) < minCorrectTrials (%d). Skipping CR residual xcorr for this session.', ...
        header, nTrials_cr, minCorrectTrials);
    psth_cr    = nan(K, T);
    X_cr_resid = [];
else
    X_cr = Hs.Y3(trI.crI, :, :);
    [psth_cr, X_cr_resid] = computeResidual_perMotif(X_cr);
end

%% -------------------- 2b) FA / Miss residuals (INTERNAL ONLY) --------------------
% Same minCorrectTrials threshold, same skip-not-downgrade check as
% Hit/CR above. FA/Miss PSTHs and xcorr are NOT saved as standalone
% outputs -- these residuals exist only to feed combinedAll below (per
% project decision: FA/Miss get no separate FA-only/Miss-only results).
nTrials_fa   = sum(trI.faI);
nTrials_miss = sum(trI.missI);

skipped_fa   = nTrials_fa   < minCorrectTrials;
skipped_miss = nTrials_miss < minCorrectTrials;

if skipped_fa
    warning('[%s] FA trial count (%d) < minCorrectTrials (%d). Excluding FA from combinedAll for this session.', ...
        header, nTrials_fa, minCorrectTrials);
    X_fa_resid = [];
else
    X_fa = Hs.Y3(trI.faI, :, :);
    [~, X_fa_resid] = computeResidual_perMotif(X_fa);   % psth_fa discarded -- internal use only
end

if skipped_miss
    warning('[%s] Miss trial count (%d) < minCorrectTrials (%d). Excluding Miss from combinedAll for this session.', ...
        header, nTrials_miss, minCorrectTrials);
    X_miss_resid = [];
else
    X_miss = Hs.Y3(trI.missI, :, :);
    [~, X_miss_resid] = computeResidual_perMotif(X_miss);   % psth_miss discarded -- internal use only
end

%% -------------------- 2c) Combined pools --------------------
% Option A: pool ALREADY-residualized per-trial-type traces -- each trial
% keeps its own condition-specific PSTH subtracted; only the RESULTING
% residuals get pooled, so a genuine Hit-vs-CR-vs-FA-vs-Miss mean-level
% difference never leaks into the pooled xcorr as spurious "coupling."

% combined_correct: Hit + CR only. Skipped entirely (not partially) if
% EITHER was skipped -- "combined correct" without one of its two
% constituents isn't a meaningful combined-correct estimate.
skipped_combinedCorrect = skipped_hit || skipped_cr;
if skipped_combinedCorrect
    X_combinedCorrect_resid = [];
else
    X_combinedCorrect_resid = cat(1, X_hit_resid, X_cr_resid);
end

% combinedAll: whichever of Hit/CR/FA/Miss individually cleared
% minCorrectTrials, pooled together (session-level PARTIAL inclusion is
% allowed here, unlike combined_correct above). included_trialTypes_combinedAll
% records exactly which types actually went in, since this composition can
% legitimately vary session to session.
includedTypes_combinedAll = {};
X_combinedAll_parts = {};
if ~skipped_hit,  X_combinedAll_parts{end+1} = X_hit_resid;  includedTypes_combinedAll{end+1} = 'hit';  end %#ok<AGROW>
if ~skipped_cr,   X_combinedAll_parts{end+1} = X_cr_resid;   includedTypes_combinedAll{end+1} = 'cr';   end %#ok<AGROW>
if ~skipped_fa,   X_combinedAll_parts{end+1} = X_fa_resid;   includedTypes_combinedAll{end+1} = 'fa';   end %#ok<AGROW>
if ~skipped_miss, X_combinedAll_parts{end+1} = X_miss_resid; includedTypes_combinedAll{end+1} = 'miss'; end %#ok<AGROW>

skipped_combinedAll = isempty(X_combinedAll_parts);
if skipped_combinedAll
    warning('[%s] All four trial types (Hit/CR/FA/Miss) skipped (below minCorrectTrials) -- combinedAll unavailable for this session.', header);
    X_combinedAll_resid = [];
else
    X_combinedAll_resid = cat(1, X_combinedAll_parts{:});
    if numel(includedTypes_combinedAll) < 2
        warning('[%s] combinedAll includes only 1 trial type (%s) this session -- effectively identical to that type alone, not a true combination.', ...
            header, includedTypes_combinedAll{1});
    end
end

%% -------------------- 3) Observed per-trial averaged xcorr on residuals --------------------
if ~isempty(X_hit_resid)
    XcorrMat_sess_hit_residual = computeMotifXcorr_perTrial_FisherZ(X_hit_resid, nBin_total, useSymmetry, doFisherZ, clipR);
else
    XcorrMat_sess_hit_residual = nan(K, K, L);
end

if ~isempty(X_cr_resid)
    XcorrMat_sess_cr_residual = computeMotifXcorr_perTrial_FisherZ(X_cr_resid, nBin_total, useSymmetry, doFisherZ, clipR);
else
    XcorrMat_sess_cr_residual = nan(K, K, L);
end

if ~isempty(X_combinedCorrect_resid)
    XcorrMat_sess_combinedCorrect_residual = computeMotifXcorr_perTrial_FisherZ(X_combinedCorrect_resid, nBin_total, useSymmetry, doFisherZ, clipR);
else
    XcorrMat_sess_combinedCorrect_residual = nan(K, K, L);
end

if ~isempty(X_combinedAll_resid)
    XcorrMat_sess_combinedAll_residual = computeMotifXcorr_perTrial_FisherZ(X_combinedAll_resid, nBin_total, useSymmetry, doFisherZ, clipR);
else
    XcorrMat_sess_combinedAll_residual = nan(K, K, L);
end

% Pooled pos-lag summary uses ONLY 0-posLagPoolWindowSec
xcorrPosLagMat_hit_residual = squeeze(mean(XcorrMat_sess_hit_residual(:, :, posLagMask_pool), 3, 'omitnan'));
xcorrPosLagMat_cr_residual  = squeeze(mean(XcorrMat_sess_cr_residual(:, :,  posLagMask_pool), 3, 'omitnan'));

xcorrPosLagMat_combinedCorrect_residual = squeeze(mean(XcorrMat_sess_combinedCorrect_residual(:, :, posLagMask_pool), 3, 'omitnan'));
xcorrPosLagMat_combinedAll_residual     = squeeze(mean(XcorrMat_sess_combinedAll_residual(:, :,     posLagMask_pool), 3, 'omitnan'));

%% -------------------- 4) Package observed outputs --------------------
S = struct();
S.meta = struct( ...
    'header',           header, ...
    'minCorrectTrials', minCorrectTrials, ...
    'nTrials_hit',       nTrials_hit, ...
    'nTrials_cr',        nTrials_cr, ...
    'nTrials_fa',        nTrials_fa, ...
    'nTrials_miss',      nTrials_miss, ...
    'skipped_hit',       skipped_hit, ...
    'skipped_cr',        skipped_cr, ...
    'skipped_fa',        skipped_fa, ...
    'skipped_miss',      skipped_miss);

if skipped_combinedCorrect
    S.meta.nTrials_combinedCorrect = NaN;
else
    S.meta.nTrials_combinedCorrect = nTrials_hit + nTrials_cr;
end
S.meta.skipped_combinedCorrect = skipped_combinedCorrect;

if skipped_combinedAll
    S.meta.nTrials_combinedAll = NaN;
else
    S.meta.nTrials_combinedAll = size(X_combinedAll_resid, 1);
end
S.meta.skipped_combinedAll             = skipped_combinedAll;
S.meta.included_trialTypes_combinedAll = includedTypes_combinedAll;

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
S.params.nShuffle               = nShuffle;
S.params.rngSeed                = rngSeed;
S.params.zscoreHs               = doZscore;
S.params.useSymmetry            = useSymmetry;
S.params.doFisherZ              = doFisherZ;
S.params.clipR                  = clipR;
S.params.minCorrectTrials       = minCorrectTrials;

S.params.showProgress           = showProgress;
S.params.progressEvery          = progressEvery;

S.obs = struct();
S.obs.psth_hit = psth_hit;   % [K x T] canonical Hit waveform subtracted per motif
S.obs.psth_cr  = psth_cr;    % [K x T] canonical CR waveform subtracted per motif
% NOTE: psth_fa / psth_miss deliberately not saved -- internal-only, see
% header note and S.meta.included_trialTypes_combinedAll for provenance.

S.obs.XcorrMat_sess_hit_residual = XcorrMat_sess_hit_residual;
S.obs.XcorrMat_sess_cr_residual  = XcorrMat_sess_cr_residual;
S.obs.XcorrMat_sess_combinedCorrect_residual = XcorrMat_sess_combinedCorrect_residual;
S.obs.XcorrMat_sess_combinedAll_residual     = XcorrMat_sess_combinedAll_residual;

S.obs.xcorrPosLagMat_hit_residual = xcorrPosLagMat_hit_residual;
S.obs.xcorrPosLagMat_cr_residual  = xcorrPosLagMat_cr_residual;
S.obs.xcorrPosLagMat_combinedCorrect_residual = xcorrPosLagMat_combinedCorrect_residual;
S.obs.xcorrPosLagMat_combinedAll_residual     = xcorrPosLagMat_combinedAll_residual;

%% -------------------- 5) Optional shuffle null + stats --------------------
if ~doTimeShuffle
    return;
end

% pooled samples for p-values
shuf_hit = nan(K, K, nShuffle);
shuf_cr  = nan(K, K, nShuffle);
shuf_combinedCorrect = nan(K, K, nShuffle);
shuf_combinedAll     = nan(K, K, nShuffle);

% per-lag descriptives (mean/std) via sums
sum_hit  = zeros(K, K, L);   sumsq_hit = zeros(K, K, L);   cnt_hit = zeros(K, K, L);
sum_cr   = zeros(K, K, L);   sumsq_cr  = zeros(K, K, L);   cnt_cr  = zeros(K, K, L);
sum_combinedCorrect = zeros(K, K, L);   sumsq_combinedCorrect = zeros(K, K, L);   cnt_combinedCorrect = zeros(K, K, L);
sum_combinedAll     = zeros(K, K, L);   sumsq_combinedAll     = zeros(K, K, L);   cnt_combinedAll     = zeros(K, K, L);

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
    % ---------------- HIT (residual) ----------------
    if ~isempty(X_hit_resid)
        Xs_hit = withinTrialShuffle_independent(X_hit_resid, shuffleMethod);
        Xc_hit = computeMotifXcorr_perTrial_FisherZ(Xs_hit, nBin_total, useSymmetry, doFisherZ, clipR);

        shuf_hit(:, :, s) = squeeze(mean(Xc_hit(:, :, posLagMask_pool), 3, 'omitnan'));

        m = isfinite(Xc_hit);
        Xc0 = Xc_hit; Xc0(~m) = 0;
        sum_hit   = sum_hit   + Xc0;
        sumsq_hit = sumsq_hit + Xc0.^2;
        cnt_hit   = cnt_hit   + double(m);
    end

    % ---------------- CR (residual) ----------------
    if ~isempty(X_cr_resid)
        Xs_cr = withinTrialShuffle_independent(X_cr_resid, shuffleMethod);
        Xc_cr = computeMotifXcorr_perTrial_FisherZ(Xs_cr, nBin_total, useSymmetry, doFisherZ, clipR);

        shuf_cr(:, :, s) = squeeze(mean(Xc_cr(:, :, posLagMask_pool), 3, 'omitnan'));

        m = isfinite(Xc_cr);
        Xc0 = Xc_cr; Xc0(~m) = 0;
        sum_cr   = sum_cr   + Xc0;
        sumsq_cr = sumsq_cr + Xc0.^2;
        cnt_cr   = cnt_cr   + double(m);
    end

    % ---------------- combined_correct (residual) ----------------
    % Shuffling the already-pooled residual pool is equivalent to
    % shuffling each trial independently within its own type first and
    % then pooling -- withinTrialShuffle_independent only cares about the
    % trial dimension, not which type each trial came from, so no
    % special-casing is needed here.
    if ~isempty(X_combinedCorrect_resid)
        Xs_cc = withinTrialShuffle_independent(X_combinedCorrect_resid, shuffleMethod);
        Xc_cc = computeMotifXcorr_perTrial_FisherZ(Xs_cc, nBin_total, useSymmetry, doFisherZ, clipR);

        shuf_combinedCorrect(:, :, s) = squeeze(mean(Xc_cc(:, :, posLagMask_pool), 3, 'omitnan'));

        m = isfinite(Xc_cc);
        Xc0 = Xc_cc; Xc0(~m) = 0;
        sum_combinedCorrect   = sum_combinedCorrect   + Xc0;
        sumsq_combinedCorrect = sumsq_combinedCorrect + Xc0.^2;
        cnt_combinedCorrect   = cnt_combinedCorrect   + double(m);
    end

    % ---------------- combinedAll (residual) ----------------
    if ~isempty(X_combinedAll_resid)
        Xs_ca = withinTrialShuffle_independent(X_combinedAll_resid, shuffleMethod);
        Xc_ca = computeMotifXcorr_perTrial_FisherZ(Xs_ca, nBin_total, useSymmetry, doFisherZ, clipR);

        shuf_combinedAll(:, :, s) = squeeze(mean(Xc_ca(:, :, posLagMask_pool), 3, 'omitnan'));

        m = isfinite(Xc_ca);
        Xc0 = Xc_ca; Xc0(~m) = 0;
        sum_combinedAll   = sum_combinedAll   + Xc0;
        sumsq_combinedAll = sumsq_combinedAll + Xc0.^2;
        cnt_combinedAll   = cnt_combinedAll   + double(m);
    end

    % Progress ping
    if showProgress
        send(dq, s);
    end
end

% Finalize per-lag null mean/std
[mu_hit, sd_hit] = finalizeMeanStd(sum_hit, sumsq_hit, cnt_hit);
[mu_cr,  sd_cr]  = finalizeMeanStd(sum_cr,  sumsq_cr,  cnt_cr);
[mu_combinedCorrect, sd_combinedCorrect] = finalizeMeanStd(sum_combinedCorrect, sumsq_combinedCorrect, cnt_combinedCorrect);
[mu_combinedAll,     sd_combinedAll]     = finalizeMeanStd(sum_combinedAll,     sumsq_combinedAll,     cnt_combinedAll);

S.shuf = struct();

S.shuf.xcorrPosLagMat_hit_residual = shuf_hit;
S.shuf.xcorrPosLagMat_cr_residual  = shuf_cr;
S.shuf.xcorrPosLagMat_combinedCorrect_residual = shuf_combinedCorrect;
S.shuf.xcorrPosLagMat_combinedAll_residual     = shuf_combinedAll;

[S.shuf.mean_hit_residual, S.shuf.std_hit_residual, S.shuf.z_hit_residual, S.shuf.p_hit_residual] = ...
    shufStats(xcorrPosLagMat_hit_residual, shuf_hit);
[S.shuf.mean_cr_residual,  S.shuf.std_cr_residual,  S.shuf.z_cr_residual,  S.shuf.p_cr_residual]  = ...
    shufStats(xcorrPosLagMat_cr_residual,  shuf_cr);
[S.shuf.mean_combinedCorrect_residual, S.shuf.std_combinedCorrect_residual, S.shuf.z_combinedCorrect_residual, S.shuf.p_combinedCorrect_residual] = ...
    shufStats(xcorrPosLagMat_combinedCorrect_residual, shuf_combinedCorrect);
[S.shuf.mean_combinedAll_residual, S.shuf.std_combinedAll_residual, S.shuf.z_combinedAll_residual, S.shuf.p_combinedAll_residual] = ...
    shufStats(xcorrPosLagMat_combinedAll_residual, shuf_combinedAll);

S.shuf.curve = struct();
S.shuf.curve.mean_hit_residual = mu_hit;
S.shuf.curve.std_hit_residual  = sd_hit;
S.shuf.curve.mean_cr_residual  = mu_cr;
S.shuf.curve.std_cr_residual   = sd_cr;
S.shuf.curve.mean_combinedCorrect_residual = mu_combinedCorrect;
S.shuf.curve.std_combinedCorrect_residual  = sd_combinedCorrect;
S.shuf.curve.mean_combinedAll_residual     = mu_combinedAll;
S.shuf.curve.std_combinedAll_residual      = sd_combinedAll;

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
function [psth, Xresid] = computeResidual_perMotif(X)
% COMPUTERESIDUAL_PERMOTIF
%   X      : [N x K x T] single-trial motif traces for ONE trial type
%             (e.g., Hit-only, CR-only, FA-only, or Miss-only), already
%             restricted to that trial type by the caller.
%   psth   : [K x T] across-trial mean waveform per motif (the "template"
%             being removed).
%   Xresid : [N x K x T] X with psth subtracted from every trial
%             (broadcast across trials).
%
% NOTE: this is deliberately a plain across-trial mean, not median or a
% smoothed estimate -- matches the two-stage averaging convention already
% used elsewhere in this pipeline (mean within animal/session before any
% further averaging). Trial-count adequacy for a stable PSTH is enforced
% by the caller via 'minCorrectTrials' (this function does not itself
% check for that -- it assumes the caller already skipped it if too few).

psth = squeeze(mean(X, 1, 'omitnan'));  % [K x T]
if size(X, 2) == 1
    % squeeze can over-collapse a singleton K dimension; guard shape
    psth = reshape(psth, 1, []);
end

Xresid = X - reshape(psth, 1, size(psth,1), size(psth,2));
end

%% ========================================================================
function Xs = withinTrialShuffle_independent(X, method)
% Independent within-trial shuffle per trial x motif (NOT a shared shift).
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
%   Unchanged from the raw-trace version: the per-trial mean-subtraction
%   here (xi = xi - mean(xi)) is a local centering step for the xcorr
%   itself and is independent of, and does not conflict with, the
%   across-trial PSTH subtraction already applied upstream in
%   computeResidual_perMotif.

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
            xi = xi - mean(xi, 'omitnan'); % local mean-subtract per trial
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