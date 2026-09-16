function S = motifH_perTrial_xcorr_residual_posLag_withinTrialTimeShuff_func(filePath, fileKeyword, varargin)
% MOTIFH_PERTRIAL_XCORR_RESIDUAL_POSLAG_WITHINTRIALTIMESHUFFLE_FUNC
%   Per-trial motif-motif positive-lag xcorr, with an OPTIONAL PSTH-
%   subtraction step controlled by 'doPSTHSubtraction'. TWO PARALLEL null
%   distributions are now computed for every stream (controlled
%   independently by 'doTimeShuffle' and 'doTrialShuffle', both true by
%   default):
%     - WITHIN-TRIAL shuffle (circshift/permute, unchanged from before):
%       for each trial and motif independently, circularly shifts (or
%       permutes) that motif's OWN time series within that trial. The
%       trial-to-trial pairing between motifs is untouched -- this
%       destroys fine-timescale/lag-specific alignment while leaving any
%       SHARED TRIAL-LEVEL state (e.g. a trial being "high engagement"
%       for both motifs at once) completely intact. Null hypothesis:
%       "no genuine lag-specific coupling beyond whatever these two
%       motifs share via overall per-trial state."
%     - TRIAL shuffle (NEW): for each motif independently, randomly
%       permutes WHICH TRIAL its data comes from, leaving each trial's
%       own internal time series completely untouched. This is the
%       mirror image of the within-trial shuffle -- it destroys the
%       trial-to-trial pairing (and therefore any shared trial-level
%       state) while preserving each trace's own fine-timescale
%       structure. Null hypothesis: "no relationship at all between
%       these two motifs' trial-to-trial fluctuations, of any kind."
%   These are DIFFERENT, complementary tests, not a stricter-vs-looser
%   pair -- see the two functions withinTrialShuffle_independent and
%   trialShuffle_independent below for the exact mechanics. Output field
%   names for the trial-shuffle null are suffixed "_trialShuffle" on top
%   of the existing "_residual" suffix (e.g.
%   S.shuf.xcorrPosLagMat_hit_residual_trialShuffle), so nothing about
%   the existing within-trial-shuffle fields changes at all.
%
%   PSTH-SUBTRACTED ("residual") MODE -- doPSTHSubtraction = true (DEFAULT,
%   unchanged from all prior behavior/already-completed Scotty runs):
%     Each motif's PSTH (across-trial mean waveform) is computed
%     SEPARATELY for Hit, CR, FA, and Miss trials, and subtracted from
%     every individual trial's motif trace before xcorr is computed. This
%     isolates trial-to-trial coupling from the shared task-locked
%     template both motifs ride on regardless of any real interaction.
%
%   RAW-TRACE MODE -- doPSTHSubtraction = false:
%     Identical pipeline in every other respect -- the ONLY difference is
%     that the PSTH is computed (still saved, for reference/comparison)
%     but NOT subtracted; xcorr is computed on the raw per-trial motif
%     traces directly. Applies to ALL four trial types (Hit/CR/FA/Miss)
%     identically, so every pooled output below is equally well-defined in
%     either mode.
%
%   POOLED OUTPUTS -- in addition to Hit-only and CR-only, FOUR pooled
%   streams are computed:
%     - combined_correct : Hit+CR pooled. Skipped entirely (not partially)
%                           if EITHER Hit or CR was individually skipped.
%     - combinedAll       : Hit+CR+FA+Miss pooled, using WHICHEVER of the
%                           four individually cleared minCorrectTrials
%                           this session (session-level PARTIAL inclusion
%                           allowed here). Composition recorded in
%                           S.meta.included_trialTypes_combinedAll.
%     - allGo   (NEW)     : Hit+Miss pooled -- i.e. every trial where a Go
%                           cue was presented, correct or not. Skipped
%                           entirely (not partially) if EITHER Hit or Miss
%                           was individually skipped -- "all Go trials"
%                           without one of its two constituents isn't a
%                           meaningful all-Go estimate, same logic as
%                           combined_correct.
%     - allNoGo (NEW)     : CR+FA pooled -- every trial where a NoGo cue
%                           was presented, correct or not. Skipped
%                           entirely (not partially) if EITHER CR or FA
%                           was individually skipped.
%   In every pooled stream, each constituent trial type keeps its OWN
%   condition-specific PSTH subtracted (if doPSTHSubtraction=true) BEFORE
%   pooling -- a genuine mean-level difference between trial types never
%   leaks into the pooled xcorr as spurious "coupling."
%
%   FA-only and Miss-only xcorr are still deliberately NOT computed/saved
%   as standalone outputs -- their traces (X_fa_used/X_miss_used) are
%   computed once, internally, and reused as inputs to combinedAll, allGo,
%   and allNoGo alike.
%
%   OUTPUT FIELD NAMES ARE IDENTICAL REGARDLESS OF doPSTHSubtraction -- the
%   mode used is recorded in S.params.doPSTHSubtraction for provenance.
%
%   S = motifH_perTrial_xcorr_residual_posLag_withinTrialTimeShuffle_func(filePath, fileKeyword, ...)
%
% NAME-VALUE PAIRS (explicit)
%   'xcorrLagWindowSec'     : total lag window for xcorr curve in seconds (default 1.0)
%   'posLagPoolWindowSec'   : positive-lag pooling window in seconds (default 0.5)
%   'doTimeShuffle'         : true (default). Controls the WITHIN-TRIAL
%                             (circshift/permute) null -- see header.
%   'doTrialShuffle'        : true (default, NEW). Controls the TRIAL
%                             shuffle null -- see header. Independent of
%                             doTimeShuffle; either, both, or neither can
%                             be enabled. If BOTH are false, the shuffle
%                             section is skipped entirely (observed-only,
%                             same as before).
%   'nShuffle'              : integer, default 1000 (applies to whichever
%                             null(s) are enabled -- same shuffle count
%                             for both if both are on).
%   'shuffleMethod'         : 'circshift' (default) or 'permute'
%   'rngSeed'               : [] (default) or scalar seed
%   'zscoreHs'              : true (default)
%   'useSymmetry'           : true (default)
%   'doFisherZ'             : true (default)
%   'clipR'                 : scalar in (0,1), default 0.999
%   'minCorrectTrials'      : minimum trial count required, per trial type
%                             (Hit, CR, FA, AND Miss -- same threshold,
%                             same check, applied to all four), to compute
%                             a PSTH + xcorr for that type in this session
%                             (default 10). Below this, that trial type is
%                             SKIPPED (excluded/NaN'd out and logged), not
%                             computed on a thin PSTH.
%   'doPSTHSubtraction'     : true (DEFAULT, unchanged prior behavior).
%                             If false, xcorr is computed on RAW per-trial
%                             traces instead of PSTH-subtracted residuals.
%                             Applies identically to Hit/CR/FA/Miss and
%                             therefore to every pooled output.
%   'showProgress'          : true (default). Prints progress during shuffles (works w/ parfor).
%   'progressEvery'         : positive integer, default 25. Print every N shuffles.
%
% OUTPUT (S)
%   S.meta   : header, minCorrectTrials, nTrials_hit, nTrials_cr,
%              nTrials_fa, nTrials_miss, skipped_hit, skipped_cr,
%              skipped_fa, skipped_miss, nTrials_combinedCorrect,
%              skipped_combinedCorrect, nTrials_combinedAll,
%              skipped_combinedAll, included_trialTypes_combinedAll,
%              nTrials_allGo, skipped_allGo, nTrials_allNoGo, skipped_allNoGo
%   S.params : same param bookkeeping as before, PLUS doPSTHSubtraction
%   S.obs    : psth_hit, psth_cr  [K x T] canonical waveforms (ALWAYS
%              computed regardless of doPSTHSubtraction, for reference;
%              psth_fa/psth_miss computed internally but NOT saved);
%              XcorrMat_sess_hit_residual, XcorrMat_sess_cr_residual,
%              XcorrMat_sess_combinedCorrect_residual,
%              XcorrMat_sess_combinedAll_residual,
%              XcorrMat_sess_allGo_residual, XcorrMat_sess_allNoGo_residual
%              [K x K x L]; xcorrPosLagMat_hit_residual,
%              xcorrPosLagMat_cr_residual,
%              xcorrPosLagMat_combinedCorrect_residual,
%              xcorrPosLagMat_combinedAll_residual,
%              xcorrPosLagMat_allGo_residual, xcorrPosLagMat_allNoGo_residual
%              [K x K]. Field names UNCHANGED regardless of
%              doPSTHSubtraction mode -- check S.params.doPSTHSubtraction
%              for provenance.
%   S.shuf   : shuffle null distributions + stats, field names mirroring
%              S.obs (suffixed _hit_residual / _cr_residual /
%              _combinedCorrect_residual / _combinedAll_residual /
%              _allGo_residual / _allNoGo_residual), plus per-lag null
%              curves in S.shuf.curve.

%% -------------------- Parse inputs --------------------
p = inputParser;

p.addParameter('xcorrLagWindowSec',   1.0, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('posLagPoolWindowSec', 0.5, @(x) isnumeric(x) && isscalar(x) && x>0);

p.addParameter('doTimeShuffle', true, @(x) islogical(x) && isscalar(x));
% NEW: independent toggle for the trial-shuffle null (permutes WHICH
% trial each motif's data comes from, rather than shifting time within a
% trial -- see header note for the exact distinction). Defaults to true
% so both nulls are computed by default; set false to skip it and save
% roughly half the shuffle-loop compute if you only want the original
% within-trial-shuffle null.
p.addParameter('doTrialShuffle', true, @(x) islogical(x) && isscalar(x));
p.addParameter('nShuffle', 1000, @(x) isnumeric(x) && isscalar(x) && x>=1 && x==round(x));
p.addParameter('shuffleMethod', 'circshift', @(s) ischar(s) || isstring(s));
p.addParameter('rngSeed', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('zscoreHs', true, @(x) islogical(x) && isscalar(x));
p.addParameter('useSymmetry', true, @(x) islogical(x) && isscalar(x));
p.addParameter('doFisherZ', true, @(x) islogical(x) && isscalar(x));
p.addParameter('clipR', 0.999, @(x) isnumeric(x) && isscalar(x) && x>0 && x<1);

% Applied identically to Hit, CR, FA, and Miss.
p.addParameter('minCorrectTrials', 10, @(x) isnumeric(x) && isscalar(x) && x>=1 && x==round(x));

% Toggle for PSTH subtraction. Default true = prior behavior (backward
% compatible with every already-completed Scotty run and already-
% collected result under this function's name).
p.addParameter('doPSTHSubtraction', true, @(x) islogical(x) && isscalar(x));

p.addParameter('showProgress', true, @(x) islogical(x) && isscalar(x));
p.addParameter('progressEvery', 25, @(x) isnumeric(x) && isscalar(x) && x>=1 && x==round(x));

p.parse(varargin{:});

xcorrLagWindowSec   = p.Results.xcorrLagWindowSec;
posLagPoolWindowSec = p.Results.posLagPoolWindowSec;

doTimeShuffle       = p.Results.doTimeShuffle;
doTrialShuffle       = p.Results.doTrialShuffle;
nShuffle            = p.Results.nShuffle;
shuffleMethod       = lower(string(p.Results.shuffleMethod));
rngSeed             = p.Results.rngSeed;
doZscore            = p.Results.zscoreHs;
useSymmetry         = p.Results.useSymmetry;
doFisherZ           = p.Results.doFisherZ;
clipR               = p.Results.clipR;

minCorrectTrials    = p.Results.minCorrectTrials;
doPSTHSubtraction   = p.Results.doPSTHSubtraction;

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
%% elsewhere in this pipeline).
if ~isfield(trI, 'hitI') || ~isfield(trI, 'crI')
    error(['trI from trialTypeInfoAuditoryGngTbytDat is missing ''hitI'' and/or ''crI''. ' ...
           'This function requires correct-trial-only masks for Go (Hit) and NoGo (CR) trials.']);
end
if ~isfield(trI, 'faI') || ~isfield(trI, 'missI')
    error(['trI from trialTypeInfoAuditoryGngTbytDat is missing ''faI'' and/or ''missI''. ' ...
           'These are required (internally) to build the combinedAll/allGo/allNoGo pools.']);
end

%% -------------------- 1) Stack trials --------------------
Hs = stack_trials_H(tbytDat_hAligned, 'zscore', doZscore);
[~, K, T] = size(Hs.Y3);

stepSec = Hs.params.Step;

nBin_total = round(xcorrLagWindowSec / stepSec);
lags = -nBin_total:nBin_total;
L = numel(lags);

nBin_pool = round(posLagPoolWindowSec / stepSec);
posLagMask_pool = (lags > 0 & lags <= nBin_pool);

%% -------------------- 2) Build per-trial-type traces (Hit, CR) --------------------
% X_*_used is what actually feeds xcorr below -- residual or raw,
% depending on doPSTHSubtraction.
nTrials_hit = sum(trI.hitI);
nTrials_cr  = sum(trI.crI);

skipped_hit = nTrials_hit < minCorrectTrials;
skipped_cr  = nTrials_cr  < minCorrectTrials;

if skipped_hit
    warning('[%s] Hit trial count (%d) < minCorrectTrials (%d). Skipping Hit xcorr for this session.', ...
        header, nTrials_hit, minCorrectTrials);
    psth_hit   = nan(K, T);
    X_hit_used = [];
else
    X_hit = Hs.Y3(trI.hitI, :, :);
    [psth_hit, X_hit_resid] = computeResidual_perMotif(X_hit);   % PSTH always computed, for reference regardless of mode
    if doPSTHSubtraction
        X_hit_used = X_hit_resid;
    else
        X_hit_used = X_hit;   % RAW-TRACE MODE: use unsubtracted trial traces directly
    end
end

if skipped_cr
    warning('[%s] CR trial count (%d) < minCorrectTrials (%d). Skipping CR xcorr for this session.', ...
        header, nTrials_cr, minCorrectTrials);
    psth_cr   = nan(K, T);
    X_cr_used = [];
else
    X_cr = Hs.Y3(trI.crI, :, :);
    [psth_cr, X_cr_resid] = computeResidual_perMotif(X_cr);
    if doPSTHSubtraction
        X_cr_used = X_cr_resid;
    else
        X_cr_used = X_cr;
    end
end

%% -------------------- 2b) FA / Miss traces (INTERNAL ONLY) --------------------
% Same minCorrectTrials threshold, same skip logic, same doPSTHSubtraction
% behavior as Hit/CR above. Neither PSTH nor xcorr for FA/Miss is saved as
% a standalone output -- these traces exist only to feed combinedAll,
% allGo, and allNoGo below.
nTrials_fa   = sum(trI.faI);
nTrials_miss = sum(trI.missI);

skipped_fa   = nTrials_fa   < minCorrectTrials;
skipped_miss = nTrials_miss < minCorrectTrials;

if skipped_fa
    warning('[%s] FA trial count (%d) < minCorrectTrials (%d). Excluding FA from combinedAll/allNoGo for this session.', ...
        header, nTrials_fa, minCorrectTrials);
    X_fa_used = [];
else
    X_fa = Hs.Y3(trI.faI, :, :);
    [~, X_fa_resid] = computeResidual_perMotif(X_fa);   % psth_fa discarded -- internal use only
    if doPSTHSubtraction
        X_fa_used = X_fa_resid;
    else
        X_fa_used = X_fa;
    end
end

if skipped_miss
    warning('[%s] Miss trial count (%d) < minCorrectTrials (%d). Excluding Miss from combinedAll/allGo for this session.', ...
        header, nTrials_miss, minCorrectTrials);
    X_miss_used = [];
else
    X_miss = Hs.Y3(trI.missI, :, :);
    [~, X_miss_resid] = computeResidual_perMotif(X_miss);   % psth_miss discarded -- internal use only
    if doPSTHSubtraction
        X_miss_used = X_miss_resid;
    else
        X_miss_used = X_miss;
    end
end

%% -------------------- 2c) Combined pools --------------------
% Each trial type keeps its OWN condition-specific PSTH subtracted (if
% doPSTHSubtraction=true) BEFORE pooling.

% combined_correct: Hit + CR only. Skipped entirely (not partially) if
% EITHER was skipped.
skipped_combinedCorrect = skipped_hit || skipped_cr;
if skipped_combinedCorrect
    X_combinedCorrect_used = [];
else
    X_combinedCorrect_used = cat(1, X_hit_used, X_cr_used);
end

% combinedAll: whichever of Hit/CR/FA/Miss individually cleared
% minCorrectTrials, pooled together (session-level PARTIAL inclusion is
% allowed here, unlike combined_correct/allGo/allNoGo).
includedTypes_combinedAll = {};
X_combinedAll_parts = {};
if ~skipped_hit,  X_combinedAll_parts{end+1} = X_hit_used;  includedTypes_combinedAll{end+1} = 'hit';  end %#ok<AGROW>
if ~skipped_cr,   X_combinedAll_parts{end+1} = X_cr_used;   includedTypes_combinedAll{end+1} = 'cr';   end %#ok<AGROW>
if ~skipped_fa,   X_combinedAll_parts{end+1} = X_fa_used;   includedTypes_combinedAll{end+1} = 'fa';   end %#ok<AGROW>
if ~skipped_miss, X_combinedAll_parts{end+1} = X_miss_used; includedTypes_combinedAll{end+1} = 'miss'; end %#ok<AGROW>

skipped_combinedAll = isempty(X_combinedAll_parts);
if skipped_combinedAll
    warning('[%s] All four trial types (Hit/CR/FA/Miss) skipped (below minCorrectTrials) -- combinedAll unavailable for this session.', header);
    X_combinedAll_used = [];
else
    X_combinedAll_used = cat(1, X_combinedAll_parts{:});
    if numel(includedTypes_combinedAll) < 2
        warning('[%s] combinedAll includes only 1 trial type (%s) this session -- effectively identical to that type alone, not a true combination.', ...
            header, includedTypes_combinedAll{1});
    end
end

% allGo (NEW): Hit + Miss -- every trial where a Go cue was presented,
% correct or not. Skipped entirely (not partially) if EITHER Hit or Miss
% was individually skipped -- same all-or-nothing logic as
% combined_correct, since "all Go trials" without one of its two
% constituents isn't a meaningful all-Go estimate.
skipped_allGo = skipped_hit || skipped_miss;
if skipped_allGo
    X_allGo_used = [];
else
    X_allGo_used = cat(1, X_hit_used, X_miss_used);
end

% allNoGo (NEW): CR + FA -- every trial where a NoGo cue was presented,
% correct or not. Same all-or-nothing logic.
skipped_allNoGo = skipped_cr || skipped_fa;
if skipped_allNoGo
    X_allNoGo_used = [];
else
    X_allNoGo_used = cat(1, X_cr_used, X_fa_used);
end

%% -------------------- 3) Observed per-trial averaged xcorr --------------------
% NOTE: variable/field names below are UNCHANGED regardless of
% doPSTHSubtraction mode -- check S.params.doPSTHSubtraction for the
% actual mode used.
if ~isempty(X_hit_used)
    XcorrMat_sess_hit_residual = computeMotifXcorr_perTrial_FisherZ(X_hit_used, nBin_total, useSymmetry, doFisherZ, clipR);
else
    XcorrMat_sess_hit_residual = nan(K, K, L);
end

if ~isempty(X_cr_used)
    XcorrMat_sess_cr_residual = computeMotifXcorr_perTrial_FisherZ(X_cr_used, nBin_total, useSymmetry, doFisherZ, clipR);
else
    XcorrMat_sess_cr_residual = nan(K, K, L);
end

if ~isempty(X_combinedCorrect_used)
    XcorrMat_sess_combinedCorrect_residual = computeMotifXcorr_perTrial_FisherZ(X_combinedCorrect_used, nBin_total, useSymmetry, doFisherZ, clipR);
else
    XcorrMat_sess_combinedCorrect_residual = nan(K, K, L);
end

if ~isempty(X_combinedAll_used)
    XcorrMat_sess_combinedAll_residual = computeMotifXcorr_perTrial_FisherZ(X_combinedAll_used, nBin_total, useSymmetry, doFisherZ, clipR);
else
    XcorrMat_sess_combinedAll_residual = nan(K, K, L);
end

if ~isempty(X_allGo_used)
    XcorrMat_sess_allGo_residual = computeMotifXcorr_perTrial_FisherZ(X_allGo_used, nBin_total, useSymmetry, doFisherZ, clipR);
else
    XcorrMat_sess_allGo_residual = nan(K, K, L);
end

if ~isempty(X_allNoGo_used)
    XcorrMat_sess_allNoGo_residual = computeMotifXcorr_perTrial_FisherZ(X_allNoGo_used, nBin_total, useSymmetry, doFisherZ, clipR);
else
    XcorrMat_sess_allNoGo_residual = nan(K, K, L);
end

xcorrPosLagMat_hit_residual = squeeze(mean(XcorrMat_sess_hit_residual(:, :, posLagMask_pool), 3, 'omitnan'));
xcorrPosLagMat_cr_residual  = squeeze(mean(XcorrMat_sess_cr_residual(:, :,  posLagMask_pool), 3, 'omitnan'));

xcorrPosLagMat_combinedCorrect_residual = squeeze(mean(XcorrMat_sess_combinedCorrect_residual(:, :, posLagMask_pool), 3, 'omitnan'));
xcorrPosLagMat_combinedAll_residual     = squeeze(mean(XcorrMat_sess_combinedAll_residual(:, :,     posLagMask_pool), 3, 'omitnan'));

xcorrPosLagMat_allGo_residual   = squeeze(mean(XcorrMat_sess_allGo_residual(:, :,   posLagMask_pool), 3, 'omitnan'));
xcorrPosLagMat_allNoGo_residual = squeeze(mean(XcorrMat_sess_allNoGo_residual(:, :, posLagMask_pool), 3, 'omitnan'));

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
    S.meta.nTrials_combinedAll = size(X_combinedAll_used, 1);
end
S.meta.skipped_combinedAll             = skipped_combinedAll;
S.meta.included_trialTypes_combinedAll = includedTypes_combinedAll;

if skipped_allGo
    S.meta.nTrials_allGo = NaN;
else
    S.meta.nTrials_allGo = nTrials_hit + nTrials_miss;
end
S.meta.skipped_allGo = skipped_allGo;

if skipped_allNoGo
    S.meta.nTrials_allNoGo = NaN;
else
    S.meta.nTrials_allNoGo = nTrials_cr + nTrials_fa;
end
S.meta.skipped_allNoGo = skipped_allNoGo;

S.params = struct();
S.params.xcorrLagWindowSec     = xcorrLagWindowSec;
S.params.posLagPoolWindowSec   = posLagPoolWindowSec;
S.params.stepSec               = stepSec;
S.params.nBin_total            = nBin_total;
S.params.nBin_pool             = nBin_pool;
S.params.lags                  = lags;
S.params.posLagMask_pool       = posLagMask_pool;

S.params.doTimeShuffle         = doTimeShuffle;
S.params.doTrialShuffle        = doTrialShuffle;
S.params.shuffleMethod         = shuffleMethod;
S.params.nShuffle               = nShuffle;
S.params.rngSeed                = rngSeed;
S.params.zscoreHs               = doZscore;
S.params.useSymmetry            = useSymmetry;
S.params.doFisherZ              = doFisherZ;
S.params.clipR                  = clipR;
S.params.minCorrectTrials       = minCorrectTrials;
S.params.doPSTHSubtraction      = doPSTHSubtraction;   % the definitive provenance flag for this run

S.params.showProgress           = showProgress;
S.params.progressEvery          = progressEvery;

S.obs = struct();
S.obs.psth_hit = psth_hit;   % [K x T] ALWAYS computed regardless of mode, for reference/comparison
S.obs.psth_cr  = psth_cr;
% NOTE: psth_fa / psth_miss deliberately not saved -- internal-only.

S.obs.XcorrMat_sess_hit_residual = XcorrMat_sess_hit_residual;
S.obs.XcorrMat_sess_cr_residual  = XcorrMat_sess_cr_residual;
S.obs.XcorrMat_sess_combinedCorrect_residual = XcorrMat_sess_combinedCorrect_residual;
S.obs.XcorrMat_sess_combinedAll_residual     = XcorrMat_sess_combinedAll_residual;
S.obs.XcorrMat_sess_allGo_residual   = XcorrMat_sess_allGo_residual;
S.obs.XcorrMat_sess_allNoGo_residual = XcorrMat_sess_allNoGo_residual;

S.obs.xcorrPosLagMat_hit_residual = xcorrPosLagMat_hit_residual;
S.obs.xcorrPosLagMat_cr_residual  = xcorrPosLagMat_cr_residual;
S.obs.xcorrPosLagMat_combinedCorrect_residual = xcorrPosLagMat_combinedCorrect_residual;
S.obs.xcorrPosLagMat_combinedAll_residual     = xcorrPosLagMat_combinedAll_residual;
S.obs.xcorrPosLagMat_allGo_residual   = xcorrPosLagMat_allGo_residual;
S.obs.xcorrPosLagMat_allNoGo_residual = xcorrPosLagMat_allNoGo_residual;

%% -------------------- 5) Optional shuffle null(s) + stats --------------------
% doTimeShuffle and doTrialShuffle are INDEPENDENT toggles -- if both are
% false, skip this whole section (observed-only, same as the original
% single-toggle behavior). If either is true, that null's accumulators
% below are used; the other's are simply left unused (still allocated,
% for simplicity/parfor-slicing consistency, but never written to if its
% toggle is off).
if ~doTimeShuffle && ~doTrialShuffle
    return;
end

% ---- within-trial shuffle accumulators (unchanged from before) ----
shuf_hit = nan(K, K, nShuffle);
shuf_cr  = nan(K, K, nShuffle);
shuf_combinedCorrect = nan(K, K, nShuffle);
shuf_combinedAll     = nan(K, K, nShuffle);
shuf_allGo   = nan(K, K, nShuffle);
shuf_allNoGo = nan(K, K, nShuffle);

sum_hit  = zeros(K, K, L);   sumsq_hit = zeros(K, K, L);   cnt_hit = zeros(K, K, L);
sum_cr   = zeros(K, K, L);   sumsq_cr  = zeros(K, K, L);   cnt_cr  = zeros(K, K, L);
sum_combinedCorrect = zeros(K, K, L);   sumsq_combinedCorrect = zeros(K, K, L);   cnt_combinedCorrect = zeros(K, K, L);
sum_combinedAll     = zeros(K, K, L);   sumsq_combinedAll     = zeros(K, K, L);   cnt_combinedAll     = zeros(K, K, L);
sum_allGo   = zeros(K, K, L);   sumsq_allGo   = zeros(K, K, L);   cnt_allGo   = zeros(K, K, L);
sum_allNoGo = zeros(K, K, L);   sumsq_allNoGo = zeros(K, K, L);   cnt_allNoGo = zeros(K, K, L);

% ---- NEW: trial-shuffle accumulators (parallel set, "_ts" suffix locally) ----
shuf_hit_ts = nan(K, K, nShuffle);
shuf_cr_ts  = nan(K, K, nShuffle);
shuf_combinedCorrect_ts = nan(K, K, nShuffle);
shuf_combinedAll_ts     = nan(K, K, nShuffle);
shuf_allGo_ts   = nan(K, K, nShuffle);
shuf_allNoGo_ts = nan(K, K, nShuffle);

sum_hit_ts  = zeros(K, K, L);   sumsq_hit_ts = zeros(K, K, L);   cnt_hit_ts = zeros(K, K, L);
sum_cr_ts   = zeros(K, K, L);   sumsq_cr_ts  = zeros(K, K, L);   cnt_cr_ts  = zeros(K, K, L);
sum_combinedCorrect_ts = zeros(K, K, L);   sumsq_combinedCorrect_ts = zeros(K, K, L);   cnt_combinedCorrect_ts = zeros(K, K, L);
sum_combinedAll_ts     = zeros(K, K, L);   sumsq_combinedAll_ts     = zeros(K, K, L);   cnt_combinedAll_ts     = zeros(K, K, L);
sum_allGo_ts   = zeros(K, K, L);   sumsq_allGo_ts   = zeros(K, K, L);   cnt_allGo_ts   = zeros(K, K, L);
sum_allNoGo_ts = zeros(K, K, L);   sumsq_allNoGo_ts = zeros(K, K, L);   cnt_allNoGo_ts = zeros(K, K, L);

try
    pctRunOnAll rng('shuffle');
catch
end

dq = [];
tStart = tic;
if showProgress
    dq = parallel.pool.DataQueue;
    nDone = 0; % on client
    afterEach(dq, @updateProgress);
end

parfor s = 1:nShuffle
    % ---------------- HIT ----------------
    if doTimeShuffle && ~isempty(X_hit_used)
        Xs_hit = withinTrialShuffle_independent(X_hit_used, shuffleMethod);
        Xc_hit = computeMotifXcorr_perTrial_FisherZ(Xs_hit, nBin_total, useSymmetry, doFisherZ, clipR);

        shuf_hit(:, :, s) = squeeze(mean(Xc_hit(:, :, posLagMask_pool), 3, 'omitnan'));

        m = isfinite(Xc_hit);
        Xc0 = Xc_hit; Xc0(~m) = 0;
        sum_hit   = sum_hit   + Xc0;
        sumsq_hit = sumsq_hit + Xc0.^2;
        cnt_hit   = cnt_hit   + double(m);
    end
    % ---------------- HIT (NEW: trial shuffle) ----------------
    if doTrialShuffle && ~isempty(X_hit_used)
        Xs_hit_ts = trialShuffle_independent(X_hit_used);
        Xc_hit_ts = computeMotifXcorr_perTrial_FisherZ(Xs_hit_ts, nBin_total, useSymmetry, doFisherZ, clipR);

        shuf_hit_ts(:, :, s) = squeeze(mean(Xc_hit_ts(:, :, posLagMask_pool), 3, 'omitnan'));

        m = isfinite(Xc_hit_ts);
        Xc0 = Xc_hit_ts; Xc0(~m) = 0;
        sum_hit_ts   = sum_hit_ts   + Xc0;
        sumsq_hit_ts = sumsq_hit_ts + Xc0.^2;
        cnt_hit_ts   = cnt_hit_ts   + double(m);
    end

    % ---------------- CR ----------------
    if doTimeShuffle && ~isempty(X_cr_used)
        Xs_cr = withinTrialShuffle_independent(X_cr_used, shuffleMethod);
        Xc_cr = computeMotifXcorr_perTrial_FisherZ(Xs_cr, nBin_total, useSymmetry, doFisherZ, clipR);

        shuf_cr(:, :, s) = squeeze(mean(Xc_cr(:, :, posLagMask_pool), 3, 'omitnan'));

        m = isfinite(Xc_cr);
        Xc0 = Xc_cr; Xc0(~m) = 0;
        sum_cr   = sum_cr   + Xc0;
        sumsq_cr = sumsq_cr + Xc0.^2;
        cnt_cr   = cnt_cr   + double(m);
    end
    % ---------------- CR (NEW: trial shuffle) ----------------
    if doTrialShuffle && ~isempty(X_cr_used)
        Xs_cr_ts = trialShuffle_independent(X_cr_used);
        Xc_cr_ts = computeMotifXcorr_perTrial_FisherZ(Xs_cr_ts, nBin_total, useSymmetry, doFisherZ, clipR);

        shuf_cr_ts(:, :, s) = squeeze(mean(Xc_cr_ts(:, :, posLagMask_pool), 3, 'omitnan'));

        m = isfinite(Xc_cr_ts);
        Xc0 = Xc_cr_ts; Xc0(~m) = 0;
        sum_cr_ts   = sum_cr_ts   + Xc0;
        sumsq_cr_ts = sumsq_cr_ts + Xc0.^2;
        cnt_cr_ts   = cnt_cr_ts   + double(m);
    end

    % ---------------- combined_correct ----------------
    % Shuffling the already-pooled trace pool is equivalent to shuffling
    % each trial independently within its own type first and then pooling
    % -- withinTrialShuffle_independent/trialShuffle_independent only
    % care about the trial dimension, not which type each trial came
    % from. Applies to every pooled stream below too.
    if doTimeShuffle && ~isempty(X_combinedCorrect_used)
        Xs_cc = withinTrialShuffle_independent(X_combinedCorrect_used, shuffleMethod);
        Xc_cc = computeMotifXcorr_perTrial_FisherZ(Xs_cc, nBin_total, useSymmetry, doFisherZ, clipR);

        shuf_combinedCorrect(:, :, s) = squeeze(mean(Xc_cc(:, :, posLagMask_pool), 3, 'omitnan'));

        m = isfinite(Xc_cc);
        Xc0 = Xc_cc; Xc0(~m) = 0;
        sum_combinedCorrect   = sum_combinedCorrect   + Xc0;
        sumsq_combinedCorrect = sumsq_combinedCorrect + Xc0.^2;
        cnt_combinedCorrect   = cnt_combinedCorrect   + double(m);
    end
    % ---------------- combined_correct (NEW: trial shuffle) ----------------
    if doTrialShuffle && ~isempty(X_combinedCorrect_used)
        Xs_cc_ts = trialShuffle_independent(X_combinedCorrect_used);
        Xc_cc_ts = computeMotifXcorr_perTrial_FisherZ(Xs_cc_ts, nBin_total, useSymmetry, doFisherZ, clipR);

        shuf_combinedCorrect_ts(:, :, s) = squeeze(mean(Xc_cc_ts(:, :, posLagMask_pool), 3, 'omitnan'));

        m = isfinite(Xc_cc_ts);
        Xc0 = Xc_cc_ts; Xc0(~m) = 0;
        sum_combinedCorrect_ts   = sum_combinedCorrect_ts   + Xc0;
        sumsq_combinedCorrect_ts = sumsq_combinedCorrect_ts + Xc0.^2;
        cnt_combinedCorrect_ts   = cnt_combinedCorrect_ts   + double(m);
    end

    % ---------------- combinedAll ----------------
    if doTimeShuffle && ~isempty(X_combinedAll_used)
        Xs_ca = withinTrialShuffle_independent(X_combinedAll_used, shuffleMethod);
        Xc_ca = computeMotifXcorr_perTrial_FisherZ(Xs_ca, nBin_total, useSymmetry, doFisherZ, clipR);

        shuf_combinedAll(:, :, s) = squeeze(mean(Xc_ca(:, :, posLagMask_pool), 3, 'omitnan'));

        m = isfinite(Xc_ca);
        Xc0 = Xc_ca; Xc0(~m) = 0;
        sum_combinedAll   = sum_combinedAll   + Xc0;
        sumsq_combinedAll = sumsq_combinedAll + Xc0.^2;
        cnt_combinedAll   = cnt_combinedAll   + double(m);
    end
    % ---------------- combinedAll (NEW: trial shuffle) ----------------
    if doTrialShuffle && ~isempty(X_combinedAll_used)
        Xs_ca_ts = trialShuffle_independent(X_combinedAll_used);
        Xc_ca_ts = computeMotifXcorr_perTrial_FisherZ(Xs_ca_ts, nBin_total, useSymmetry, doFisherZ, clipR);

        shuf_combinedAll_ts(:, :, s) = squeeze(mean(Xc_ca_ts(:, :, posLagMask_pool), 3, 'omitnan'));

        m = isfinite(Xc_ca_ts);
        Xc0 = Xc_ca_ts; Xc0(~m) = 0;
        sum_combinedAll_ts   = sum_combinedAll_ts   + Xc0;
        sumsq_combinedAll_ts = sumsq_combinedAll_ts + Xc0.^2;
        cnt_combinedAll_ts   = cnt_combinedAll_ts   + double(m);
    end

    % ---------------- allGo ----------------
    if doTimeShuffle && ~isempty(X_allGo_used)
        Xs_ag = withinTrialShuffle_independent(X_allGo_used, shuffleMethod);
        Xc_ag = computeMotifXcorr_perTrial_FisherZ(Xs_ag, nBin_total, useSymmetry, doFisherZ, clipR);

        shuf_allGo(:, :, s) = squeeze(mean(Xc_ag(:, :, posLagMask_pool), 3, 'omitnan'));

        m = isfinite(Xc_ag);
        Xc0 = Xc_ag; Xc0(~m) = 0;
        sum_allGo   = sum_allGo   + Xc0;
        sumsq_allGo = sumsq_allGo + Xc0.^2;
        cnt_allGo   = cnt_allGo   + double(m);
    end
    % ---------------- allGo (NEW: trial shuffle) ----------------
    if doTrialShuffle && ~isempty(X_allGo_used)
        Xs_ag_ts = trialShuffle_independent(X_allGo_used);
        Xc_ag_ts = computeMotifXcorr_perTrial_FisherZ(Xs_ag_ts, nBin_total, useSymmetry, doFisherZ, clipR);

        shuf_allGo_ts(:, :, s) = squeeze(mean(Xc_ag_ts(:, :, posLagMask_pool), 3, 'omitnan'));

        m = isfinite(Xc_ag_ts);
        Xc0 = Xc_ag_ts; Xc0(~m) = 0;
        sum_allGo_ts   = sum_allGo_ts   + Xc0;
        sumsq_allGo_ts = sumsq_allGo_ts + Xc0.^2;
        cnt_allGo_ts   = cnt_allGo_ts   + double(m);
    end

    % ---------------- allNoGo ----------------
    if doTimeShuffle && ~isempty(X_allNoGo_used)
        Xs_ang = withinTrialShuffle_independent(X_allNoGo_used, shuffleMethod);
        Xc_ang = computeMotifXcorr_perTrial_FisherZ(Xs_ang, nBin_total, useSymmetry, doFisherZ, clipR);

        shuf_allNoGo(:, :, s) = squeeze(mean(Xc_ang(:, :, posLagMask_pool), 3, 'omitnan'));

        m = isfinite(Xc_ang);
        Xc0 = Xc_ang; Xc0(~m) = 0;
        sum_allNoGo   = sum_allNoGo   + Xc0;
        sumsq_allNoGo = sumsq_allNoGo + Xc0.^2;
        cnt_allNoGo   = cnt_allNoGo   + double(m);
    end
    % ---------------- allNoGo (NEW: trial shuffle) ----------------
    if doTrialShuffle && ~isempty(X_allNoGo_used)
        Xs_ang_ts = trialShuffle_independent(X_allNoGo_used);
        Xc_ang_ts = computeMotifXcorr_perTrial_FisherZ(Xs_ang_ts, nBin_total, useSymmetry, doFisherZ, clipR);

        shuf_allNoGo_ts(:, :, s) = squeeze(mean(Xc_ang_ts(:, :, posLagMask_pool), 3, 'omitnan'));

        m = isfinite(Xc_ang_ts);
        Xc0 = Xc_ang_ts; Xc0(~m) = 0;
        sum_allNoGo_ts   = sum_allNoGo_ts   + Xc0;
        sumsq_allNoGo_ts = sumsq_allNoGo_ts + Xc0.^2;
        cnt_allNoGo_ts   = cnt_allNoGo_ts   + double(m);
    end

    if showProgress
        send(dq, s);
    end
end

S.shuf = struct();

% ---- within-trial shuffle results (unchanged fields/behavior) ----
if doTimeShuffle
    [mu_hit, sd_hit] = finalizeMeanStd(sum_hit, sumsq_hit, cnt_hit);
    [mu_cr,  sd_cr]  = finalizeMeanStd(sum_cr,  sumsq_cr,  cnt_cr);
    [mu_combinedCorrect, sd_combinedCorrect] = finalizeMeanStd(sum_combinedCorrect, sumsq_combinedCorrect, cnt_combinedCorrect);
    [mu_combinedAll,     sd_combinedAll]     = finalizeMeanStd(sum_combinedAll,     sumsq_combinedAll,     cnt_combinedAll);
    [mu_allGo,   sd_allGo]   = finalizeMeanStd(sum_allGo,   sumsq_allGo,   cnt_allGo);
    [mu_allNoGo, sd_allNoGo] = finalizeMeanStd(sum_allNoGo, sumsq_allNoGo, cnt_allNoGo);

    S.shuf.xcorrPosLagMat_hit_residual = shuf_hit;
    S.shuf.xcorrPosLagMat_cr_residual  = shuf_cr;
    S.shuf.xcorrPosLagMat_combinedCorrect_residual = shuf_combinedCorrect;
    S.shuf.xcorrPosLagMat_combinedAll_residual     = shuf_combinedAll;
    S.shuf.xcorrPosLagMat_allGo_residual   = shuf_allGo;
    S.shuf.xcorrPosLagMat_allNoGo_residual = shuf_allNoGo;

    [S.shuf.mean_hit_residual, S.shuf.std_hit_residual, S.shuf.z_hit_residual, S.shuf.p_hit_residual] = ...
        shufStats(xcorrPosLagMat_hit_residual, shuf_hit);
    [S.shuf.mean_cr_residual,  S.shuf.std_cr_residual,  S.shuf.z_cr_residual,  S.shuf.p_cr_residual]  = ...
        shufStats(xcorrPosLagMat_cr_residual,  shuf_cr);
    [S.shuf.mean_combinedCorrect_residual, S.shuf.std_combinedCorrect_residual, S.shuf.z_combinedCorrect_residual, S.shuf.p_combinedCorrect_residual] = ...
        shufStats(xcorrPosLagMat_combinedCorrect_residual, shuf_combinedCorrect);
    [S.shuf.mean_combinedAll_residual, S.shuf.std_combinedAll_residual, S.shuf.z_combinedAll_residual, S.shuf.p_combinedAll_residual] = ...
        shufStats(xcorrPosLagMat_combinedAll_residual, shuf_combinedAll);
    [S.shuf.mean_allGo_residual, S.shuf.std_allGo_residual, S.shuf.z_allGo_residual, S.shuf.p_allGo_residual] = ...
        shufStats(xcorrPosLagMat_allGo_residual, shuf_allGo);
    [S.shuf.mean_allNoGo_residual, S.shuf.std_allNoGo_residual, S.shuf.z_allNoGo_residual, S.shuf.p_allNoGo_residual] = ...
        shufStats(xcorrPosLagMat_allNoGo_residual, shuf_allNoGo);

    S.shuf.curve.mean_hit_residual = mu_hit;
    S.shuf.curve.std_hit_residual  = sd_hit;
    S.shuf.curve.mean_cr_residual  = mu_cr;
    S.shuf.curve.std_cr_residual   = sd_cr;
    S.shuf.curve.mean_combinedCorrect_residual = mu_combinedCorrect;
    S.shuf.curve.std_combinedCorrect_residual  = sd_combinedCorrect;
    S.shuf.curve.mean_combinedAll_residual     = mu_combinedAll;
    S.shuf.curve.std_combinedAll_residual      = sd_combinedAll;
    S.shuf.curve.mean_allGo_residual   = mu_allGo;
    S.shuf.curve.std_allGo_residual    = sd_allGo;
    S.shuf.curve.mean_allNoGo_residual = mu_allNoGo;
    S.shuf.curve.std_allNoGo_residual  = sd_allNoGo;
end

% ---- NEW: trial-shuffle results (parallel set, "_trialShuffle" suffix on
% top of the existing "_residual" suffix, so field names never collide
% with the within-trial-shuffle ones above) ----
if doTrialShuffle
    [mu_hit_ts, sd_hit_ts] = finalizeMeanStd(sum_hit_ts, sumsq_hit_ts, cnt_hit_ts);
    [mu_cr_ts,  sd_cr_ts]  = finalizeMeanStd(sum_cr_ts,  sumsq_cr_ts,  cnt_cr_ts);
    [mu_combinedCorrect_ts, sd_combinedCorrect_ts] = finalizeMeanStd(sum_combinedCorrect_ts, sumsq_combinedCorrect_ts, cnt_combinedCorrect_ts);
    [mu_combinedAll_ts,     sd_combinedAll_ts]     = finalizeMeanStd(sum_combinedAll_ts,     sumsq_combinedAll_ts,     cnt_combinedAll_ts);
    [mu_allGo_ts,   sd_allGo_ts]   = finalizeMeanStd(sum_allGo_ts,   sumsq_allGo_ts,   cnt_allGo_ts);
    [mu_allNoGo_ts, sd_allNoGo_ts] = finalizeMeanStd(sum_allNoGo_ts, sumsq_allNoGo_ts, cnt_allNoGo_ts);

    S.shuf.xcorrPosLagMat_hit_residual_trialShuffle = shuf_hit_ts;
    S.shuf.xcorrPosLagMat_cr_residual_trialShuffle  = shuf_cr_ts;
    S.shuf.xcorrPosLagMat_combinedCorrect_residual_trialShuffle = shuf_combinedCorrect_ts;
    S.shuf.xcorrPosLagMat_combinedAll_residual_trialShuffle     = shuf_combinedAll_ts;
    S.shuf.xcorrPosLagMat_allGo_residual_trialShuffle   = shuf_allGo_ts;
    S.shuf.xcorrPosLagMat_allNoGo_residual_trialShuffle = shuf_allNoGo_ts;

    [S.shuf.mean_hit_residual_trialShuffle, S.shuf.std_hit_residual_trialShuffle, S.shuf.z_hit_residual_trialShuffle, S.shuf.p_hit_residual_trialShuffle] = ...
        shufStats(xcorrPosLagMat_hit_residual, shuf_hit_ts);
    [S.shuf.mean_cr_residual_trialShuffle,  S.shuf.std_cr_residual_trialShuffle,  S.shuf.z_cr_residual_trialShuffle,  S.shuf.p_cr_residual_trialShuffle]  = ...
        shufStats(xcorrPosLagMat_cr_residual,  shuf_cr_ts);
    [S.shuf.mean_combinedCorrect_residual_trialShuffle, S.shuf.std_combinedCorrect_residual_trialShuffle, S.shuf.z_combinedCorrect_residual_trialShuffle, S.shuf.p_combinedCorrect_residual_trialShuffle] = ...
        shufStats(xcorrPosLagMat_combinedCorrect_residual, shuf_combinedCorrect_ts);
    [S.shuf.mean_combinedAll_residual_trialShuffle, S.shuf.std_combinedAll_residual_trialShuffle, S.shuf.z_combinedAll_residual_trialShuffle, S.shuf.p_combinedAll_residual_trialShuffle] = ...
        shufStats(xcorrPosLagMat_combinedAll_residual, shuf_combinedAll_ts);
    [S.shuf.mean_allGo_residual_trialShuffle, S.shuf.std_allGo_residual_trialShuffle, S.shuf.z_allGo_residual_trialShuffle, S.shuf.p_allGo_residual_trialShuffle] = ...
        shufStats(xcorrPosLagMat_allGo_residual, shuf_allGo_ts);
    [S.shuf.mean_allNoGo_residual_trialShuffle, S.shuf.std_allNoGo_residual_trialShuffle, S.shuf.z_allNoGo_residual_trialShuffle, S.shuf.p_allNoGo_residual_trialShuffle] = ...
        shufStats(xcorrPosLagMat_allNoGo_residual, shuf_allNoGo_ts);

    S.shuf.curve.mean_hit_residual_trialShuffle = mu_hit_ts;
    S.shuf.curve.std_hit_residual_trialShuffle  = sd_hit_ts;
    S.shuf.curve.mean_cr_residual_trialShuffle  = mu_cr_ts;
    S.shuf.curve.std_cr_residual_trialShuffle   = sd_cr_ts;
    S.shuf.curve.mean_combinedCorrect_residual_trialShuffle = mu_combinedCorrect_ts;
    S.shuf.curve.std_combinedCorrect_residual_trialShuffle  = sd_combinedCorrect_ts;
    S.shuf.curve.mean_combinedAll_residual_trialShuffle     = mu_combinedAll_ts;
    S.shuf.curve.std_combinedAll_residual_trialShuffle      = sd_combinedAll_ts;
    S.shuf.curve.mean_allGo_residual_trialShuffle   = mu_allGo_ts;
    S.shuf.curve.std_allGo_residual_trialShuffle    = sd_allGo_ts;
    S.shuf.curve.mean_allNoGo_residual_trialShuffle = mu_allNoGo_ts;
    S.shuf.curve.std_allNoGo_residual_trialShuffle  = sd_allNoGo_ts;
end

if showProgress
    fprintf('[%s] Completed %d/%d shuffles in %.1f sec (doPSTHSubtraction=%d, doTimeShuffle=%d, doTrialShuffle=%d)\n', ...
        header, nShuffle, nShuffle, toc(tStart), doPSTHSubtraction, doTimeShuffle, doTrialShuffle);
end

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
% Plain across-trial mean per motif (the PSTH/template), and the residual
% after subtracting it from every trial. ALWAYS computed regardless of
% doPSTHSubtraction -- the caller decides whether Xresid actually gets
% used downstream, or whether the raw X is used instead.
psth = squeeze(mean(X, 1, 'omitnan'));
if size(X, 2) == 1
    psth = reshape(psth, 1, []);
end
Xresid = X - reshape(psth, 1, size(psth,1), size(psth,2));
end

%% ========================================================================
function Xs = trialShuffle_independent(X)
% TRIALSHUFFLE_INDEPENDENT
%   X: [N x K x T]. For EACH MOTIF independently, randomly permutes WHICH
%   TRIAL its data comes from -- the mirror image of
%   withinTrialShuffle_independent below. Each trial's own internal time
%   series is left completely untouched (no shift, no permutation of
%   time points within a trial); what changes is which trial index gets
%   paired with which for a given motif. This destroys the trial-to-trial
%   PAIRING between motifs (and therefore any shared trial-level state,
%   e.g. a trial being "high engagement" for many motifs at once) while
%   preserving each individual trace's own fine-timescale structure --
%   the opposite of what circular/permute shifting destroys and
%   preserves. See the header-comment null-hypothesis discussion for the
%   distinction this makes in practice.
[N, K, T] = size(X); %#ok<ASGLU> % T unused here but kept for symmetry/readability with the sibling function
Xs = X;
for k = 1:K
    permIdx = randperm(N);
    Xs(:, k, :) = X(permIdx, k, :);
end
end

%% ========================================================================
function Xs = withinTrialShuffle_independent(X, method)
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
v(v < 0 & v > -1e-12) = 0;
sd(valid) = sqrt(v(valid));
end

%% ------------------------------------------------------------------------
function [mu, sig, z, p] = shufStats(obsMat, shufMat)
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
    Xn = squeeze(HsY3(n, :, :));
    if any(~isfinite(Xn(:)))
        Xn(~isfinite(Xn)) = 0;
    end

    if useSymmetry
        for i = 1:K
            xi = Xn(i, :).';
            xi = xi - mean(xi, 'omitnan');
            for j = i:K
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