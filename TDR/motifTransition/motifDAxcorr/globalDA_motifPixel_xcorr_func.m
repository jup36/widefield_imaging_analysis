function S = globalDA_motifPixel_xcorr_func(fileDA, fileTbytDA, fileH, fileW, varargin)
%GLOBALDA_MOTIFPIXEL_XCORR_FUNC
%   Cross-correlation between cortical dopamine and calcium motif activity,
%   computed PER PIXEL and then reduced to one correlogram per motif by a
%   footprint-weighted spatial average.
%
%   WHY PER-PIXEL THEN AVERAGE (not average then correlate)
%   -------------------------------------------------------
%   For pixel-indexed signals A and B,
%       mean_p(A_p B_p) = mean_p(A) * mean_p(B) + cov_p(A_p, B_p)
%   Correlating first and averaging second KEEPS that covariance term --
%   the part that is large only when DA's spatial pattern at a given lag
%   matches the motif's. Averaging first destroys it. That term is the
%   whole reason to work in pixel space, and it is what makes output (1)
%   below differ from output (2).
%
%   THREE CORRELOGRAMS PER MOTIF AND STREAM, all [K x L]
%     XcorrMat_pixel_<stream>   per-pixel xcorr, footprint-weighted mean
%     XcorrMat_recon_<stream>   footprint-weighted spatial mean of the
%                               reconstruction, vs. spatially averaged DA
%     XcorrMat_H_<stream>       raw H vs. spatially averaged DA
%   All three come from the SAME trials, the same PSTH subtraction, and the
%   same shuffle draws, so any difference between them is attributable to
%   the spatial treatment and not to a different sample or seed.
%
%   RECONSTRUCTION AND LAG UNITS
%   ----------------------------
%   X_k(p,t) = sum_l W(p,k,l) H_k(t-l). W's lag index l is in units of H's
%   NATIVE frame period, so the convolution is performed on H's own grid
%   and only the RESULT is resampled to the trial grid. Convolving after
%   resampling would silently reinterpret each lag as 50 ms.
%
%   NULLS
%   -----
%   trial shuffle (primary): permute which trial the DA comes from. There
%     are only N^2 distinct trial pairings, so the full cross-trial tensor
%     C[nDA, nX, lag] is computed ONCE and every draw is a sum over one
%     permutation's entries. nShuffle is therefore free -- 1000 costs the
%     same as 10. (Fisher-z is nonlinear, so the z-transform happens per
%     pairing inside the tensor, not afterwards.)
%   circshift (check): independently shift each signal in time within each
%     trial. Cannot be cached -- every draw is a new signal -- so it is run
%     at lower nShuffle and, by default, only on a subset of motifs. Its
%     job is to catch coupling driven by residual task locking that
%     survives PSTH subtraction, which the trial shuffle cannot detect.
%     All pixels share one shift per signal per trial, preserving DA's
%     spatial pattern so this null stays comparable to the scalar case.
%
%   S = globalDA_motifPixel_xcorr_func(fileDA, fileH, fileW, ...)
%
% INPUTS (paths; bucket-style paths are fine on the cluster)
%   fileDA     : <header>_green_dff_smCollect.mat -- dffsmCell only.
%                NOTE this file is written WITHOUT -append by
%                dffPostprocess_auditory_gng_dual, so re-running the
%                postprocess wipes anything appended to it later (this is
%                exactly what removed DAglobalC/globalDA_trI). Nothing
%                here relies on appended variables for that reason.
%   fileTbytDA : <header>_green_tbytDat_dff.mat -- tbytDat with frameT.
%                frameT is written by dffCombinedProcessingDual in the same
%                pass that produced dffsmCell, as the LED times of exactly
%                the frames that went into dffsm, so the two cannot drift.
%   fileH      : refit file -- hC, tbytDat
%   fileW      : common basis file -- W_basis [P x K x L], nanpxs
%
%   trI is derived here from the refit file's tbytDat (see local
%   trialTypeInfoFromTbytDat) rather than loaded, so no fifth path and no
%   dependence on a previously appended variable.
%
% NAME-VALUE
%   'hRow', 'tRow'        : rows of hC for H and timestamps (1, 3)
%   'imageSize'           : [64 64]
%   'alignWin', 'alignDt' : [-0.9 5], 0.05   (50 ms: H's native step)
%   'maxLagSec'           : 2.0
%   'evtField'            : 'evtOn'
%   'pixelSelect'         : 'footprint' (default) | 'all' | 'topFrac'
%   'topFrac'             : 0.20
%   'doPSTHSubtraction'   : true
%   'doFisherZ', 'clipR'  : true, 0.999
%   'minCorrectTrials'    : 10
%   'nShuffleTrial'       : 1000   (free -- see above)
%   'doCircShuffle'       : true
%   'nShuffleCirc'        : 200
%   'circMotifs'          : []  -> all motifs; pass e.g. [goIdx nogoIdx]
%   'rngSeed'             : 1
%   'useParfor'           : true
%   'verbose'             : true

%% -------------------- parse --------------------
p = inputParser;
p.addParameter('hRow', 1, @(x) isnumeric(x) && isscalar(x));
p.addParameter('tRow', 3, @(x) isnumeric(x) && isscalar(x));
p.addParameter('imageSize', [64 64], @(x) isnumeric(x) && numel(x)==2);
p.addParameter('alignWin', [-0.9 5], @(x) isnumeric(x) && numel(x)==2);
p.addParameter('alignDt', 0.05, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('maxLagSec', 2.0, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('evtField', 'evtOn', @(s) ischar(s) || isstring(s));
p.addParameter('pixelSelect', 'footprint', @(s) any(strcmpi(string(s), ["all","footprint","topfrac"])));
p.addParameter('topFrac', 0.20, @(x) isnumeric(x) && isscalar(x) && x>0 && x<=1);
p.addParameter('doPSTHSubtraction', true, @(x) islogical(x) && isscalar(x));
p.addParameter('doFisherZ', true, @(x) islogical(x) && isscalar(x));
p.addParameter('clipR', 0.999, @(x) isnumeric(x) && isscalar(x) && x>0 && x<1);
p.addParameter('minCorrectTrials', 10, @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('nShuffleTrial', 1000, @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('doCircShuffle', true, @(x) islogical(x) && isscalar(x));
p.addParameter('nShuffleCirc', 200, @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('circMotifs', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('rngSeed', 1, @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('useParfor', true, @(x) islogical(x) && isscalar(x));
p.addParameter('verbose', true, @(x) islogical(x) && isscalar(x));
p.addParameter('saveDir', '', @(s) ischar(s) || isstring(s));
p.addParameter('saveKeyword', 'globalDA_motifPixel_xcorr', @(s) ischar(s) || isstring(s));
p.parse(varargin{:});
opt = p.Results;
pixelSelect = lower(string(opt.pixelSelect));
if ~isempty(opt.rngSeed), rng(opt.rngSeed); end

tAll = tic;

%% -------------------- load --------------------
vprintf(opt.verbose, '\n=== globalDA_motifPixel_xcorr_func ===\n');
vprintf(opt.verbose, 'DA : %s\nH  : %s\nW  : %s\n', fileDA, fileH, fileW);

dA = load(fileDA, 'dffsmCell');
assert(isfield(dA, 'dffsmCell'), 'fileDA is missing "dffsmCell".');
dT = load(fileTbytDA, 'tbytDat');
assert(isfield(dT, 'tbytDat'), 'fileTbytDA is missing "tbytDat".');
assert(isfield(dT.tbytDat, 'frameT'), ...
    ['fileTbytDA has no frameT field -- pass the <header>_green_tbytDat_dff.mat ' ...
     'written by dffPostprocess, not the parseGng tbytDat.']);
dH = load(fileH, 'hC', 'tbytDat');
assert(isfield(dH,'hC') && isfield(dH,'tbytDat'), 'fileH must contain hC and tbytDat.');
dW = load(fileW, 'W_basis', 'nanpxs');
assert(isfield(dW,'W_basis') && isfield(dW,'nanpxs'), 'fileW must contain W_basis and nanpxs.');

dffsmCell  = dA.dffsmCell;
tbytDA     = dT.tbytDat;      % green channel: frameT for the DA frames
hC         = dH.hC;
tbytDat    = dH.tbytDat;      % refit file: event times for the motif side
W_basis    = dW.W_basis;
nanpxs     = dW.nanpxs;

trI = trialTypeInfoFromTbytDat(tbytDat);

[P, K, L] = size(W_basis);
assert(P == prod(opt.imageSize), 'W_basis has %d pixels; imageSize implies %d.', P, prod(opt.imageSize));
nTrials = numel(tbytDat);
assert(numel(dffsmCell) == nTrials, 'dffsmCell has %d trials, tbytDat %d.', numel(dffsmCell), nTrials);
assert(numel(tbytDA) == nTrials, 'green tbytDat has %d trials, refit tbytDat %d.', numel(tbytDA), nTrials);

% The two tbytDat copies come from different channels' postprocess runs.
% They should describe the same session, so a disagreement in event times
% means the files are mismatched -- catch it here rather than producing a
% silently misaligned result.
evtH  = arrayfun(@(s) numScalar(s.(char(opt.evtField))), tbytDat(:));
evtDA = arrayfun(@(s) numScalar(s.(char(opt.evtField))), tbytDA(:));
bothOK = isfinite(evtH) & isfinite(evtDA);
maxDiff = max(abs(evtH(bothOK) - evtDA(bothOK)));
assert(isempty(maxDiff) || maxDiff < 1e-6, ...
    ['%s differs between the refit and green tbytDat by up to %.4g s -- ' ...
     'these files describe different sessions or different alignments.'], ...
    char(opt.evtField), maxDiff);

validPix = true(P,1); validPix(nanpxs(:)) = false;
pv = find(validPix); Pv = numel(pv);

tint  = opt.alignWin(1):opt.alignDt:opt.alignWin(2);
nTime = numel(tint);
maxLag = round(opt.maxLagSec / opt.alignDt);
lagSec = (-maxLag:maxLag) * opt.alignDt;
Lg = 2*maxLag + 1;
assert(maxLag < nTime, 'maxLagSec exceeds the trial window.');

vprintf(opt.verbose, 'K=%d motifs, L=%d W-lags, %d/%d valid pixels, %d bins @ %g s, %d lags.\n', ...
    K, L, Pv, P, nTime, opt.alignDt, Lg);

%% -------------------- footprint weights --------------------
wSum = squeeze(sum(abs(W_basis), 3));     % P x K
Wt = zeros(Pv, K);
for k = 1:K
    w = zeros(P,1);
    switch pixelSelect
        case "all",       w(validPix) = 1;
        case "footprint", w = wSum(:,k) .* validPix;
        case "topfrac"
            fp = wSum(:,k); fp(~validPix) = -inf;
            nTop = max(1, round(opt.topFrac * Pv));
            [~, ord] = sort(fp,'descend'); w(ord(1:nTop)) = 1;
    end
    if sum(w) > 0, w = w / sum(w); end
    Wt(:,k) = w(pv);
end

%% -------------------- align pixel-space DA once --------------------
% Reused by every motif, so it is built once. [Pv x nTime x nTrials] single
% is ~350 MB at Pv=2410, nTime=119, nTrials=300.
vprintf(opt.verbose, 'Aligning pixel-space DA...\n');
DA3 = nan(Pv, nTime, nTrials, 'single');
nFrameMismatch = 0;
for tr = 1:nTrials
    X = dffsmCell{tr};
    tAbs = numVector(tbytDA(tr).frameT);          % LED times of exactly these frames
    if isempty(X) || isempty(tAbs), continue; end
    nUse = min(size(X,3), numel(tAbs));
    if size(X,3) ~= numel(tAbs), nFrameMismatch = nFrameMismatch + 1; end
    if nUse < 2, continue; end
    evt = numScalar(tbytDA(tr).(char(opt.evtField)));
    if ~isfinite(evt), continue; end

    M = reshape(double(X(:,:,1:nUse)), P, nUse);
    M = M(pv, :);
    DA3(:,:,tr) = single(interpToGrid(M, tAbs(1:nUse) - evt, tint));   % already trial-local
end
if nFrameMismatch > 0
    warning('motifPixelXcorr:FrameTMismatch', ...
        ['%d/%d trials have size(dffsm,3) ~= numel(frameT) (truncated to the shorter). ' ...
         'These should match exactly -- check that fileDA and fileTbytDA came from the same run.'], ...
        nFrameMismatch, nTrials);
end

%% -------------------- reconstruct motifs on H's native grid ------------
% Convolution FIRST on H's own frame period (W's lag unit), resampling
% SECOND. Blocks are concatenated on absolute time and sorted once.
vprintf(opt.verbose, 'Building motif reconstruction (native grid)...\n');
blocks = find(~cellfun(@isempty, hC(opt.hRow, :)));
Hcat = []; Tcat = [];
for b = blocks
    Hb = double(hC{opt.hRow, b});
    tb = double(hC{opt.tRow, b}(:)');
    n = min(size(Hb,2), numel(tb));
    Hcat = [Hcat, Hb(:,1:n)];  Tcat = [Tcat, tb(1:n)];  %#ok<AGROW>
end
assert(~isempty(Tcat), 'No usable H frames.');
[Tcat, si] = sort(Tcat); Hcat = Hcat(:, si);
[Tcat, ui] = unique(Tcat, 'stable'); Hcat = Hcat(:, ui);

%% -------------------- streams --------------------
streams = {'hit','cr'};
masks   = {logical(trI.hitI(:)'), logical(trI.crI(:)')};

S = struct();
S.meta = struct('K',K,'L',L,'Lg',Lg,'lagSec',lagSec,'tint',tint,'nTime',nTime, ...
    'nTrialsTotal',nTrials,'nValidPix',Pv,'imageSize',opt.imageSize, ...
    'direction','positive lag = motif leads DA');
S.params = opt;
S.params.fileDA = fileDA; S.params.fileH = fileH; S.params.fileW = fileW;
S.obs = struct(); S.shuf = struct();

circMotifs = opt.circMotifs; if isempty(circMotifs), circMotifs = 1:K; end

% Only Lg of the nfft inverse-transform outputs are ever used, so the
% inverse FFT is replaced by an [Lg x nfft] DFT matrix multiply. BLAS gemm
% on a skinny matrix beats a full ifft when Lg << nfft (81 vs 256 here).
nfft = 2^nextpow2(2*nTime - 1);
DFT  = buildLagDFT(nfft, maxLag);

for si2 = 1:numel(streams)
    st = streams{si2};
    trIdx = find(masks{si2});
    nTr = numel(trIdx);
    S.meta.(['nTrials_' st]) = nTr;
    S.meta.(['skipped_' st]) = nTr < opt.minCorrectTrials;

    if nTr < opt.minCorrectTrials
        vprintf(opt.verbose, '[%s] n=%d < %d -- skipped.\n', st, nTr, opt.minCorrectTrials);
        S.obs.(['XcorrMat_pixel_' st]) = nan(K, Lg);
        S.obs.(['XcorrMat_recon_' st]) = nan(K, Lg);
        S.obs.(['XcorrMat_H_' st])     = nan(K, Lg);
        continue;
    end
    vprintf(opt.verbose, '\n[%s] %d trials\n', st, nTr);

    DAs = double(DA3(:,:,trIdx));                       % Pv x nTime x nTr
    DAg = squeeze(mean(DAs, 1, 'omitnan'))';            % nTr x nTime (unweighted global)
    if opt.doPSTHSubtraction
        DAs = DAs - mean(DAs, 3, 'omitnan');
        DAg = DAg - mean(DAg, 1, 'omitnan');
    end

    Xp = nan(K, Lg); Xr = nan(K, Lg); Xh = nan(K, Lg);
    % Only the SUMMARY stats are kept, never the raw draws: the trial null
    % alone would be [K x Lg x nShuffle] = 24 x 81 x 1000, and it is not
    % needed once mean/std/z/p are computed. Accumulated per motif below.
    trMean = nan(K, Lg); trStd = nan(K, Lg); trZ = nan(K, Lg); trP = nan(K, Lg);
    ciMean = nan(K, Lg); ciStd = nan(K, Lg); ciZ = nan(K, Lg); ciP = nan(K, Lg);

    for k = 1:K
        % ---- reconstruction on native grid, then resample per trial ----
        Xk_native = reconOneMotif(Hcat, W_basis, k, pv, L);      % Pv x nFrames
        Xk = nan(Pv, nTime, nTr, 'single');
        for ii = 1:nTr
            tr = trIdx(ii);
            evt = numScalar(tbytDat(tr).(char(opt.evtField)));
            if ~isfinite(evt), continue; end
            % Window to this trial's own frames BEFORE interpolating. The
            % concatenated series is ~8600 frames; passing all of it to
            % interp1 to extract 119 points was the single largest hidden
            % cost in the previous version.
            sel = frameWindow(Tcat, evt + tint(1), evt + tint(end));
            if numel(sel) < 2, continue; end
            Xk(:,:,ii) = single(interpToGrid(Xk_native(:, sel), Tcat(sel) - evt, tint));
        end
        if opt.doPSTHSubtraction
            Xk = Xk - mean(Xk, 3, 'omitnan');
        end

        wk = Wt(:,k);

        % ---- (1) per-pixel xcorr, footprint-weighted: observed = diagonal
        %      of the cross-trial tensor; every trial-shuffle draw is a sum
        %      over one permutation's off-diagonal entries.
        Ctens = crossTrialTensor(DAs, Xk, wk, maxLag, opt.doFisherZ, opt.clipR, ...
                                 opt.useParfor, DFT, nfft);
        Xp(k,:) = diagMean(Ctens);

        % Draws are generated, summarized, and discarded -- nothing of size
        % [Lg x nShuffle] is retained.
        drawsT = zeros(Lg, opt.nShuffleTrial);
        for s = 1:opt.nShuffleTrial
            drawsT(:,s) = permMean(Ctens, randperm(nTr))';
        end
        [trMean(k,:), trStd(k,:), trZ(k,:), trP(k,:)] = summarizeDraws(Xp(k,:), drawsT);

        % ---- (2) reconstruction spatial mean vs global DA (scalars) ----
        Xr(k,:) = scalarXcorr(DAg, squeeze(sum(Xk .* wk, 1, 'omitnan'))', maxLag, opt.doFisherZ, opt.clipR);

        % ---- (3) raw H vs global DA (scalars) ----
        Hk = nan(nTr, nTime);
        for ii = 1:nTr
            evt = numScalar(tbytDat(trIdx(ii)).(char(opt.evtField)));
            if ~isfinite(evt), continue; end
            sel = frameWindow(Tcat, evt + tint(1), evt + tint(end));
            if numel(sel) < 2, continue; end
            Hk(ii,:) = interpToGrid(Hcat(k, sel), Tcat(sel) - evt, tint);
        end
        if opt.doPSTHSubtraction, Hk = Hk - mean(Hk, 1, 'omitnan'); end
        Xh(k,:) = scalarXcorr(DAg, Hk, maxLag, opt.doFisherZ, opt.clipR);

        % ---- circshift null (pixel path only, subset of motifs) ----
        if opt.doCircShuffle && ismember(k, circMotifs)
            drawsC = circShuffleNull(DAs, Xk, wk, maxLag, opt.doFisherZ, opt.clipR, ...
                opt.nShuffleCirc, opt.useParfor, DFT, nfft);     % Lg x nShuffleCirc
            [ciMean(k,:), ciStd(k,:), ciZ(k,:), ciP(k,:)] = summarizeDraws(Xp(k,:), drawsC);
        end

        if opt.verbose && mod(k, 4) == 0
            fprintf('  motif %d/%d (%.1f min elapsed)\n', k, K, toc(tAll)/60);
        end
    end

    S.obs.(['XcorrMat_pixel_' st]) = Xp;
    S.obs.(['XcorrMat_recon_' st]) = Xr;
    S.obs.(['XcorrMat_H_' st])     = Xh;

    S.shuf.(['mean_pixel_' st '_trialShuffle']) = trMean;
    S.shuf.(['std_pixel_'  st '_trialShuffle']) = trStd;
    S.shuf.(['z_pixel_'    st '_trialShuffle']) = trZ;
    S.shuf.(['p_pixel_'    st '_trialShuffle']) = trP;

    if opt.doCircShuffle
        S.shuf.(['mean_pixel_' st '_circShuffle']) = ciMean;
        S.shuf.(['std_pixel_'  st '_circShuffle']) = ciStd;
        S.shuf.(['z_pixel_'    st '_circShuffle']) = ciZ;
        S.shuf.(['p_pixel_'    st '_circShuffle']) = ciP;
        S.meta.circMotifs = circMotifs;
    end
end

S.meta.elapsedSec = toc(tAll);

%% -------------------- save --------------------
% Only summary stats are in S -- the raw shuffle draws were summarized and
% discarded per motif, so this file is a few hundred KB rather than tens of MB.
if strlength(strtrim(string(opt.saveDir))) > 0
    outDir = char(string(opt.saveDir));
    if exist(outDir, 'dir') ~= 7, mkdir(outDir); end

    [~, daName] = fileparts(fileDA);
    header = regexprep(daName, '_green_dff_smCollect$', '');

    % pixelSelect and the shuffle counts are in the name: two runs differing
    % only in weighting must not overwrite each other.
    outName = sprintf('%s_%s_%s_nT%d_nC%d_%s.mat', header, ...
        char(string(opt.saveKeyword)), char(pixelSelect), ...
        opt.nShuffleTrial, opt.nShuffleCirc * opt.doCircShuffle, ...
        char(datetime('today','Format','MMddyy')));
    outPath = fullfile(outDir, outName);

    S.meta.savedTo = outPath;
    result = S;   %#ok<NASGU>   'result' matches the array-task convention
    save(outPath, 'result', '-v7.3');
    vprintf(opt.verbose, 'Saved:\n  %s\n', outPath);
end

vprintf(opt.verbose, '\nDone in %.1f min.\n', S.meta.elapsedSec/60);
end

%% ========================================================================
function Xk = reconOneMotif(Hcat, W_basis, k, pv, L)
% X_k(p,t) = sum_l W(p,k,l) H_k(t-l), on H's NATIVE grid (W's lag unit).
% Edge effect: the first L-1 frames of the concatenated series are
% incomplete. Left uncorrected -- they are a few frames out of thousands
% and only affect trials at the very start of a block.
T = size(Hcat, 2);
Hlag = zeros(L, T);
for l = 0:L-1
    Hlag(l+1, l+1:T) = Hcat(k, 1:T-l);
end
Xk = squeeze(W_basis(pv, k, :)) * Hlag;    % Pv x T
end

%% ========================================================================
function Y = interpToGrid(M, tRel, tint)
% Linear interpolation of [nRow x nSample] onto tint, no extrapolation.
Y = nan(size(M,1), numel(tint));
ok = isfinite(tRel);
if sum(ok) < 2, return; end
t = tRel(ok); Msub = M(:, ok);
[t, si] = sort(t); Msub = Msub(:, si);
[t, ui] = unique(t, 'stable'); Msub = Msub(:, ui);
if numel(t) < 2, return; end
inR = tint >= t(1) & tint <= t(end);
if ~any(inR), return; end
Y(:, inR) = interp1(t(:), Msub', tint(inR)', 'linear')';
end

%% ========================================================================
function C = crossTrialTensor(DAs, Xk, wk, maxLag, doFisherZ, clipR, useParfor, DFT, nfft)
% C(a,b,:) = footprint-weighted mean over pixels of the Fisher-z lag
% correlation between DA of trial a and motif of trial b.
%
% The forward FFTs are hoisted OUT of the pair loop: F_DA(a) depends only
% on a and F_X(b) only on b, so 2*N transforms suffice for all N^2 pairs.
% The previous version recomputed both inside every pair -- 2*N^2 forward
% transforms, the dominant cost of the whole function.
[Pv, ~, nTr] = size(DAs);

FDA  = complex(zeros(Pv, nfft, nTr, 'single'));
FX   = complex(zeros(Pv, nfft, nTr, 'single'));
nrmA = zeros(Pv, nTr, 'single');
nrmB = zeros(Pv, nTr, 'single');

for n = 1:nTr
    a = DAs(:,:,n); a(~isfinite(a)) = 0; a = a - mean(a, 2);
    b = Xk(:,:,n);  b(~isfinite(b)) = 0; b = b - mean(b, 2);
    FDA(:,:,n)  = fft(single(a), nfft, 2);
    FX(:,:,n)   = fft(single(b), nfft, 2);
    nrmA(:,n)   = sqrt(sum(a.^2, 2));
    nrmB(:,n)   = sqrt(sum(b.^2, 2));
end

Lg = 2*maxLag + 1;
C  = nan(nTr, nTr, Lg);

if useParfor
    parfor a = 1:nTr
        C(a,:,:) = tensorRow(FDA(:,:,a), nrmA(:,a), FX, nrmB, wk, DFT, doFisherZ, clipR);
    end
else
    for a = 1:nTr
        C(a,:,:) = tensorRow(FDA(:,:,a), nrmA(:,a), FX, nrmB, wk, DFT, doFisherZ, clipR);
    end
end
end

%% ========================================================================
function R = tensorRow(Fa, na, FX, nrmB, wk, DFT, doFisherZ, clipR)
nTr = size(FX, 3);
Lg  = size(DFT, 1);
R   = nan(1, nTr, Lg);
for b = 1:nTr
    r = lagCorrFromFFT(Fa, FX(:,:,b), na, nrmB(:,b), DFT);
    r = max(min(r, clipR), -clipR);
    if doFisherZ, r = atanh(r); end
    R(1,b,:) = wk' * double(r);
end
end

%% ========================================================================
function r = lagCorrFromFFT(FA, FB, nrmA, nrmB, DFT)
% Normalized lag correlation from precomputed transforms.
% Positive lag = B leads A, matching xcorr(A,B) directly (see buildLagDFT).
% Silent pixels (zero norm) return 0 rather than NaN.
S = FA .* conj(FB);                 % Pv x nfft
cc = real(S * DFT.');               % Pv x Lg -- only the lags we need
nrm = nrmA .* nrmB;
nrm(nrm == 0) = eps;
r = cc ./ nrm;
end

%% ========================================================================
function DFT = buildLagDFT(nfft, maxLag)
% Rows of the inverse-DFT matrix for exactly the lags -maxLag:maxLag.
% Circular index of lag tau in an ifft output is mod(tau, nfft).
%
% NO FLIP. real(ifft(FA .* conj(FB))) is already xcorr(A,B) in MATLAB's
% convention: c(tau) = sum_t A(t+tau) B(t), so positive tau means the
% SECOND argument leads the first. With A = DA and B = motif, positive lag
% therefore means motif leads DA -- the convention used throughout this
% project. An earlier version applied fliplr here and mirrored every
% correlogram; verified against a known ground-truth shift, not against a
% reversed reference.
lags = (-maxLag:maxLag);
idx  = mod(lags, nfft);                               % 0-based circular index
n    = 0:nfft-1;
DFT  = exp(2i*pi*(idx(:)*n)/nfft) / nfft;             % Lg x nfft
end

%% ========================================================================
function sel = frameWindow(Tcat, tLo, tHi)
% Indices of Tcat covering [tLo tHi], padded one frame each side so
% interp1 can bracket the first and last bins.
i0 = find(Tcat <= tLo, 1, 'last');  if isempty(i0), i0 = 1; end
i1 = find(Tcat >= tHi, 1, 'first'); if isempty(i1), i1 = numel(Tcat); end
sel = i0:i1;
end

%% ========================================================================
function v = diagMean(C)
nTr = size(C,1);
idx = sub2ind([nTr nTr], (1:nTr)', (1:nTr)');
Cr = reshape(C, nTr*nTr, []);
v = tanh(mean(Cr(idx,:), 1, 'omitnan'));
end

function v = permMean(C, perm)
nTr = size(C,1);
idx = sub2ind([nTr nTr], perm(:), (1:nTr)');
Cr = reshape(C, nTr*nTr, []);
v = tanh(mean(Cr(idx,:), 1, 'omitnan'));
end

%% ========================================================================
function D = circShuffleNull(DAs, Xk, wk, maxLag, doFisherZ, clipR, nShuf, useParfor, DFT, nfft)
% Independent within-trial circular shift of each signal, one shift per
% signal per trial SHARED ACROSS PIXELS, so DA's spatial pattern survives
% and only the DA-vs-motif timing is broken.
%
% Unlike the trial shuffle this cannot be cached -- each draw is a new
% signal -- but the shift is applied in the FREQUENCY domain as a phase
% ramp, so the forward transforms are still hoisted and each draw costs
% only the phase multiply plus the DFT reduction.
[Pv, nTime, nTr] = size(DAs);
Lg = 2*maxLag + 1;

FDA  = complex(zeros(Pv, nfft, nTr, 'single'));
FX   = complex(zeros(Pv, nfft, nTr, 'single'));
nrmA = zeros(Pv, nTr, 'single');
nrmB = zeros(Pv, nTr, 'single');
for n = 1:nTr
    a = DAs(:,:,n); a(~isfinite(a)) = 0; a = a - mean(a, 2);
    b = Xk(:,:,n);  b(~isfinite(b)) = 0; b = b - mean(b, 2);
    FDA(:,:,n) = fft(single(a), nfft, 2);
    FX(:,:,n)  = fft(single(b), nfft, 2);
    nrmA(:,n)  = sqrt(sum(a.^2, 2));
    nrmB(:,n)  = sqrt(sum(b.^2, 2));
end

kvec = single(0:nfft-1);
D = zeros(Lg, nShuf);
if useParfor
    parfor s = 1:nShuf
        D(:,s) = oneCircDraw(FDA, FX, nrmA, nrmB, wk, DFT, doFisherZ, clipR, nTime, nTr, kvec, nfft);
    end
else
    for s = 1:nShuf
        D(:,s) = oneCircDraw(FDA, FX, nrmA, nrmB, wk, DFT, doFisherZ, clipR, nTime, nTr, kvec, nfft);
    end
end
end

%% ========================================================================
function v = oneCircDraw(FDA, FX, nrmA, nrmB, wk, DFT, doFisherZ, clipR, nTime, nTr, kvec, nfft)
% A circular shift by m is a multiplication by exp(-2i*pi*k*m/nfft) in the
% frequency domain. NOTE the shift is circular over the nfft-padded length,
% which for the ZERO-PADDED region is equivalent to shifting within the
% padded frame -- acceptable here because the null only needs to destroy
% DA-motif alignment, not preserve exact within-trial wraparound.
Lg  = size(DFT, 1);
acc = zeros(1, Lg); n = 0;
for b = 1:nTr
    mA = randi(nTime) - 1;
    mB = randi(nTime) - 1;
    phA = exp(-2i*pi*kvec*mA/nfft);
    phB = exp(-2i*pi*kvec*mB/nfft);
    Fa = FDA(:,:,b) .* phA;
    Fb = FX(:,:,b)  .* phB;
    r = lagCorrFromFFT(Fa, Fb, nrmA(:,b), nrmB(:,b), DFT);
    r = max(min(r, clipR), -clipR);
    if doFisherZ, r = atanh(r); end
    acc = acc + wk' * double(r);
    n = n + 1;
end
v = tanh(acc / max(n,1))';
end

%% ========================================================================
function v = scalarXcorr(A, B, maxLag, doFisherZ, clipR)
% Per-trial xcorr of two [nTr x nTime] matrices, Fisher-z averaged.
nTr = size(A,1);
nfft = 2^nextpow2(2*size(A,2) - 1);
DFT = buildLagDFT(nfft, maxLag);
acc = zeros(1, 2*maxLag+1); n = 0;
for i = 1:nTr
    a = A(i,:); b = B(i,:);
    a(~isfinite(a)) = 0; b(~isfinite(b)) = 0;
    a = a - mean(a); b = b - mean(b);
    if all(a == 0) || all(b == 0), continue; end
    r = lagCorrFromFFT(fft(a, nfft), fft(b, nfft), norm(a), norm(b), DFT);
    r = max(min(r, clipR), -clipR);
    if doFisherZ, r = atanh(r); end
    acc = acc + r; n = n + 1;
end
if n == 0, v = nan(1, 2*maxLag+1); return; end
v = acc / n;
if doFisherZ, v = tanh(v); end
end

%% ========================================================================
function [mu, sd, z, pv] = summarizeDraws(obs, draws)
% draws : [Lg x nShuffle]. Returns row vectors; the draws are discarded by
% the caller, so this is the only trace the null leaves in the output.
mu = mean(draws, 2, 'omitnan')';
sd = std(draws, 0, 2, 'omitnan')'; sd(sd == 0) = eps;
z  = (obs - mu) ./ sd;
nS = size(draws, 2);
pv = (1 + sum(abs(draws' - mu) >= abs(obs - mu), 1, 'omitnan')) ./ (nS + 1);
end

%% ========================================================================
function trI = trialTypeInfoFromTbytDat(tbytDat)
% Same logic as trialTypeInfoAuditoryGngExactPath, but taking the struct
% directly instead of a path -- so trI is derived from the tbytDat already
% loaded here rather than re-read from disk or depended on as a previously
% appended variable.
trI.waterI   = cellfun(@(a) ~isempty(a), {tbytDat.water});
trI.lickI    = cellfun(@(a) ~isempty(a), {tbytDat.Lick});
trI.airpuffI = cellfun(@(a) ~isempty(a), {tbytDat.airpuff});
trI.goI      = [tbytDat.rewardTrI]' == 1;
trI.nogoI    = [tbytDat.punishTrI]' == 1;
trI.hitI     = cellfun(@(a) ~isempty(a), {tbytDat.hitLicks})' & trI.waterI';
trI.missI    = cell2mat({tbytDat.rewardTrI})' & cellfun(@(a) isempty(a), {tbytDat.water})';
trI.crI      = cell2mat({tbytDat.punishTrI})' & cellfun(@(a) isempty(a), {tbytDat.airpuff})';
if isfield(tbytDat, 'faLicks')
    trI.faI  = cellfun(@(a) ~isempty(a), {tbytDat.faLicks})' & trI.airpuffI';
else
    trI.faI  = cell2mat({tbytDat.punishTrI})' & trI.airpuffI';
end
end

%% ========================================================================
function v = numVector(x)
if iscell(x), x = cell2mat(x); end
if isempty(x), v = []; else, v = double(x(:)'); end
end

%% ========================================================================
function x = numScalar(x)
if iscell(x), x = cell2mat(x); end
if isempty(x), x = NaN; else, x = double(x(1)); end
end

function vprintf(tf, varargin)
if tf, fprintf(varargin{:}); end
end