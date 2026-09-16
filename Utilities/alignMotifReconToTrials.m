function [tbytDat_XAligned, info] = alignMotifReconToTrials(hC, W_basis, nanpxs, tbytDat, varargin)
%ALIGNMOTIFRECONTOTRIALS
%   Build tbytDat_XAligned: the per-trial, interpolated, pixel-averaged
%   reconstruction of every motif, on the same time grid used for the
%   global DA alignment.
%
%   tbytDat_XAligned{1, tr} : [K x numel(tint)] motif traces for trial tr
%   tbytDat_XAligned{2, tr} : tint, seconds from the event
%
%   This matches tbytDat_hAligned in shape and grid, so it can be dropped
%   straight into globalDA_motifH_perTrial_xcorr_func in place of H.
%
%   MEMORY
%   ------
%   The [P x K x T] pixel reconstruction is NEVER formed. Spatial averaging
%   commutes with the convolution, so the pixel-averaged trace is
%
%       Xmean_k(t) = sum_l wbar_k(l) * H_k(t-l),   wbar_k(l) = mean_p W(p,k,l)
%
%   i.e. H convolved with a [K x L] kernel. Per block the working set is
%   [K x T] (24 x 860 here, ~0.2 MB); concatenated across blocks it is
%   [K x sum(T)], a few MB. That identity is the one verified by
%   reconstructMotifPixelActivity with 'returnFull', true -- run that check
%   once per new W_basis before trusting this function.
%
%   BLOCKS
%   ------
%   All blocks are concatenated on their own absolute timestamps and sorted
%   before alignment, so a trial whose window straddles a block boundary is
%   filled from both blocks rather than truncated at the seam.
%
%   [tbytDat_XAligned, info] = alignMotifReconToTrials(hC, W_basis, nanpxs, tbytDat, ...)
%
% NAME-VALUE
%   'hRow'        : row of hC holding H (default 2). Row 1 is typically the
%                   raw fit and row 2 the time-corrected H -- confirm which
%                   your file uses, since the two sit on different time bases.
%   'tRow'        : row of hC holding timestamps (default 3)
%   'blocks'      : which columns of hC to use (default: all non-empty)
%   'motifs'      : which motifs to include (default 1:K). Rows of the
%                   output are in this order; info.motifs records it.
%   'imageSize'   : [64 64]
%   'alignWin'    : [-0.9 5]
%   'alignDt'     : 0.01
%   'evtField'    : 'evtOn'
%   'pixelSelect' : 'all' (default) | 'footprint' | 'topFrac' -- see
%                   reconstructMotifPixelActivity. 'all' averages over every
%                   valid pixel, including those where a motif has no
%                   footprint, which dilutes motif-specific signal.
%   'topFrac'     : 0.20
%   'ledPrefix'   : 'lime' (default). NOT used for the interpolation --
%                   hC carries its own frame timestamps. Used only as a
%                   coverage cross-check: if tbytDat has <ledPrefix>LED,
%                   each trial's LED time span is compared against the
%                   frames actually found, and a systematic shortfall is
%                   reported. Pass '' to skip the check.
%   'verbose'     : true
%
% OUTPUT (info)
%   .tint, .motifs, .kernel [K x L], .nTrialsAligned, .fracCovered per trial,
%   .nFramesUsed per trial, .ledCheck, .params

%% -------------------- parse --------------------
p = inputParser;
p.addParameter('hRow', 2, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('tRow', 3, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('blocks', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('motifs', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('imageSize', [64 64], @(x) isnumeric(x) && numel(x) == 2);
p.addParameter('alignWin', [-0.9 5], @(x) isnumeric(x) && numel(x) == 2 && x(2) > x(1));
p.addParameter('alignDt', 0.01, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('evtField', 'evtOn', @(s) ischar(s) || isstring(s));
p.addParameter('pixelSelect', 'all', @(s) any(strcmpi(string(s), ["all","footprint","topfrac"])));
p.addParameter('topFrac', 0.20, @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 1);
p.addParameter('ledPrefix', 'lime', @(s) ischar(s) || isstring(s));
p.addParameter('verbose', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;
pixelSelect = lower(string(opt.pixelSelect));
evtField = char(string(opt.evtField));

[P, K, L] = size(W_basis);
assert(P == prod(opt.imageSize), ...
    'W_basis has %d pixel rows but imageSize implies %d.', P, prod(opt.imageSize));

motifs = opt.motifs;
if isempty(motifs), motifs = 1:K; end
motifs = unique(motifs(:)', 'stable');
assert(all(motifs >= 1 & motifs <= K), 'motifs must lie in 1..%d.', K);
nM = numel(motifs);

blocks = opt.blocks;
if isempty(blocks)
    blocks = find(~cellfun(@isempty, hC(opt.hRow, :)));
end
assert(~isempty(blocks), 'No non-empty blocks found in hC row %d.', opt.hRow);

%% -------------------- pixel weights (once, not per block) --------------------
validPix = true(P, 1);
if ~isempty(nanpxs)
    nanpxs = nanpxs(:);
    assert(all(nanpxs >= 1 & nanpxs <= P), 'nanpxs has indices outside 1..%d.', P);
    validPix(nanpxs) = false;
end
assert(any(validPix), 'nanpxs excludes every pixel.');

wSum = squeeze(sum(abs(W_basis), 3));     % P x K footprint
if K == 1, wSum = wSum(:); end

kernel = zeros(nM, L);
for ii = 1:nM
    k = motifs(ii);
    w = zeros(P, 1);
    switch pixelSelect
        case "all"
            w(validPix) = 1;
        case "footprint"
            w = wSum(:,k) .* validPix;
        case "topfrac"
            fp = wSum(:,k); fp(~validPix) = -inf;
            nTop = max(1, round(opt.topFrac * sum(validPix)));
            [~, ord] = sort(fp, 'descend');
            w(ord(1:nTop)) = 1;
    end
    s = sum(w);
    if s > 0, w = w / s; end
    kernel(ii, :) = w' * squeeze(W_basis(:, k, :));   % 1 x L
end

if opt.verbose
    fprintf('Pixels: %d total, %d valid. Kernel: [%d x %d] (pixelSelect = %s).\n', ...
        P, sum(validPix), nM, L, pixelSelect);
end

%% -------------------- per block: Xmean via the kernel, then concatenate ----
Xall = [];
Tall = [];
for b = blocks
    H = double(hC{opt.hRow, b});
    if isempty(H), continue; end
    assert(size(H,1) == K, 'Block %d: H has %d rows but W_basis has K=%d.', b, size(H,1), K);
    T = size(H, 2);

    assert(size(hC,1) >= opt.tRow && ~isempty(hC{opt.tRow, b}), ...
        'Block %d: no timestamps in hC row %d.', b, opt.tRow);
    tBlk = double(hC{opt.tRow, b}(:)');
    if numel(tBlk) ~= T
        n = min(numel(tBlk), T);
        warning('alignMotifRecon:TimeLength', ...
            'Block %d: %d timestamps vs %d frames -- truncating to %d.', b, numel(tBlk), T, n);
        tBlk = tBlk(1:n); H = H(:, 1:n); T = n;
    end

    Xb = zeros(nM, T);
    for ii = 1:nM
        hRowK = H(motifs(ii), :);
        Hlag = zeros(L, T);
        for l = 0:L-1
            Hlag(l+1, l+1:T) = hRowK(1:T-l);
        end
        Xb(ii, :) = kernel(ii, :) * Hlag;
    end

    Xall = [Xall, Xb];   %#ok<AGROW>  -- [nM x sum(T)], a few MB at most
    Tall = [Tall, tBlk]; %#ok<AGROW>
end

assert(~isempty(Tall), 'No usable frames found across the requested blocks.');

% Sort and de-duplicate once, globally: blocks may not be stored in
% chronological order, and interp1 needs strictly increasing x.
[Tall, si] = sort(Tall, 'ascend');
Xall = Xall(:, si);
[Tall, ui] = unique(Tall, 'stable');
Xall = Xall(:, ui);

if opt.verbose
    fprintf('Concatenated %d blocks -> [%d x %d] (%.1f MB), t = [%.2f %.2f] s.\n', ...
        numel(blocks), nM, numel(Tall), numel(Xall)*8/1e6, Tall(1), Tall(end));
end

%% -------------------- trial alignment --------------------
tint = opt.alignWin(1):opt.alignDt:opt.alignWin(2);
nTime = numel(tint);
nTrials = numel(tbytDat);

tbytDat_XAligned = cell(2, nTrials);
for tr = 1:nTrials
    tbytDat_XAligned{1, tr} = NaN(nM, nTime);
    tbytDat_XAligned{2, tr} = tint;
end

fracCovered = zeros(1, nTrials);
nFramesUsed = zeros(1, nTrials);

for tr = 1:nTrials
    evt = local_numScalar(tbytDat(tr).(evtField));
    if ~isfinite(evt), continue; end

    % Pull only the frames this trial needs. One extra frame of padding on
    % each side so interp1 can bracket the first and last bins rather than
    % leaving them NaN for want of a neighbour.
    lo = evt + tint(1);  hi = evt + tint(end);
    i0 = find(Tall <= lo, 1, 'last');   if isempty(i0), i0 = 1; end
    i1 = find(Tall >= hi, 1, 'first');  if isempty(i1), i1 = numel(Tall); end
    idx = i0:i1;
    if numel(idx) < 2, continue; end

    tRel = Tall(idx) - evt;
    Y = Xall(:, idx);                    % nM x nSel
    nFramesUsed(tr) = numel(idx);

    % No extrapolation: bins outside this trial's own coverage stay NaN, so
    % a short trial is visibly short rather than padded with edge values.
    inRange = tint >= tRel(1) & tint <= tRel(end);
    if ~any(inRange), continue; end

    % interp1 columnwise: [nSel x nM] -> [nIn x nM], then transpose back.
    tbytDat_XAligned{1, tr}(:, inRange) = interp1(tRel(:), Y', tint(inRange)', 'linear')';
    fracCovered(tr) = sum(inRange) / nTime;
end

nAligned = sum(cellfun(@(x) any(isfinite(x(:))), tbytDat_XAligned(1, :)));

%% -------------------- shape asserts (mirror the DA pipeline) --------------
assert(isequal(size(tbytDat_XAligned), [2 nTrials]));
for tr = 1:nTrials
    assert(isequal(size(tbytDat_XAligned{1,tr}), [nM nTime]));
    assert(isequal(tbytDat_XAligned{2,tr}, tint));
end

if opt.verbose
    fprintf('Aligned %d/%d trials to [%.2f %.2f] s at %g s (%d bins); median coverage %.0f%%.\n', ...
        nAligned, nTrials, opt.alignWin(1), opt.alignWin(2), opt.alignDt, nTime, ...
        100*median(fracCovered(fracCovered > 0)));
end

%% -------------------- optional LED coverage cross-check --------------------
% The interpolation does not use the LED times -- hC carries its own frame
% timestamps. This only asks whether the frames found for each trial span
% as much of the trial as the LED train says they should; a systematic
% shortfall means hC and tbytDat disagree about when trials happened.
ledCheck = struct('checked', false);
ledField = [char(string(opt.ledPrefix)) 'LED'];
if strlength(string(opt.ledPrefix)) > 0 && isfield(tbytDat, ledField)
    ledSpan = nan(1, nTrials); frmSpan = nan(1, nTrials);
    for tr = 1:nTrials
        tL = local_numVector(tbytDat(tr).(ledField));
        if numel(tL) < 2, continue; end
        evt = local_numScalar(tbytDat(tr).(evtField));
        if ~isfinite(evt), continue; end
        ledSpan(tr) = min(max(tL) - evt, tint(end)) - max(min(tL) - evt, tint(1));
        frmSpan(tr) = fracCovered(tr) * (tint(end) - tint(1));
    end
    ok = isfinite(ledSpan) & isfinite(frmSpan) & ledSpan > 0;
    shortfall = (ledSpan(ok) - frmSpan(ok)) ./ ledSpan(ok);
    ledCheck = struct('checked', true, 'field', ledField, ...
        'medianShortfallFrac', median(shortfall), 'nChecked', sum(ok));
    if opt.verbose
        fprintf('LED coverage check (%s): median shortfall %.1f%% over %d trials.\n', ...
            ledField, 100*median(shortfall), sum(ok));
    end
    if median(shortfall) > 0.10
        warning('alignMotifRecon:Coverage', ...
            ['Frames cover %.0f%% less of each trial than the %s train implies. ' ...
             'hC timestamps and tbytDat may disagree about trial timing.'], ...
            100*median(shortfall), ledField);
    end
elseif strlength(string(opt.ledPrefix)) > 0 && opt.verbose
    fprintf('LED coverage check skipped: tbytDat has no field "%s".\n', ledField);
end

%% -------------------- info --------------------
info = struct();
info.tint           = tint;
info.motifs         = motifs;
info.kernel         = kernel;
info.nTrialsAligned = nAligned;
info.fracCovered    = fracCovered;
info.nFramesUsed    = nFramesUsed;
info.blocksUsed     = blocks;
info.validPix       = validPix;
info.ledCheck       = ledCheck;
info.params         = opt;
end

%% ========================================================================
function v = local_numVector(x)
if iscell(x), x = cell2mat(x); end
if isempty(x), v = []; else, v = double(x(:)'); end
end

function x = local_numScalar(x)
if iscell(x), x = cell2mat(x); end
if isempty(x), x = NaN; else, x = double(x(1)); end
end