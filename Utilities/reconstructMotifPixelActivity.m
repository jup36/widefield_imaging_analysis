function out = reconstructMotifPixelActivity(hC, W_basis, nanpxs, varargin)
%RECONSTRUCTMOTIFPIXELACTIVITY
%   Reconstruct per-motif pixel-space activity from a seqNMF factorization
%   (Mackevicius et al., equation 1):
%
%       X_k(p,t) = sum_{l=0}^{L-1} W(p,k,l) * H_k(t-l)
%
%   and return the spatial mean over valid pixels, Xmean [K x T], which is
%   the quantity intended for cross-correlation against the global DA
%   transient.
%
%   WHY THE FULL TENSOR IS OPTIONAL
%   -------------------------------
%   Spatial averaging COMMUTES with the convolution: averaging X_k over
%   pixels equals convolving H_k with the pixel-averaged kernel
%   wbar_k(l) = mean_p W(p,k,l). So Xmean can be computed from a [K x L]
%   kernel and H alone -- no [P x K x T] array is ever formed. That matters
%   at scale: P=4096, K=24, T=860 is 676 MB per block in double, and
%   trial-aligned (T ~ 177,000) would be ~70 GB.
%
%   'returnFull' therefore defaults to false. Set it true for ONE block to
%   verify the reconstruction visually; the function then computes Xmean
%   BOTH ways and asserts they agree to numerical precision, which is the
%   check that licenses using the cheap path everywhere else.
%
%   WHAT THIS DOES AND DOES NOT FIX
%   -------------------------------
%   Reconstruction puts the motif signal into dF/F-like units with
%   calcium-like temporal kinetics (W carries the temporal profile), which
%   is a real improvement over correlating a sparse activation coefficient
%   against a continuous dF/F trace. But note that Xmean is a LINEAR FIR
%   filter of H_k -- it smooths and re-times, it cannot add information.
%   Expect the H-vs-DA and Xmean-vs-DA correlograms to differ mainly in
%   peak lag and width, not in whether a relationship exists.
%
%   If the goal is motif-SPECIFIC spatial content, averaging over all
%   pixels partly defeats it: pixels where motif k has no footprint
%   contribute zeros and dilute the trace. 'pixelSelect' offers a
%   footprint-weighted or top-fraction alternative. Averaging DA over the
%   SAME pixels would make the comparison fully apples-to-apples.
%
%   out = reconstructMotifPixelActivity(hC, W_basis, nanpxs, ...)
%
% INPUTS
%   hC       : cell array; hC{hRow, block} = [K x T] H, hC{tRow, block} =
%              [1 x T] timestamps.
%   W_basis  : [P x K x L] motif basis in FULL pixel space (P = prod(imageSize)).
%   nanpxs   : indices of invalid (NaN/vasculature/mask) pixels, excluded
%              from the spatial mean. Pass [] for none.
%
% NAME-VALUE
%   'block'       : which column of hC (default 1)
%   'hRow'        : row of hC holding H (default 2 -- time-corrected H in
%                   this project; row 1 is typically the raw fit)
%   'tRow'        : row of hC holding timestamps (default 3)
%   'motifs'      : which motifs to reconstruct (default 1:K)
%   'imageSize'   : [64 64]
%   'returnFull'  : false (default). true -> also return X [P x K x T].
%   'pixelSelect' : 'all' (default) | 'footprint' | 'topFrac'
%                     all       : unweighted mean over valid pixels
%                     footprint : mean weighted by each motif's own summed
%                                 |W| footprint (motif-specific, still linear)
%                     topFrac   : unweighted mean over the top 'topFrac'
%                                 fraction of pixels by footprint
%   'topFrac'     : 0.20 (used only by 'topFrac')
%   'doPlot'      : false. Diagnostic figure: W footprint per lag, the
%                   reconstructed frame at each motif's peak, and H vs Xmean.
%   'plotMotif'   : motif to feature in the diagnostic figure (default: first)
%
% OUTPUT (out)
%   .Xmean       [K x T]        spatial mean per motif (always)
%   .X           [P x K x T]    full reconstruction (only if returnFull)
%   .kernel      [K x L]        wbar_k(l), the pixel-averaged temporal kernel
%   .H           [K x T]        the H used
%   .timeX       [1 x T]        timestamps for those frames
%   .validPix    [P x 1] logical
%   .params      settings, for provenance

%% -------------------- parse --------------------
p = inputParser;
p.addParameter('block', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('hRow', 2, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('tRow', 3, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('motifs', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('imageSize', [64 64], @(x) isnumeric(x) && numel(x) == 2);
p.addParameter('returnFull', false, @(x) islogical(x) && isscalar(x));
p.addParameter('pixelSelect', 'all', @(s) any(strcmpi(string(s), ["all","footprint","topfrac"])));
p.addParameter('topFrac', 0.20, @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 1);
p.addParameter('doPlot', false, @(x) islogical(x) && isscalar(x));
p.addParameter('plotMotif', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
p.parse(varargin{:});
opt = p.Results;
pixelSelect = lower(string(opt.pixelSelect));

%% -------------------- pull H and W --------------------
assert(size(hC,1) >= opt.hRow, 'hC has %d rows; hRow=%d requested.', size(hC,1), opt.hRow);
assert(size(hC,2) >= opt.block, 'hC has %d blocks; block=%d requested.', size(hC,2), opt.block);

H = double(hC{opt.hRow, opt.block});
assert(~isempty(H), 'hC{%d,%d} is empty.', opt.hRow, opt.block);

[P, K, L] = size(W_basis);
assert(size(H,1) == K, 'H has %d rows but W_basis has K=%d motifs.', size(H,1), K);
T = size(H,2);

Pexp = prod(opt.imageSize);
assert(P == Pexp, ['W_basis has %d pixel rows but imageSize implies %d. ' ...
    'This function expects W_basis already in FULL pixel space.'], P, Pexp);

timeX = [];
if size(hC,1) >= opt.tRow && ~isempty(hC{opt.tRow, opt.block})
    timeX = double(hC{opt.tRow, opt.block}(:)');
    if numel(timeX) ~= T
        warning('reconstructMotif:TimeLength', ...
            'Timestamps have %d entries but H has %d frames -- truncating to the shorter.', numel(timeX), T);
        n = min(numel(timeX), T);
        timeX = timeX(1:n); H = H(:,1:n); T = n;
    end
end

motifs = opt.motifs;
if isempty(motifs), motifs = 1:K; end
motifs = unique(motifs(:)', 'stable');
assert(all(motifs >= 1 & motifs <= K), 'motifs must lie in 1..%d.', K);

%% -------------------- valid pixels --------------------
validPix = true(P,1);
if ~isempty(nanpxs)
    nanpxs = nanpxs(:);
    assert(all(nanpxs >= 1 & nanpxs <= P), 'nanpxs contains indices outside 1..%d.', P);
    validPix(nanpxs) = false;
end
assert(any(validPix), 'nanpxs excludes every pixel.');
fprintf('Pixels: %d total, %d valid (%d excluded by nanpxs).\n', P, sum(validPix), sum(~validPix));

%% -------------------- lag matrix --------------------
% Hlag(l+1, t) = H_k(t-l). Building this once turns the convolution into a
% single matrix product per motif: X_k = W_k(:,:) * Hlag  ([P x L]*[L x T]).
    function Hlag = local_lagMatrix(hRow)
        Hlag = zeros(L, T);
        for l = 0:L-1
            Hlag(l+1, l+1:T) = hRow(1:T-l);
        end
    end

%% -------------------- spatial weights per motif --------------------
% wSum(:,k) is motif k's footprint: total |W| across lags at each pixel.
wSum = squeeze(sum(abs(W_basis), 3));       % P x K
if K == 1, wSum = wSum(:); end

weights = zeros(P, K);
for k = 1:K
    w = zeros(P,1);
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
    if s > 0, weights(:,k) = w / s; else, weights(:,k) = 0; end
end

%% -------------------- reconstruct --------------------
Xmean  = nan(K, T);
kernel = nan(K, L);

if opt.returnFull
    X = nan(P, K, T, 'single');   % single: halves the footprint, ample precision here
    fprintf('Allocating full reconstruction: %d x %d x %d single (%.2f GB).\n', ...
        P, numel(motifs), T, P*numel(motifs)*T*4/1e9);
end

for k = motifs
    Hlag = local_lagMatrix(H(k,:));
    Wk   = squeeze(W_basis(:,k,:));           % P x L
    if L == 1, Wk = Wk(:); end

    % cheap path: pixel-weighted kernel, then one [1 x L]*[L x T] product
    kernel(k,:) = (weights(:,k)' * Wk);       % 1 x L
    Xmean(k,:)  = kernel(k,:) * Hlag;

    if opt.returnFull
        Xk = Wk * Hlag;                       % P x T  -- the full reconstruction
        Xk(~validPix, :) = NaN;               % mask invalid pixels for display
        X(:,k,:) = single(Xk);

        % Verify the identity that licenses the cheap path everywhere else.
        % Index to valid pixels explicitly: masked entries are NaN, and
        % 0*NaN = NaN in MATLAB, so a weighted product over all pixels
        % would be NaN regardless of the zero weight.
        vp  = validPix;
        chk = weights(vp,k)' * Xk(vp,:);
        ref = Xmean(k,:);
        tol = 1e-6 * max(1, max(abs(ref)));
        assert(max(abs(chk - ref)) < tol, ...
            'Motif %d: spatial-mean identity failed (max diff %.3g) -- check nanpxs/weights.', ...
            k, max(abs(chk - ref)));
    end
end

if opt.returnFull
    fprintf('Full reconstruction verified against the kernel path for %d motif(s).\n', numel(motifs));
end

%% -------------------- diagnostic figure --------------------
if opt.doPlot
    kPlot = opt.plotMotif;
    if isempty(kPlot), kPlot = motifs(1); end
    assert(ismember(kPlot, motifs), 'plotMotif %d was not reconstructed.', kPlot);
    local_diagnosticFigure(W_basis, H, Xmean, X_or_empty(opt.returnFull), ...
        kPlot, opt.imageSize, validPix, timeX, L);
end

    function Xarg = X_or_empty(tf)
        if tf, Xarg = X; else, Xarg = []; end
    end

%% -------------------- output --------------------
out = struct();
out.Xmean    = Xmean;
out.kernel   = kernel;
out.H        = H;
out.timeX    = timeX;
out.validPix = validPix;
out.motifs   = motifs;
out.params   = opt;
if opt.returnFull, out.X = X; end
end

%% ========================================================================
function local_diagnosticFigure(W_basis, H, Xmean, X, k, imageSize, validPix, timeX, L)
% Three checks in one figure:
%   row 1 : motif k's W footprint at each lag -- is the basis sensible?
%   row 2 : reconstructed frames around the peak of H_k -- does the
%           reconstruction reproduce that spatial pattern in time?
%   row 3 : H_k vs. Xmean_k -- the temporal relationship the cheap path
%           will actually feed into the cross-correlation.
figure('Color','w', 'Position', [100 100 1400 800]);

nShow = min(L, 10);
for l = 1:nShow
    subplot(3, nShow, l);
    img = reshape(W_basis(:,k,l), imageSize);
    img(~reshape(validPix, imageSize)) = NaN;
    imagesc(img, 'AlphaData', ~isnan(img)); axis image off;
    title(sprintf('W lag %d', l-1), 'FontSize', 8);
end

[~, tPk] = max(H(k,:));
offs = round(linspace(0, L-1, nShow));
for i = 1:nShow
    subplot(3, nShow, nShow + i);
    tt = min(size(H,2), tPk + offs(i));
    if isempty(X)
        axis off;
        if i == 1, text(0, 0.5, 'returnFull=false', 'FontSize', 9); end
        continue;
    end
    img = reshape(double(X(:,k,tt)), imageSize);
    imagesc(img, 'AlphaData', ~isnan(img)); axis image off;
    title(sprintf('t=peak+%d', offs(i)), 'FontSize', 8);
end

subplot(3, 1, 3);
if isempty(timeX), x = 1:size(H,2); xl = 'Frame'; else, x = timeX; xl = 'Time (s)'; end
yyaxis left;
plot(x, H(k,:), '-', 'LineWidth', 1.2); ylabel(sprintf('H_{%d}', k));
yyaxis right;
plot(x, Xmean(k,:), '-', 'LineWidth', 1.6); ylabel('spatial mean of reconstruction');
xlabel(xl); set(gca, 'TickDir', 'out'); grid on; box off;
title(sprintf('Motif %d: H vs. pixel-averaged reconstruction (a linear FIR filter of H)', k));
end