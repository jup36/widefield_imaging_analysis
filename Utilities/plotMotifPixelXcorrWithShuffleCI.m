function hFig = plotMotifPixelXcorrWithShuffleCI(R, motifK, varargin)
%PLOTMOTIFPIXELXCORRWITHSHUFFLECI
%   Plot one motif's pixel-space DA-motif cross-correlogram from a
%   globalDA_motifPixel_xcorr_func result, with the trial-shuffle null
%   drawn as a band and significant lags marked.
%
%   Repurposed from plotGlobalDAMotifXcorrWithShuffleCI for the pixel-space
%   output, which differs in three ways:
%
%     1. Three observed correlograms per motif (XcorrMat_pixel_*,
%        XcorrMat_recon_*, XcorrMat_H_*), not one.
%     2. Only the PIXEL variant carries a null. recon and H are computed
%        from the same trials but without shuffles, so they are drawn as
%        reference lines only.
%     3. The null is stored as SUMMARY stats (mean/std/z/p per motif x lag),
%        the raw draws having been discarded to keep files small. The band is
%        therefore mean +/- z*std -- the empirical-percentile option of the
%        old function is not possible here and is not offered.
%
%   Units: observed and null are both in correlation (r) space -- both were
%   averaged in Fisher-z and tanh'd back -- so the band and the curves share
%   a scale and can be overlaid directly.
%
%   Why overlay recon: the pixel and recon correlograms differ ONLY by the
%   spatial covariance term (correlate-then-average vs average-then-
%   correlate), so the gap between them is exactly how much DA's spatial
%   pattern matches the motif's footprint. If they overlay, DA is spatially
%   homogeneous with respect to that motif.
%
%   DIRECTION: positive lag = MOTIF LEADS DA.
%
%   hFig = plotMotifPixelXcorrWithShuffleCI(R, motifK, ...)
%
% INPUTS
%   R      : the saved struct -- either the loaded file (with field .result)
%            or result itself
%   motifK : motif row to plot (1..K)
%
% NAME-VALUE
%   'stream'      : 'hit' (default) | 'cr'
%   'showRecon'   : true  -- overlay XcorrMat_recon_* (reference, no band)
%   'showH'       : false -- overlay XcorrMat_H_* (reference, no band)
%   'alpha'       : 0.05  -- band level and significance threshold
%   'markSig'     : true  -- mark lags with stored trial-shuffle p < alpha
%   'lineWidth'   : 2
%   'yLim'        : [] (auto) or [lo hi]
%   'titleStr'    : '' (auto)
%   'ax'          : axes to draw into; default new figure
%   'figSaveDir'  : '' -- set to save a PDF
%   'header'      : '' -- session header for the filename
%   'reprint'     : false
%
% OUTPUT
%   hFig : figure handle

%% ---- parse ----
p = inputParser;
p.addParameter('stream', 'hit', @(s) any(strcmpi(string(s), ["hit","cr"])));
p.addParameter('showRecon', true, @(x) islogical(x) && isscalar(x));
p.addParameter('showH', false, @(x) islogical(x) && isscalar(x));
p.addParameter('alpha', 0.05, @(x) isnumeric(x) && isscalar(x) && x>0 && x<1);
p.addParameter('markSig', true, @(x) islogical(x) && isscalar(x));
p.addParameter('lineWidth', 2, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('yLim', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2 && x(2)>x(1)));
p.addParameter('titleStr', '', @(s) ischar(s) || isstring(s));
p.addParameter('ax', [], @(h) isempty(h) || isgraphics(h, 'axes'));
p.addParameter('figSaveDir', compatiblepath('Z:\Rodent Data\dualImaging_parkj\collectFigure\motifTDR\xcorr_motif_DA'), @(s) ischar(s) || isstring(s));
p.addParameter('header', '', @(s) ischar(s) || isstring(s));
p.addParameter('reprint', false, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;
st = lower(char(string(opt.stream)));

%% ---- accept either the loaded file or result itself ----
if isfield(R, 'result') && ~isfield(R, 'obs')
    R = R.result;
end
assert(isfield(R, 'obs') && isfield(R, 'shuf') && isfield(R, 'meta'), ...
    'Input must be a globalDA_motifPixel_xcorr_func result (fields obs, shuf, meta).');

fPix = ['XcorrMat_pixel_' st];
assert(isfield(R.obs, fPix), 'R.obs.%s not found.', fPix);

K = size(R.obs.(fPix), 1);
assert(isnumeric(motifK) && isscalar(motifK) && motifK >= 1 && motifK <= K, ...
    'motifK must be in 1..%d.', K);

if isfield(R.meta, ['skipped_' st]) && R.meta.(['skipped_' st])
    warning('Stream "%s" was skipped (n=%d < minCorrectTrials) -- curves are NaN.', ...
        st, R.meta.(['nTrials_' st]));
end

x = R.meta.lagSec(:);
yPix = R.obs.(fPix)(motifK, :)';

%% ---- null: summary stats only ----
fMu = ['mean_pixel_' st '_trialShuffle'];
fSd = ['std_pixel_'  st '_trialShuffle'];
fP  = ['p_pixel_'    st '_trialShuffle'];
hasNull = isfield(R.shuf, fMu) && isfield(R.shuf, fSd);
if ~hasNull
    warning('No trial-shuffle summaries for stream "%s" -- plotting observed only.', st);
end

zcrit = abs(norminv(opt.alpha/2, 0, 1));

%% ---- plot ----
if isempty(opt.ax)
    hFig = figure('Color', 'w');
    ax = axes(hFig);
else
    ax = opt.ax;
    hFig = ancestor(ax, 'figure');
end
hold(ax, 'on');

legH = gobjects(0); legL = {};
nullCol = [0.20 0.50 0.75];

if hasNull
    mu = R.shuf.(fMu)(motifK, :)';
    sd = R.shuf.(fSd)(motifK, :)';
    lo = mu - zcrit .* sd;
    hi = mu + zcrit .* sd;
    ok = isfinite(lo) & isfinite(hi);
    hB = patch(ax, [x(ok); flipud(x(ok))], [lo(ok); flipud(hi(ok))], nullCol, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.20);
    hM = plot(ax, x, mu, '-', 'Color', nullCol, 'LineWidth', 1.1);
    legH(end+1) = hB; legL{end+1} = sprintf('trial shuffle %d%% band', round(100*(1-opt.alpha)));
    legH(end+1) = hM; legL{end+1} = 'trial shuffle mean';
end

% references first so the pixel curve draws on top
if opt.showH && isfield(R.obs, ['XcorrMat_H_' st])
    yH = R.obs.(['XcorrMat_H_' st])(motifK, :)';
    legH(end+1) = plot(ax, x, yH, ':', 'Color', [0.93 0.69 0.13], 'LineWidth', 1.6);
    legL{end+1} = 'raw H vs global DA';
end
if opt.showRecon && isfield(R.obs, ['XcorrMat_recon_' st])
    yR = R.obs.(['XcorrMat_recon_' st])(motifK, :)';
    legH(end+1) = plot(ax, x, yR, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.6);
    legL{end+1} = 'recon spatial mean vs global DA';
end

legH(end+1) = plot(ax, x, yPix, '-', 'Color', [0.85 0.20 0.10], 'LineWidth', opt.lineWidth);
legL{end+1} = 'pixel-space (footprint-weighted)';

xline(ax, 0, '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
yline(ax, 0, '-',  'Color', [0.82 0.82 0.82], 'HandleVisibility', 'off');

if ~isempty(opt.yLim), ylim(ax, opt.yLim); end

% Significance marks from the STORED p (computed against the full null in
% the analysis function), placed after the y-limits are final so they sit
% inside the axes rather than being clipped.
if opt.markSig && hasNull && isfield(R.shuf, fP)
    pv = R.shuf.(fP)(motifK, :);
    sigLags = x(pv(:) < opt.alpha);
    if ~isempty(sigLags)
        yl = ylim(ax);
        yMark = yl(1) + 0.04 * diff(yl);
        plot(ax, sigLags, yMark * ones(size(sigLags)), 's', ...
            'MarkerSize', 4, 'MarkerFaceColor', nullCol, 'MarkerEdgeColor', 'none', ...
            'HandleVisibility', 'off');
    end
end

xlabel(ax, 'Lag (s)   \leftarrow DA leads      motif leads \rightarrow');
ylabel(ax, 'DA-motif xcorr (coeff)');

if strlength(string(opt.titleStr)) > 0
    title(ax, char(opt.titleStr), 'Interpreter', 'none');
else
    psthTag = 'PSTH-subtracted';
    if isfield(R.params, 'doPSTHSubtraction') && ~R.params.doPSTHSubtraction
        psthTag = 'raw';
    end
    nTr = NaN;
    if isfield(R.meta, ['nTrials_' st]), nTr = R.meta.(['nTrials_' st]); end
    title(ax, sprintf('Motif %d | %s trials (n=%d) | %s', ...
        motifK, upper(st), nTr, psthTag), 'Interpreter', 'none');
end

grid(ax, 'on'); box(ax, 'off'); set(ax, 'TickDir', 'out');
xlim(ax, [x(1) x(end)]);
legend(ax, legH, legL, 'Location', 'best');
hold(ax, 'off');

%% ---- optional save ----
if strlength(strtrim(string(opt.figSaveDir))) > 0
    outDir = char(string(opt.figSaveDir));
    if exist(outDir, 'dir') ~= 7, mkdir(outDir); end
    hdr = strtrim(char(string(opt.header)));
    if isempty(hdr), hdr = 'session'; end
    outFile = fullfile(outDir, sprintf('%s_pixelXcorr_motif%02d_%s.pdf', hdr, motifK, st));
    if ~(exist(outFile, 'file') && ~opt.reprint)
        set(hFig, 'InvertHardcopy', 'off');
        print(hFig, outFile, '-dpdf', '-painters', '-bestfit');
        fprintf('Saved figure:\n  %s\n', outFile);
    end
end
end