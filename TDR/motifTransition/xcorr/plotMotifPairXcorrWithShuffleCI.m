function hFig = plotMotifPairXcorrWithShuffleCI(S, iMotif, jMotif, varargin)
% PLOTMOTIFPAIRXCORRWITHSHUFFLECI
%   Plot observed motif-pair xcorr curve across full lags (e.g., -1..1 s)
%   with bin-by-bin shuffle-based confidence interval (upper/lower bounds).
%
%   hFig = plotMotifPairXcorrWithShuffleCI(S, iMotif, jMotif, ...)
%
% INPUTS
%   S      : struct output from motifH_perTrial_crosscorr_posLag_withinTrialTimeShuffle_func
%   iMotif : source motif index (row)
%   jMotif : target motif index (col)
%
% NAME–VALUE PAIRS
%   'condition'  : 'all' (default) | 'go' | 'nogo'
%   'alpha'      : CI level; default 0.05 for 95% CI
%   'useZ'       : true (default) -> use z*std for bounds (fast)
%   'lineWidth'  : default 2
%   'doGrid'     : default true
%   'titleStr'   : default '' (auto if empty)
%
%   'figSaveDir' : base directory for saving
%   'header'     : session header string for subfolder + filename (e.g., "m1045_122424")
%   'reprint'    : if file exists, re-save anyway (default false)
%
% OUTPUT
%   hFig : figure handle

%% ---- parse ----
p = inputParser;
p.addParameter('condition', 'all', @(s) ischar(s) || isstring(s));
p.addParameter('alpha', 0.05, @(x) isnumeric(x) && isscalar(x) && x>0 && x<1);
p.addParameter('useZ', true, @(x) islogical(x) && isscalar(x));
p.addParameter('lineWidth', 2, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('doGrid', true, @(x) islogical(x) && isscalar(x));
p.addParameter('titleStr', '', @(s) ischar(s) || isstring(s));

p.addParameter('figSaveDir', "Z:\Rodent Data\dualImaging_parkj\collectFigure\motifTDR\xcorr_mds\motif_pair_xcorrelograms", ...
    @(s) ischar(s) || isstring(s));
p.addParameter('header', "", @(s) ischar(s) || isstring(s));
p.addParameter('reprint', false, @(x) islogical(x) && isscalar(x));

p.parse(varargin{:});

cond       = lower(string(p.Results.condition));
alpha      = p.Results.alpha;
useZ       = p.Results.useZ;
lw         = p.Results.lineWidth;
doGrid     = p.Results.doGrid;
titleStr   = string(p.Results.titleStr);

figSaveDir = string(p.Results.figSaveDir);
header     = string(p.Results.header);
reprint    = p.Results.reprint;

if ~ismember(cond, ["all","go","nogo"])
    error('condition must be ''all'', ''go'', or ''nogo''.');
end

% z-value for two-sided CI
zcrit = abs(norminv(alpha/2, 0, 1));  % e.g. 1.96 for 95%

%% ---- pull observed curve + shuffle curve stats ----
lagsSec = S.params.lags .* S.params.stepSec;   % seconds

switch cond
    case "all"
        obsCurve  = squeeze(S.obs.XcorrMat_sess(iMotif, jMotif, :));
        muCurve   = squeeze(S.shuf.curve.mean_all(iMotif, jMotif, :));
        sdCurve   = squeeze(S.shuf.curve.std_all(iMotif, jMotif, :));
        condLabel = "all";
    case "go"
        obsCurve  = squeeze(S.obs.XcorrMat_sess_go(iMotif, jMotif, :));
        muCurve   = squeeze(S.shuf.curve.mean_go(iMotif, jMotif, :));
        sdCurve   = squeeze(S.shuf.curve.std_go(iMotif, jMotif, :));
        condLabel = "go";
    case "nogo"
        obsCurve  = squeeze(S.obs.XcorrMat_sess_nogo(iMotif, jMotif, :));
        muCurve   = squeeze(S.shuf.curve.mean_nogo(iMotif, jMotif, :));
        sdCurve   = squeeze(S.shuf.curve.std_nogo(iMotif, jMotif, :));
        condLabel = "nogo";
end

% bounds: mean ± z*std (per-lag)
lo = muCurve - zcrit .* sdCurve;
hi = muCurve + zcrit .* sdCurve;

%% ---- plot ----
hFig = figure('Color','w'); hold on;

% CI patch
x = lagsSec(:);
yLo = lo(:);
yHi = hi(:);
patch([x; flipud(x)], [yLo; flipud(yHi)], 1, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.25);

% shuffle mean
plot(lagsSec, muCurve, '-', 'LineWidth', 1);

% observed
plot(lagsSec, obsCurve, '-', 'LineWidth', lw);

% aesthetics
xline(0, '--');
yline(0, '--');

xlabel('Lag (s)');
ylabel('xcorr (coeff)');

% title (derived)
ciPct = round(100 * (1 - alpha));
titleStr = sprintf('%s | Motif %d → %d | observed xcorr with shuffle %d%% CI', ...
    upper(condLabel), iMotif, jMotif, ciPct);
title(titleStr, 'Interpreter','none');

if doGrid, grid on; end
axis tight;

legend({'Shuffle CI','Shuffle mean','Observed','lag=0','xcorr=0'}, 'Location','best');
hold off;

%% ---- auto-save ----
if strlength(strtrim(header)) > 0
    header = strtrim(header);

    % subfolder = figSaveDir/header
    subDir = fullfile(char(figSaveDir), char(header));
    if ~exist(subDir, 'dir')
        mkdir(subDir);
    end

    % filename now includes condition to avoid overwrite:
    % pairwise_xcorr_<header>_<cond>_<iMotif>_<jMotif>.pdf
    baseName = sprintf('pairwise_xcorr_%s_%s_%d_%d', char(header), char(condLabel), iMotif, jMotif);
    outPath  = fullfile(subDir, [baseName '.pdf']);

    if exist(outPath, 'file') && ~reprint
        % do nothing
    else
        print(hFig, outPath, '-dpdf', '-painters', '-bestfit');
    end
end

end
