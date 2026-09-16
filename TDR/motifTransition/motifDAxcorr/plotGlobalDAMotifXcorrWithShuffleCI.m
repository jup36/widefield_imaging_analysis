function hFig = plotGlobalDAMotifXcorrWithShuffleCI(S, motifK, varargin)
%PLOTGLOBALDAMOTIFXCORRWITHSHUFFLECI
%   Plot the observed global-DA vs. motif-H cross-correlogram for ONE
%   motif, with shuffle-based confidence bands drawn for BOTH nulls
%   stored by globalDA_motifH_perTrial_xcorr_func:
%
%     within-trial circshift null  -> S.shuf.XcorrMat_<stream>
%     across-trial permutation null-> S.shuf.XcorrMat_<stream>_trialShuffle
%
%   Showing both matters here: they answer different questions and can
%   disagree. The within-trial null asks whether the LAG STRUCTURE beats
%   chance given each signal's own autocorrelation; the trial-shuffle null
%   asks whether the DA-motif TRIAL PAIRING matters at all. A peak that
%   clears the first but not the second reflects consistent within-trial
%   timing without trial-specific coupling, and vice versa.
%
%   DIRECTION: positive lag = H LEADS DA (see the analysis function's
%   header). The x-label states this so a saved figure is unambiguous.
%
%   hFig = plotGlobalDAMotifXcorrWithShuffleCI(S, motifK, ...)
%
% INPUTS
%   S      : output struct from globalDA_motifH_perTrial_xcorr_func
%   motifK : which motif row of the [K x L] correlogram to plot
%
% NAME-VALUE PAIRS
%   'stream'      : 'cr' (default) | 'hit' -- which trial-type stream
%   'whichNull'   : 'both' (default) | 'withinTrial' | 'trialShuffle' | 'none'
%   'alpha'       : CI level, default 0.05 -> 95%
%   'useZ'        : true (default) -> bounds are mean +/- z*std of the null.
%                   false -> empirical percentiles of the null draws, which
%                   makes no normality assumption and is preferable when
%                   the null is skewed; costs nothing since the draws are
%                   stored.
%   'markSig'     : true (default) -- shade lags where the stored per-lag
%                   permutation p is below alpha (uses S.shuf.p_*).
%   'lineWidth'   : default 2
%   'doGrid'      : default true
%   'titleStr'    : '' (auto)
%   'ax'          : axes to draw into; default new figure
%   'figSaveDir'  : base directory for saving
%   'header'      : session header for subfolder + filename (e.g. "m1045_122424")
%   'reprint'     : re-save if the file exists (default false)
%
% OUTPUT
%   hFig : figure handle

%% ---- parse ----
p = inputParser;
p.addParameter('stream', 'cr', @(s) any(strcmpi(string(s), ["hit","cr"])));
p.addParameter('whichNull', 'both', @(s) any(strcmpi(string(s), ["both","withintrial","trialshuffle","none"])));
p.addParameter('alpha', 0.05, @(x) isnumeric(x) && isscalar(x) && x>0 && x<1);
p.addParameter('useZ', true, @(x) islogical(x) && isscalar(x));
p.addParameter('markSig', true, @(x) islogical(x) && isscalar(x));
p.addParameter('lineWidth', 2, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('doGrid', true, @(x) islogical(x) && isscalar(x));
p.addParameter('titleStr', '', @(s) ischar(s) || isstring(s));
p.addParameter('ax', [], @(h) isempty(h) || isgraphics(h, 'axes'));
p.addParameter('figSaveDir', "Z:\Rodent Data\dualImaging_parkj\collectFigure\motifTDR\globalDA_motif_xcorrelograms", ...
    @(s) ischar(s) || isstring(s));
p.addParameter('header', "", @(s) ischar(s) || isstring(s));
p.addParameter('reprint', false, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});

stream    = lower(char(string(p.Results.stream)));
whichNull = lower(string(p.Results.whichNull));
alpha     = p.Results.alpha;
useZ      = p.Results.useZ;
markSig   = p.Results.markSig;
lw        = p.Results.lineWidth;

%% ---- validate + pull ----
obsField = ['XcorrMat_' stream];
assert(isfield(S, 'obs') && isfield(S.obs, obsField), ...
    'S.obs.%s not found -- was this stream computed?', obsField);

XobsAll = S.obs.(obsField);                 % K x L
K = size(XobsAll, 1);
assert(motifK >= 1 && motifK <= K, 'motifK must be in 1..%d (K = %d).', K, K);

if isfield(S.meta, ['skipped_' stream]) && S.meta.(['skipped_' stream])
    warning('Stream "%s" was skipped (n=%d < minCorrectTrials) -- the curve is all NaN.', ...
        stream, S.meta.(['nTrials_' stream]));
end

lagsSec = S.meta.lagSec(:);
obsCurve = XobsAll(motifK, :)';

% Null specs: {suffix, display name, colour}
nullSpec = {};
switch whichNull
    case "both"
        nullSpec = {{'',              'within-trial shuffle', [0.35 0.35 0.35]}, ...
                    {'_trialShuffle', 'trial shuffle',        [0.20 0.50 0.75]}};
    case "withintrial"
        nullSpec = {{'',              'within-trial shuffle', [0.35 0.35 0.35]}};
    case "trialshuffle"
        nullSpec = {{'_trialShuffle', 'trial shuffle',        [0.20 0.50 0.75]}};
    case "none"
        nullSpec = {};
end

% Drop any null that wasn't actually computed, with a note rather than an error.
keep = true(1, numel(nullSpec));
for i = 1:numel(nullSpec)
    f = ['XcorrMat_' stream nullSpec{i}{1}];
    if ~isfield(S, 'shuf') || ~isfield(S.shuf, f)
        fprintf('  (%s not present in S.shuf -- skipping that band.)\n', nullSpec{i}{2});
        keep(i) = false;
    end
end
nullSpec = nullSpec(keep);

zcrit = abs(norminv(alpha/2, 0, 1));
pctl  = 100 * [alpha/2, 1 - alpha/2];

%% ---- plot ----
if isempty(p.Results.ax)
    hFig = figure('Color','w');
    ax = axes(hFig);
else
    ax = p.Results.ax;
    hFig = ancestor(ax, 'figure');
end
hold(ax, 'on');

x = lagsSec;
legH = gobjects(0); legL = {};

for i = 1:numel(nullSpec)
    sfx = nullSpec{i}{1}; nm = nullSpec{i}{2}; col = nullSpec{i}{3};
    nullStack = S.shuf.(['XcorrMat_' stream sfx]);   % K x L x nShuffle
    draws = squeeze(nullStack(motifK, :, :));        % L x nShuffle

    if useZ
        mu = mean(draws, 2, 'omitnan');
        sd = std(draws, 0, 2, 'omitnan');
        lo = mu - zcrit .* sd;
        hi = mu + zcrit .* sd;
    else
        % Empirical percentiles -- no normality assumption. Preferred when
        % the null is skewed, which correlation nulls often are near +/-1.
        q  = prctile(draws, pctl, 2);
        lo = q(:,1); hi = q(:,2);
        mu = median(draws, 2, 'omitnan');
    end

    hP = patch(ax, [x; flipud(x)], [lo(:); flipud(hi(:))], col, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.20);
    hM = plot(ax, x, mu, '-', 'Color', col, 'LineWidth', 1.1);

    legH(end+1) = hP; %#ok<AGROW>
    legL{end+1} = sprintf('%s %d%% band', nm, round(100*(1-alpha))); %#ok<AGROW>
    legH(end+1) = hM; %#ok<AGROW>
    legL{end+1} = sprintf('%s mean', nm); %#ok<AGROW>

    % significance shading from the STORED per-lag p (computed against the
    % full null in the analysis function, not re-derived here)
    pField = ['p_' stream sfx];
    if markSig && isfield(S.shuf, pField)
        pv = S.shuf.(pField)(motifK, :);
        sigLags = x(pv(:) < alpha);
        if ~isempty(sigLags)
            yl = ylim(ax);
            yMark = yl(1) + (0.03 + 0.03*(i-1)) * diff(yl);
            plot(ax, sigLags, yMark * ones(size(sigLags)), '.', ...
                'Color', col, 'MarkerSize', 10, 'HandleVisibility', 'off');
        end
    end
end

hO = plot(ax, x, obsCurve, '-', 'Color', [0.85 0.20 0.10], 'LineWidth', lw);
legH(end+1) = hO; legL{end+1} = 'observed';

xline(ax, 0, '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
yline(ax, 0, '-',  'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off');

xlabel(ax, 'Lag (s)   \leftarrow DA leads      H leads \rightarrow');
ylabel(ax, 'DA-motif xcorr (coeff)');

if strlength(string(p.Results.titleStr)) > 0
    title(ax, char(p.Results.titleStr), 'Interpreter', 'none');
else
    psthTag = 'PSTH-subtracted';
    if isfield(S.params, 'doPSTHSubtraction') && ~S.params.doPSTHSubtraction
        psthTag = 'raw (no PSTH subtraction)';
    end
    nTr = NaN;
    if isfield(S.meta, ['nTrials_' stream]), nTr = S.meta.(['nTrials_' stream]); end
    title(ax, sprintf('Global DA vs. motif %d | %s trials (n=%d) | %s', ...
        motifK, upper(stream), nTr, psthTag), 'Interpreter', 'none');
end

if p.Results.doGrid, grid(ax, 'on'); end
box(ax, 'off'); set(ax, 'TickDir', 'out');
xlim(ax, [x(1) x(end)]);
legend(ax, legH, legL, 'Location', 'best');
hold(ax, 'off');

%% ---- auto-save ----
header = strtrim(string(p.Results.header));
if strlength(header) > 0
    subDir = fullfile(char(string(p.Results.figSaveDir)), char(header));
    if ~exist(subDir, 'dir'), mkdir(subDir); end
    % stream and PSTH state in the filename: a subtracted and an
    % unsubtracted run of the same motif must not overwrite each other.
    psthTagFile = 'psthSub';
    if isfield(S.params, 'doPSTHSubtraction') && ~S.params.doPSTHSubtraction
        psthTagFile = 'raw';
    end
    baseName = sprintf('globalDA_motif%02d_%s_%s_%s', motifK, stream, psthTagFile, char(header));
    outPath  = fullfile(subDir, [baseName '.pdf']);
    if ~(exist(outPath, 'file') && ~p.Results.reprint)
        set(hFig, 'InvertHardcopy', 'off');
        print(hFig, outPath, '-dpdf', '-painters', '-bestfit');
        fprintf('Saved figure:\n  %s\n', outPath);
    end
end
end