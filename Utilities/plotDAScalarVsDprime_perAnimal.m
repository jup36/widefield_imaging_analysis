function out = plotDAScalarVsDprime_perAnimal(rezC, sessT, varargin)
%PLOTDASCALARVSDPRIME_PERANIMAL
%   DA-vs-d' correlation computed WITHIN each animal across its own
%   sessions, then summarized across animals. Two figures:
%
%   FIGURE 1 -- per-animal fits: the usual 2 x 2 scatter (rows Hit/CR,
%     columns cue/response), with one regression line per animal drawn over
%     that animal's own sessions. Shows directly whether animals agree.
%
%   FIGURE 2 -- per-animal r: one dot per animal per panel, with the
%     across-animal mean and a group-level test.
%
%   WHY THIS OVER A POOLED CORRELATION
%   ----------------------------------
%   Each animal's r uses only that animal's sessions, so between-animal
%   offsets in DA or d' cannot contribute -- by construction, not by
%   correction. The animal then becomes the unit of inference, matching the
%   two-stage (sessions within animal, then animals) convention used
%   throughout this project.
%
%   GROUP TEST: mean Fisher-z of the per-animal r's, tested by EXACT
%   sign-flip permutation over all 2^nAnimals sign assignments (512 for 9
%   mice). Under H0 -- no consistent within-animal relationship -- each
%   animal's r is symmetric about zero, so flipping signs is exchangeable.
%   Exact, assumption-light, and valid at n = 9. The smallest attainable p
%   is 2/2^n (two-sided), ~0.004 for 9 animals.
%
%   PER-ANIMAL SIGNIFICANCE is reported alongside each animal's r, with two
%   cautions built in. (1) With ~5-10 sessions, an individual test has
%   little power: a true within-animal r of 0.5 is detected well under half
%   the time at n = 8. A non-significant animal is therefore NOT evidence of
%   no relationship. (2) 9 animals x 4 panels is 36 tests, so each panel also
%   carries BH-FDR q-values across its animals. The GROUP test remains the
%   primary inference; per-animal p's say which animals drive it.
%
%   Per-animal p comes from 'perAnimalTest':
%     'parametric'  (default) -- MATLAB corr()'s t-based p (Pearson) or its
%                    exact/approximate p (Spearman)
%     'permutation' -- two-sided, shuffling DA across that animal's own
%                    sessions. Assumption-free, but with n sessions the
%                    finest attainable p is 1/(nPermAnimal+1).
%
%   GROUP-LABEL RELABELING TEST ('groupPermTest', default true)
%   ------------------------------------------------------------
%   Asks whether fast and slow learners differ in their per-animal r.
%   Statistic: mean Fisher-z r of fast animals minus that of slow animals.
%   Null: every way of choosing which nFast of the n animals carry the
%   "fast" label, enumerated exhaustively -- C(7,2) = 21 for this cohort --
%   so the test is exact. Two-sided p counts partitions whose |difference|
%   is at least the observed one; the one-sided p's are also returned.
%
%   READ THE FLOOR BEFORE THE p. The smallest attainable p is
%   1/nPartitions -- 1/21 = 0.048 here -- so even the most extreme possible
%   split barely clears 0.05. The DA cohort has only 2 fast learners (m1092
%   and m1094 were imaged for calcium only, without a DA sensor), so this
%   limit is structural, not a data-processing loss. Each panel title shows
%   the floor next to the p for that reason.
%
%   And if this contrast was chosen because a panel looked extreme, its p
%   is exploratory: a test selected for fitting the data will tend to pass.
%
%   CAVEAT (same as the pooled version): d' rises with training, so within
%   an animal, d' and session order are strongly correlated. Set
%   'xVar','sessOrder' to run the identical analysis on time in task -- if
%   that relationship is as strong, time is the simpler explanation.
%
%   out = plotDAScalarVsDprime_perAnimal(rezC, sessT, ...)
%
% NAME-VALUE
%   'animals'       : {} (default, all) or a cellstr/string of mouse IDs,
%                     e.g. 'm1045' or {'m1045','m1049'}. With a SINGLE
%                     animal, its r and p go straight into Figure 1's panel
%                     titles and Figure 2 is skipped -- a group test on one
%                     animal is undefined. BH q-values are computed over
%                     the animals SELECTED, so for one animal q equals p.
%   'signal'        : 'global' (default) or motif number
%   'xVar'          : 'dprime' (default) | 'sessOrder'
%   'corrType'      : 'Pearson' (default) | 'Spearman'
%   'minSessions'   : 4   -- an animal needs this many finite pairs for an r
%   'perAnimalTest' : 'parametric' (default) | 'permutation'
%   'nPermAnimal'   : 10000 (permutation mode only)
%   'alpha'         : 0.05 -- threshold for marking significant animals
%   'markOn'        : 'q' (default, FDR-corrected) | 'p' (uncorrected)
%   'groupPermTest' : true -- exact fast-vs-slow relabeling test per panel
%   'maxExactPart'  : 20000 -- above this many partitions, Monte Carlo
%   'nPermGroup'    : 10000 -- Monte Carlo draws when not exact
%   'cueWin'        : [0 2] -- cue-period window, s from tone onset
%   'respWin'       : [2 4] -- response-period window, e.g. [2 5]
%                     These build the two columns of the figure. Passing
%                     'windows' explicitly (an N x 2 matrix) overrides both
%                     and allows more than two windows.
%   'trialTypes', 'windows', 'windowNames', 'fastLearners', 'slowLearners'
%                   : as in buildDAScalarTable
%   'fastColor'     : [0.85 0.33 0.10]
%   'slowColor'     : [0.00 0.45 0.74]
%   'figSaveDir'    : ''
%
% OUTPUT
%   .rTable   one row per animal x stream x window: r, p, q (BH within
%             panel, across animals), n sessions, group
%   .stats    one row per panel: mean r, sign-flip p, n animals; the same
%             separately for fast and slow; and the relabeling test:
%             diffZ_fastMinusSlow, p_group (two-sided), p_fastLower,
%             p_fastHigher, nPartitions, pFloor_group, groupTestExact
%   .fig1, .fig2

%% ---- parse ----
p = inputParser;
p.addParameter('animals', {}, @(c) isempty(c) || ischar(c) || iscellstr(c) || isstring(c));
p.addParameter('signal', 'global', @(x) (ischar(x) || isstring(x)) || (isnumeric(x) && isscalar(x)));
p.addParameter('xVar', 'dprime', @(s) any(strcmpi(string(s), ["dprime","sessOrder"])));
p.addParameter('corrType', 'Pearson', @(s) any(strcmpi(string(s), ["Pearson","Spearman"])));
p.addParameter('minSessions', 4, @(x) isnumeric(x) && isscalar(x) && x >= 3);
p.addParameter('perAnimalTest', 'parametric', @(s) any(strcmpi(string(s), ["parametric","permutation"])));
p.addParameter('nPermAnimal', 10000, @(x) isnumeric(x) && isscalar(x) && x >= 100);
p.addParameter('alpha', 0.05, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
p.addParameter('markOn', 'q', @(s) any(strcmpi(string(s), ["q","p"])));
p.addParameter('groupPermTest', true, @(x) islogical(x) && isscalar(x));
p.addParameter('maxExactPart', 20000, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('nPermGroup', 10000, @(x) isnumeric(x) && isscalar(x) && x >= 100);
p.addParameter('trialTypes', {'hit','cr'}, @(c) iscellstr(c) || isstring(c));
p.addParameter('cueWin',  [0 2], @(x) isnumeric(x) && numel(x) == 2 && x(2) > x(1));
p.addParameter('respWin', [2 4], @(x) isnumeric(x) && numel(x) == 2 && x(2) > x(1));
p.addParameter('windows', [0 2; 2 4], @(x) isnumeric(x) && size(x,2) == 2);
p.addParameter('windowNames', {'cue','response'}, @(c) iscellstr(c) || isstring(c));
p.addParameter('fastLearners', {'m1044','m1045','m1092','m1094'}, @iscell);
p.addParameter('slowLearners', {'m1048','m1049','m1613','m1859','m1873'}, @iscell);
p.addParameter('fastColor', [0.85 0.33 0.10], @(x) isnumeric(x) && numel(x)==3);
p.addParameter('slowColor', [0.00 0.45 0.74], @(x) isnumeric(x) && numel(x)==3);
p.addParameter('figSaveDir', '', @(s) ischar(s) || isstring(s));
p.parse(varargin{:});
opt = p.Results;

trialTypes  = cellstr(opt.trialTypes);
windowNames = cellstr(opt.windowNames);
% 'windows' given explicitly wins; otherwise build it from cueWin/respWin.
if ismember('windows', p.UsingDefaults)
    W = [opt.cueWin(:)'; opt.respWin(:)'];
else
    W = opt.windows;
end
assert(size(W,1) == numel(windowNames), ...
    '%d windows but %d windowNames -- pass matching ''windowNames''.', size(W,1), numel(windowNames));
xVar = char(string(opt.xVar));
if strcmpi(xVar, 'dprime'), xLabel = 'd'''; else, xLabel = 'session order (within animal)'; end

T = buildDAScalarTable(rezC, sessT, 'signal', opt.signal, 'trialTypes', trialTypes, ...
    'windows', W, 'windowNames', windowNames, ...
    'fastLearners', opt.fastLearners, 'slowLearners', opt.slowLearners);

if isnumeric(opt.signal), sigLabel = sprintf('motif %d DA', opt.signal); sigTag = sprintf('motif%02d', opt.signal);
else, sigLabel = 'global DA'; sigTag = 'global'; end

if ~isempty(opt.animals)
    want = cellstr(string(opt.animals));
    unknown = setdiff(want, unique(T.animal));
    assert(isempty(unknown), ['Not found (or in neither fast/slow group): %s. ' ...
        'Available: %s'], strjoin(unknown, ', '), strjoin(unique(T.animal)', ', '));
    T = T(ismember(T.animal, want), :);
end
animals = unique(T.animal, 'stable');
singleAnimal = numel(animals) == 1;
nR = numel(trialTypes); nC = size(W,1);

%% ---- per-animal r ----
rRows = struct('animal',{},'group',{},'trialType',{},'window',{},'r',{},'p',{},'nSessions',{});
for t = 1:nR
    for w = 1:nC
        for a = 1:numel(animals)
            m = strcmp(T.animal, animals{a}) & strcmp(T.trialType, trialTypes{t}) & ...
                strcmp(T.window, windowNames{w}) & isfinite(T.(xVar)) & isfinite(T.DA);
            n = sum(m);
            r = NaN; pv = NaN;
            if n >= opt.minSessions && std(T.(xVar)(m)) > 0 && std(T.DA(m)) > 0
                xa = T.(xVar)(m); ya = T.DA(m);
                [r, pv] = corr(xa, ya, 'Type', char(opt.corrType));
                if strcmpi(opt.perAnimalTest, 'permutation')
                    pv = permCorrP(xa, ya, r, char(opt.corrType), opt.nPermAnimal);
                end
            end
            g = T.group(find(strcmp(T.animal, animals{a}), 1));
            rRows(end+1) = struct('animal', animals{a}, 'group', g{1}, ...
                'trialType', trialTypes{t}, 'window', windowNames{w}, ...
                'r', r, 'p', pv, 'nSessions', n); %#ok<AGROW>
        end
    end
end
rTable = struct2table(rRows);

% BH-FDR WITHIN each panel, across its animals: that is the family a reader
% scans when asking "which animals show it". Correcting across all 36 tests
% would be defensible too, but would conflate unrelated panels.
rTable.q = nan(height(rTable), 1);
for t = 1:nR
    for w = 1:nC
        m = strcmp(rTable.trialType, trialTypes{t}) & strcmp(rTable.window, windowNames{w}) & isfinite(rTable.p);
        rTable.q(m) = bhFDR(rTable.p(m));
    end
end
rTable = movevars(rTable, 'q', 'After', 'p');

%% ---- figure 1: per-animal fit lines over the scatter ----
fig1 = figure('Color', 'w', 'Position', [80 80 390*nC 360*nR]);
ax1 = gobjects(nR, nC);
for t = 1:nR
    for w = 1:nC
        ax1(t,w) = subplot(nR, nC, (t-1)*nC + w); hold(ax1(t,w), 'on');
        for a = 1:numel(animals)
            m = strcmp(T.animal, animals{a}) & strcmp(T.trialType, trialTypes{t}) & ...
                strcmp(T.window, windowNames{w}) & isfinite(T.(xVar)) & isfinite(T.DA);
            if ~any(m), continue; end
            Ta = T(m, :);
            c = opt.slowColor; if strcmp(Ta.group{1}, 'fast'), c = opt.fastColor; end
            pre = ~Ta.isDay4;
            scatter(ax1(t,w), Ta.(xVar)(pre),  Ta.DA(pre),  28, c, 'LineWidth', 1);
            scatter(ax1(t,w), Ta.(xVar)(~pre), Ta.DA(~pre), 28, c, 'filled', 'MarkerFaceAlpha', 0.8);
            if height(Ta) >= opt.minSessions && std(Ta.(xVar)) > 0
                b = polyfit(Ta.(xVar), Ta.DA, 1);
                xf = [min(Ta.(xVar)) max(Ta.(xVar))];
                plot(ax1(t,w), xf, polyval(b, xf), '-', 'Color', [c 0.8], 'LineWidth', 1.6);
            end
        end
        xline(ax1(t,w), 0, ':', 'Color', [0.7 0.7 0.7]);
        yline(ax1(t,w), 0, ':', 'Color', [0.7 0.7 0.7]);
        ttl = sprintf('%s | %s (%g-%g s)', upperType(trialTypes{t}), windowNames{w}, W(w,1), W(w,2));
        if singleAnimal
            k = strcmp(rTable.trialType, trialTypes{t}) & strcmp(rTable.window, windowNames{w});
            if any(k) && isfinite(rTable.r(k))
                ttl = sprintf('%s\nr = %.2f, p = %.3g (n = %d sessions)', ttl, ...
                    rTable.r(k), rTable.p(k), rTable.nSessions(k));
            else
                ttl = sprintf('%s\n(< %d usable sessions)', ttl, opt.minSessions);
            end
        end
        title(ax1(t,w), ttl, 'FontWeight', 'normal');
        xlabel(ax1(t,w), xLabel, 'FontWeight', 'bold');
        ylabel(ax1(t,w), sprintf('mean %s', sigLabel), 'FontWeight', 'bold');
        set(ax1(t,w), 'TickDir', 'out', 'Box', 'off');
    end
end
if singleAnimal
    sgtitle(fig1, sprintf('%s (%s) | %s vs %s | open = pre-day 4, filled = day 4+', ...
        animals{1}, T.group{1}, sigLabel, xLabel), 'FontWeight', 'bold', 'FontSize', 11);
else
    sgtitle(fig1, sprintf('%s vs %s | one line per animal (open = pre-day 4, filled = day 4+)', sigLabel, xLabel), ...
        'FontWeight', 'bold', 'FontSize', 11);
end

%% ---- figure 2: per-animal r + group test ----
stats = table();
fig2 = [];
if singleAnimal
    fprintf('Single animal selected -- Figure 2 (group test) skipped; r and p are in Figure 1.\n');
else
fig2 = figure('Color', 'w', 'Position', [120 120 320*nC 340*nR]);
ax2 = gobjects(nR, nC);
for t = 1:nR
    for w = 1:nC
        ax2(t,w) = subplot(nR, nC, (t-1)*nC + w); hold(ax2(t,w), 'on');
        m = strcmp(rTable.trialType, trialTypes{t}) & strcmp(rTable.window, windowNames{w}) & isfinite(rTable.r);
        R = rTable(m, :);

        [mAll, pAll] = signFlipMean(R.r);
        isF = strcmp(R.group, 'fast');
        [mF, pF] = signFlipMean(R.r(isF));
        [mS, pS] = signFlipMean(R.r(~isF));

        G = struct('diff', NaN, 'p', NaN, 'pLow', NaN, 'pHigh', NaN, ...
                   'nPart', NaN, 'floor', NaN, 'exact', false);
        if opt.groupPermTest
            G = groupRelabelTest(R.r, isF, opt.maxExactPart, opt.nPermGroup);
        end

        stats = [stats; table(string(trialTypes{t}), string(windowNames{w}), ...
            mAll, pAll, height(R), mF, pF, sum(isF), mS, pS, sum(~isF), ...
            G.diff, G.p, G.pLow, G.pHigh, G.nPart, G.floor, G.exact, ...
            'VariableNames', {'trialType','window','meanR','p_signFlip','nAnimals', ...
                              'meanR_fast','p_fast','nFast','meanR_slow','p_slow','nSlow', ...
                              'diffZ_fastMinusSlow','p_group','p_fastLower','p_fastHigher', ...
                              'nPartitions','pFloor_group','groupTestExact'})]; %#ok<AGROW>

        jit = (rand(height(R),1) - 0.5) * 0.25;
        if strcmpi(opt.markOn, 'q'), sigV = R.q; else, sigV = R.p; end
        isSig = isfinite(sigV) & sigV < opt.alpha;
        scatter(ax2(t,w), 1 + jit(isF),  R.r(isF),  60, opt.fastColor, 'filled');
        scatter(ax2(t,w), 1 + jit(~isF), R.r(~isF), 60, opt.slowColor, 'filled');
        % significant animals: black ring plus ID label, so the figure shows
        % WHICH animals clear the threshold, not just how many
        if any(isSig)
            scatter(ax2(t,w), 1 + jit(isSig), R.r(isSig), 110, 'k', 'LineWidth', 1.3);
            text(ax2(t,w), 1 + jit(isSig) + 0.07, R.r(isSig), R.animal(isSig), ...
                'FontSize', 7, 'VerticalAlignment', 'middle');
        end
        % Summaries, left to right: fast, slow, all -- between the animal
        % scatter and the grand summary. All three are mean +/- SEM taken
        % in Fisher-z and back-transformed, so they are on the same footing
        % as the group tests (and the bars are slightly asymmetric in r).
        xF = 1.25; xS = 1.38; xA = 1.55;
        drawZSummary(ax2(t,w), xF, R.r(isF),  opt.fastColor);
        drawZSummary(ax2(t,w), xS, R.r(~isF), opt.slowColor);
        drawZSummary(ax2(t,w), xA, R.r,       [0 0 0]);
        yline(ax2(t,w), 0, '--', 'Color', [0.5 0.5 0.5]);
        xlim(ax2(t,w), [0.72 1.68]); ylim(ax2(t,w), [-1 1]);
        set(ax2(t,w), 'XTick', [1 xF xS xA], ...
            'XTickLabel', {'animals', sprintf('fast\n(n=%d)', sum(isF)), ...
                           sprintf('slow\n(n=%d)', sum(~isF)), 'all'}, ...
            'TickDir', 'out', 'Box', 'off');
        ylabel(ax2(t,w), sprintf('per-animal r (%s)', lower(char(opt.corrType))), 'FontWeight', 'bold');
        ttl = sprintf('%s | %s\nmean r = %.2f, p = %.3g (n = %d mice) | %d sig.', ...
            upperType(trialTypes{t}), windowNames{w}, mAll, pAll, height(R), sum(isSig));
        if opt.groupPermTest && isfinite(G.p)
            % the floor is printed beside p so it can't be misread as small
            ttl = sprintf('%s\nfast vs slow: p = %.3g (floor %.3g, %d relabelings)', ...
                ttl, G.p, G.floor, G.nPart);
        end
        title(ax2(t,w), ttl, 'FontWeight', 'normal');
    end
end
sgtitle(fig2, sprintf(['Within-animal r: %s vs %s | group p: exact sign-flip | ' ...
    'ringed = %s < %.2f (%s)'], sigLabel, xLabel, lower(char(opt.markOn)), opt.alpha, ...
    lower(char(opt.perAnimalTest))), 'FontWeight', 'bold', 'FontSize', 10);
end

nDrop = sum(~isfinite(rTable.r));
if nDrop > 0
    fprintf('%d animal x panel combination(s) had < %d usable sessions (or no variance) and no r.\n', ...
        nDrop, opt.minSessions);
end
fprintf('\n---- Per-animal correlations (%s, %s p; q = BH-FDR within panel) ----\n', ...
    char(opt.corrType), lower(char(opt.perAnimalTest)));
disp(sortrows(rTable, {'trialType','window','group','animal'}));
if ~singleAnimal
    fprintf('---- Group level (exact sign-flip on Fisher-z r) ----\n');
    disp(stats);
end

out = struct('rTable', rTable, 'stats', stats, 'T', T, 'fig1', fig1, 'fig2', fig2);

if strlength(strtrim(string(opt.figSaveDir))) > 0
    d = char(string(opt.figSaveDir));
    if exist(d, 'dir') ~= 7, mkdir(d); end
    if singleAnimal, aTag = ['_' animals{1}]; else, aTag = ''; end
    for pair = {{fig1, 'perAnimalFits'}, {fig2, 'perAnimalR'}}
        if isempty(pair{1}{1}), continue; end
        wTag = strjoin(arrayfun(@(i) sprintf('%g-%g', W(i,1), W(i,2)), 1:size(W,1), ...
            'UniformOutput', false), '_');
        f = fullfile(d, sprintf('DAScalarVs%s_%s_%s%s_win%s.pdf', xVar, pair{1}{2}, sigTag, aTag, wTag));
        set(pair{1}{1}, 'InvertHardcopy', 'off');
        print(pair{1}{1}, f, '-dpdf', '-painters', '-bestfit');
        fprintf('Saved:\n  %s\n', f);
    end
end
end

%% ========================================================================
function [meanR, p] = signFlipMean(r)
% Mean of Fisher-z per-animal r's, tested by EXACT sign-flip permutation
% over all 2^n sign assignments. Returns meanR back in r units.
r = r(isfinite(r));
n = numel(r);
if n < 2, meanR = NaN; p = NaN; return; end
z = atanh(clip(r));
obs = mean(z);
signs = 1 - 2 * (dec2bin(0:2^n - 1, n) - '0');     % 2^n x n, entries +/-1
null = signs * z(:) / n;
p = sum(abs(null) >= abs(obs) - 1e-12) / size(signs, 1);
meanR = tanh(obs);
end

function drawZSummary(ax, x, r, col)
% Mean +/- SEM in Fisher-z, back-transformed to r. With one animal there is
% no SEM, so only the point is drawn -- an error bar from n = 1 would be
% meaningless. With n = 2 the SEM is computable but crude; it is drawn,
% and the n in the tick label is the reader's cue to discount it.
r = r(isfinite(r));
if isempty(r), return; end
z = atanh(clip(r(:)));
mz = mean(z);
if numel(z) > 1
    se = std(z) / sqrt(numel(z));
    errorbar(ax, x, tanh(mz), tanh(mz) - tanh(mz - se), tanh(mz + se) - tanh(mz), 'o', ...
        'Color', col, 'MarkerFaceColor', col, 'MarkerEdgeColor', col, ...
        'MarkerSize', 7, 'LineWidth', 1.4, 'CapSize', 8);
else
    plot(ax, x, tanh(mz), 'o', 'Color', col, 'MarkerFaceColor', col, 'MarkerSize', 7);
end
end

function G = groupRelabelTest(r, isF, maxExact, nPerm)
% Exact test of fast vs slow per-animal r by relabeling. Statistic: mean
% Fisher-z of the "fast" set minus mean of the rest. The null enumerates
% every choice of nFast animals out of n (exact) when that count is at most
% maxExact, otherwise samples nPerm random choices.
G = struct('diff', NaN, 'p', NaN, 'pLow', NaN, 'pHigh', NaN, ...
           'nPart', NaN, 'floor', NaN, 'exact', false);
ok = isfinite(r);
r = r(ok); isF = logical(isF(ok));
n = numel(r); nF = sum(isF);
if nF < 1 || nF >= n, return; end             % need both groups present

z = atanh(clip(r(:)));
st = @(m) mean(z(m)) - mean(z(~m));
obs = st(isF);

nPart = nchoosek(n, nF);
if nPart <= maxExact
    C = nchoosek(1:n, nF);
    null = zeros(size(C,1), 1);
    for i = 1:size(C,1)
        m = false(n,1); m(C(i,:)) = true;
        null(i) = st(m);
    end
    tol = 1e-12;
    G.p     = mean(abs(null) >= abs(obs) - tol);
    G.pLow  = mean(null <= obs + tol);         % fast LOWER than slow
    G.pHigh = mean(null >= obs - tol);         % fast HIGHER than slow
    G.floor = 1 / nPart;
    G.exact = true;
else
    null = zeros(nPerm, 1);
    for i = 1:nPerm
        m = false(n,1); m(randperm(n, nF)) = true;
        null(i) = st(m);
    end
    G.p     = (1 + sum(abs(null) >= abs(obs))) / (nPerm + 1);
    G.pLow  = (1 + sum(null <= obs)) / (nPerm + 1);
    G.pHigh = (1 + sum(null >= obs)) / (nPerm + 1);
    G.floor = 1 / (nPerm + 1);
end
G.diff  = obs;
G.nPart = nPart;
end

function p = permCorrP(x, y, rObs, corrType, nPerm)
% Two-sided permutation p for one animal: shuffle y across that animal's
% own sessions. Exhaustive when n! <= nPerm (exact), else Monte Carlo.
n = numel(x);
if factorial(n) <= nPerm
    P = perms(1:n);
    rN = arrayfun(@(i) corr(x, y(P(i,:)), 'Type', corrType), 1:size(P,1));
    p = mean(abs(rN) >= abs(rObs) - 1e-12);
else
    rN = zeros(nPerm, 1);
    for s = 1:nPerm
        rN(s) = corr(x, y(randperm(n)), 'Type', corrType);
    end
    p = (1 + sum(abs(rN) >= abs(rObs) - 1e-12)) / (nPerm + 1);
end
end

function q = bhFDR(p)
p = p(:); n = numel(p);
if n == 0, q = p; return; end
[ps, si] = sort(p);
qs = ps .* n ./ (1:n)';
qs = flipud(cummin(flipud(qs)));
qs = min(qs, 1);
q = nan(n, 1); q(si) = qs;
end

function r = clip(r)
r = max(min(r, 0.999), -0.999);
end

function s = upperType(tt)
switch lower(tt)
    case 'cr', s = 'CR';
    case 'fa', s = 'FA';
    otherwise, s = [upper(tt(1)) tt(2:end)];
end
end