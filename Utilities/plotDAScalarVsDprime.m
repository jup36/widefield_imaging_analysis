function out = plotDAScalarVsDprime(rezC, sessT, varargin)
%PLOTDASCALARVSDPRIME
%   Session-level scatter of mean DA activity against behavioural d', in a
%   2 x 2 grid: rows = trial stream (Hit, CR), columns = window (cue 0-2 s,
%   response 2-4 s). One point per session.
%
%   ENCODING
%     colour : fast (orange) vs slow (blue) learner
%     fill   : pre-day-4 session = open circle, day-4-onward = filled
%
%   THE SCALAR is the session's mean PETH averaged over the window --
%   equivalently, the mean over trials of each trial's window mean, since
%   both averages are linear. It uses the same normalized traces as every
%   other DA figure (session-level z-score by default), so values are
%   comparable across sessions.
%
%   STATISTICS -- why not just corr()
%   ---------------------------------
%   Sessions are nested in animals, so a pooled Pearson r across ~90
%   sessions treats non-independent points as independent. Worse, it mixes
%   two sources of covariation: WITHIN an animal (does DA rise as that
%   animal's d' rises?) and BETWEEN animals (do animals with higher d'
%   simply have different DA baselines?). Only the first speaks to learning.
%
%   So each panel reports the ANIMAL-CENTERED r -- the correlation after
%   subtracting each animal's own mean d' and mean DA -- tested against a
%   WITHIN-ANIMAL permutation null (d' shuffled only among each animal's own
%   sessions). Centering is what makes the reported number mean what the
%   test tests: in simulation, animals differing in both d' and DA but with
%   NO within-animal link gave a pooled r of 0.98 and a centered r near 0.
%   Printing the pooled value would badly overstate the effect. The pooled r
%   is still returned in .stats, for reference only. The permutation was
%   checked for calibration (false-positive rate <= 0.05 with no
%   within-animal effect).
%
%   CAVEAT: d' rises with training, so a DA-d' correlation cannot by itself
%   separate "DA tracks performance" from "DA changes with time in task".
%   Comparing against the same plot with session order on x would help.
%
%   out = plotDAScalarVsDprime(rezC, sessT, ...)
%
% NAME-VALUE
%   'signal'        : 'global' (default) or motif number 1..K
%   'trialTypes'    : {'hit','cr'}   -- rows
%   'windows'       : [0 2; 2 4]     -- columns
%   'windowNames'   : {'cue','response'}
%   'fastLearners'  : {'m1044','m1045','m1092','m1094'}
%   'slowLearners'  : {'m1048','m1049','m1613','m1859','m1873'}
%   'fastColor'     : [0.85 0.33 0.10]
%   'slowColor'     : [0.00 0.45 0.74]
%   'markerSize'    : 55
%   'nPerm'         : 5000   within-animal permutation draws
%   'showFit'       : true   pooled least-squares line (descriptive only)
%   'rngSeed'       : 1
%   'figSaveDir'    : ''
%
% OUTPUT
%   .T      table: one row per session x stream x window, with d', DA,
%           group, isDay4, animal -- everything plotted, for further stats
%   .stats  table: r, p_withinAnimalPerm, n sessions, n animals per panel
%   .fig, .ax

%% ---- parse ----
p = inputParser;
p.addParameter('signal', 'global', @(x) (ischar(x) || isstring(x)) || (isnumeric(x) && isscalar(x)));
p.addParameter('trialTypes', {'hit','cr'}, @(c) iscellstr(c) || isstring(c));
p.addParameter('windows', [0 2; 2 4], @(x) isnumeric(x) && size(x,2) == 2);
p.addParameter('windowNames', {'cue','response'}, @(c) iscellstr(c) || isstring(c));
p.addParameter('fastLearners', {'m1044','m1045','m1092','m1094'}, @iscell);
p.addParameter('slowLearners', {'m1048','m1049','m1613','m1859','m1873'}, @iscell);
p.addParameter('fastColor', [0.85 0.33 0.10], @(x) isnumeric(x) && numel(x)==3);
p.addParameter('slowColor', [0.00 0.45 0.74], @(x) isnumeric(x) && numel(x)==3);
p.addParameter('markerSize', 55, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('nPerm', 5000, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('showFit', true, @(x) islogical(x) && isscalar(x));
p.addParameter('rngSeed', 1, @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('figSaveDir', '', @(s) ischar(s) || isstring(s));
p.parse(varargin{:});
opt = p.Results;
if ~isempty(opt.rngSeed), rng(opt.rngSeed); end

trialTypes  = cellstr(opt.trialTypes);
windowNames = cellstr(opt.windowNames);
W = opt.windows;
assert(size(W,1) == numel(windowNames), 'windows and windowNames must match in number.');
assert(ismember('isDay4', sessT.Properties.VariableNames), ...
    'sessT has no isDay4 column -- rerun batch_DAPeth_globalAndMotifWeighted with day4MarkC.');

%% ---- signal selector ----
if isnumeric(opt.signal)
    k = opt.signal;
    getTrace = @(pe) pe.motifMean(k, :);
    sigLabel = sprintf('motif %d DA', k);  sigTag = sprintf('motif%02d', k);
else
    assert(strcmpi(opt.signal, 'global'), 'signal must be ''global'' or a motif number.');
    getTrace = @(pe) pe.globalMean;
    sigLabel = 'global DA';  sigTag = 'global';
end

%% ---- build the long table: session x stream x window ----
tint = rezC{sessT.a(1), sessT.s(1)}.meta.tint;
winI = arrayfun(@(i) tint >= W(i,1) & tint < W(i,2), 1:size(W,1), 'UniformOutput', false);
for i = 1:numel(winI)
    assert(any(winI{i}), 'Window %s selects no bins of the time grid.', windowNames{i});
end

normTag = '';
r0 = rezC{sessT.a(1), sessT.s(1)};
if isfield(r0.meta, 'normMode') && startsWith(r0.meta.normMode, 'zscore'), normTag = ' (z)'; end

animal = {}; group = {}; tt_ = {}; win_ = {}; dp = []; da = []; d4 = []; hdr = {};
for i = 1:height(sessT)
    r = rezC{sessT.a(i), sessT.s(i)};
    aId = char(sessT.animal{i});
    if ismember(aId, opt.fastLearners), g = 'fast';
    elseif ismember(aId, opt.slowLearners), g = 'slow';
    else, continue;                       % not in either group (e.g. m1237)
    end
    for t = 1:numel(trialTypes)
        tt = trialTypes{t};
        if ~isfield(r.peth, tt) || r.peth.(tt).n == 0, continue; end
        y = getTrace(r.peth.(tt));
        for w = 1:size(W,1)
            animal{end+1,1} = aId;               %#ok<AGROW>
            group{end+1,1}  = g;                 %#ok<AGROW>
            tt_{end+1,1}    = tt;                %#ok<AGROW>
            win_{end+1,1}   = windowNames{w};    %#ok<AGROW>
            hdr{end+1,1}    = char(sessT.header{i}); %#ok<AGROW>
            dp(end+1,1)     = sessT.dprime(i);   %#ok<AGROW>
            da(end+1,1)     = mean(y(winI{w}), 'omitnan'); %#ok<AGROW>
            d4(end+1,1)     = sessT.isDay4(i);   %#ok<AGROW>
        end
    end
end
T = table(animal, group, hdr, tt_, win_, dp, da, logical(d4), ...
    'VariableNames', {'animal','group','header','trialType','window','dprime','DA','isDay4'});

nNaN = sum(~isfinite(T.dprime)) / (numel(trialTypes) * size(W,1));
fprintf('%d session(s) have no d'' and are excluded from every panel.\n', round(nNaN));

%% ---- plot ----
nR = numel(trialTypes); nC = size(W,1);
fig = figure('Color', 'w', 'Position', [100 100 390*nC 360*nR]);
ax = gobjects(nR, nC);
stats = table();

for t = 1:nR
    for w = 1:nC
        ax(t,w) = subplot(nR, nC, (t-1)*nC + w);
        hold(ax(t,w), 'on');
        m = strcmp(T.trialType, trialTypes{t}) & strcmp(T.window, windowNames{w}) & ...
            isfinite(T.dprime) & isfinite(T.DA);
        Tp = T(m, :);

        % draw pre-day-4 (open) first so filled day-4 points sit on top
        for d4Val = [false true]
            for gName = {'slow','fast'}
                sel = strcmp(Tp.group, gName{1}) & Tp.isDay4 == d4Val;
                if ~any(sel), continue; end
                c = opt.slowColor; if strcmp(gName{1}, 'fast'), c = opt.fastColor; end
                if d4Val
                    scatter(ax(t,w), Tp.dprime(sel), Tp.DA(sel), opt.markerSize, c, 'filled', ...
                        'MarkerEdgeColor', c, 'MarkerFaceAlpha', 0.85, 'HandleVisibility', 'off');
                else
                    scatter(ax(t,w), Tp.dprime(sel), Tp.DA(sel), opt.markerSize, c, ...
                        'LineWidth', 1.2, 'HandleVisibility', 'off');
                end
            end
        end

        % stats: pooled r, within-animal permutation null
        [rObs, pPerm, rPooled] = withinAnimalPerm(Tp.dprime, Tp.DA, Tp.animal, opt.nPerm);
        nA = numel(unique(Tp.animal));
        stats = [stats; table(string(trialTypes{t}), string(windowNames{w}), ...
            rObs, pPerm, rPooled, height(Tp), nA, ...
            'VariableNames', {'trialType','window','r_withinAnimal','p_withinAnimalPerm', ...
                              'r_pooled_reference','nSessions','nAnimals'})]; %#ok<AGROW>

        if opt.showFit && height(Tp) >= 3
            b = polyfit(Tp.dprime, Tp.DA, 1);
            xf = linspace(min(Tp.dprime), max(Tp.dprime), 50);
            plot(ax(t,w), xf, polyval(b, xf), '-', 'Color', [0.35 0.35 0.35], ...
                'LineWidth', 1.2, 'HandleVisibility', 'off');
        end

        xline(ax(t,w), 0, ':', 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
        yline(ax(t,w), 0, ':', 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');

        title(ax(t,w), sprintf('%s | %s (%g-%g s)\nr_{within} = %.2f, p_{perm} = %.3g (n = %d sess, %d mice)', ...
            upperType(trialTypes{t}), windowNames{w}, W(w,1), W(w,2), rObs, pPerm, height(Tp), nA), ...
            'FontWeight', 'normal', 'FontSize', 10);
        xlabel(ax(t,w), 'd''', 'FontWeight', 'bold');
        ylabel(ax(t,w), sprintf('mean %s%s', sigLabel, normTag), 'FontWeight', 'bold');
        set(ax(t,w), 'TickDir', 'out', 'Box', 'off', 'FontSize', 10);
    end
end

% shared y within each row: windows of one stream are compared on one scale
for t = 1:nR
    yl = cell2mat(arrayfun(@(h) ylim(h), ax(t,:)', 'UniformOutput', false));
    set(ax(t,:), 'YLim', [min(yl(:,1)), max(yl(:,2))]);
end
xl = cell2mat(arrayfun(@(h) xlim(h), ax(:), 'UniformOutput', false));
set(ax(:), 'XLim', [min(xl(:,1)), max(xl(:,2))]);

% legend: explicit proxy handles, so it documents the encoding rather than
% whichever scatter happened to draw first
hold(ax(1,1), 'on');
h1 = scatter(ax(1,1), nan, nan, opt.markerSize, opt.fastColor, 'filled');
h2 = scatter(ax(1,1), nan, nan, opt.markerSize, opt.slowColor, 'filled');
h3 = scatter(ax(1,1), nan, nan, opt.markerSize, [0.4 0.4 0.4], 'LineWidth', 1.2);
h4 = scatter(ax(1,1), nan, nan, opt.markerSize, [0.4 0.4 0.4], 'filled');
legend(ax(1,1), [h1 h2 h3 h4], {'fast learner','slow learner','pre-day 4','day 4 onward'}, ...
    'Location', 'best', 'Box', 'off', 'FontSize', 9);

sgtitle(fig, sprintf('%s vs d'' | one point per session | p from within-animal permutation', sigLabel), ...
    'FontWeight', 'bold', 'FontSize', 11);

disp(stats);
out = struct('T', T, 'stats', stats, 'fig', fig, 'ax', ax);

if strlength(strtrim(string(opt.figSaveDir))) > 0
    d = char(string(opt.figSaveDir));
    if exist(d, 'dir') ~= 7, mkdir(d); end
    f = fullfile(d, sprintf('DAScalarVsDprime_%s.pdf', sigTag));
    set(fig, 'InvertHardcopy', 'off');
    print(fig, f, '-dpdf', '-painters', '-bestfit');
    fprintf('Saved:\n  %s\n', f);
end
end

%% ========================================================================
function [rObs, p, rPooled] = withinAnimalPerm(x, y, animal, nPerm)
% Animal-centered Pearson r, tested by shuffling x only WITHIN each animal.
% Centering removes every animal's mean from both variables, so rObs
% measures within-animal covariation alone -- the quantity the permutation
% null is built to test. The null is centered the same way on every draw.
x = x(:); y = y(:);
rPooled = NaN; rObs = NaN; p = NaN;
if numel(x) < 3, return; end
rPooled = corr(x, y);
[~, ~, g] = unique(animal);
grpIdx = accumarray(g, (1:numel(x))', [], @(v) {v});

xc = centerWithin(x, grpIdx);
yc = centerWithin(y, grpIdx);
if all(xc == 0) || all(yc == 0), return; end   % no within-animal variance
rObs = corr(xc, yc);

rNull = nan(nPerm, 1);
for s = 1:nPerm
    xs = x;
    for k = 1:numel(grpIdx)
        ix = grpIdx{k};
        if numel(ix) > 1, xs(ix) = x(ix(randperm(numel(ix)))); end
    end
    rNull(s) = corr(centerWithin(xs, grpIdx), yc);
end
p = (1 + sum(abs(rNull) >= abs(rObs))) / (nPerm + 1);
end

function v = centerWithin(v, grpIdx)
for k = 1:numel(grpIdx)
    ix = grpIdx{k};
    v(ix) = v(ix) - mean(v(ix));
end
end

function s = upperType(tt)
switch lower(tt)
    case 'cr', s = 'CR';
    case 'fa', s = 'FA';
    otherwise, s = [upper(tt(1)) tt(2:end)];
end
end