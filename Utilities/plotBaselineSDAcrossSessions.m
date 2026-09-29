function out = plotBaselineSDAcrossSessions(rezC, sessT, mouseId, varargin)
%PLOTBASELINESDACROSSSESSIONS
%   Per-session pre-tone baseline SD in dF/F -- the denominator of the
%   z-scored PETHs -- plotted in chronological order for one animal.
%
%   WHY THIS MATTERS
%   ----------------
%   Under 'zscoreSession' the PETH is (dF/F - baseline) / baselineSD, with
%   one SD per signal pooled over that session's baseline bins. If baselineSD
%   FALLS across training -- a calmer animal, fewer spontaneous transients,
%   more stable imaging -- then z rises even with no change in dopamine
%   release. That would reproduce the apparent growth of the response-epoch
%   DA peak across sessions without any of it being real.
%
%   So the test is simple: if baselineSD declines with session order while
%   the z-scored response grows, the effect is at least partly the
%   denominator. If baselineSD is flat, the denominator is not the
%   explanation. The right panel plots peak response z against 1/baselineSD
%   directly, since that is the form the artifact would take -- an exact
%   1/SD relationship is what a pure denominator effect looks like.
%
%   The bottom panel converts back: peak z * baselineSD is the peak in
%   dF/F units, which is the quantity that should be reported if the
%   denominator turns out to be moving.
%
%   out = plotBaselineSDAcrossSessions(rezC, sessT, mouseId, ...)
%
% NAME-VALUE
%   'signal'     : 'global' (default) or motif number 1..K
%   'trialType'  : 'hit' (default) -- which stream's peak to relate SD to
%   'respWin'    : [2 4] s -- window for the peak used in panels 2 and 3
%   'markDay4'   : true -- open circles before day 4, filled from day 4 on
%                  (needs sessT.isDay4; ignored if absent)
%   'color'      : [0.15 0.15 0.15]
%   'figSaveDir' : ''
%
% OUTPUT
%   .T    per-session table: header, date, sessOrder, isDay4, baselineSD,
%         peakZ, peakDff (= peakZ * baselineSD), invSD
%   .fig, .ax

p = inputParser;
p.addParameter('signal', 'global', @(x) (ischar(x) || isstring(x)) || (isnumeric(x) && isscalar(x)));
p.addParameter('trialType', 'hit', @(s) ischar(s) || isstring(s));
p.addParameter('respWin', [2 4], @(x) isnumeric(x) && numel(x) == 2 && x(2) > x(1));
p.addParameter('markDay4', true, @(x) islogical(x) && isscalar(x));
p.addParameter('color', [0.15 0.15 0.15], @(x) isnumeric(x) && numel(x) == 3);
p.addParameter('figSaveDir', '', @(s) ischar(s) || isstring(s));
p.parse(varargin{:});
opt = p.Results;

mouseId = char(string(mouseId));
tt = lower(char(string(opt.trialType)));

rows = find(string(sessT.animal) == string(mouseId));
assert(~isempty(rows), 'mouseId %s not found in sessT.', mouseId);
S = sortrows(sessT(rows, :), 'date');

%% ---- signal selector ----
if isnumeric(opt.signal)
    k = opt.signal;
    getSD    = @(r) r.meta.baselineSD.motif(k);
    getTrace = @(pe) pe.motifMean(k, :);
    sigLabel = sprintf('motif %d DA', k);  sigTag = sprintf('motif%02d', k);
else
    assert(strcmpi(opt.signal, 'global'), 'signal must be ''global'' or a motif number.');
    getSD    = @(r) r.meta.baselineSD.global;
    getTrace = @(pe) pe.globalMean;
    sigLabel = 'global DA';  sigTag = 'global';
end

r0 = rezC{S.a(1), S.s(1)};
assert(isfield(r0.meta, 'baselineSD'), ...
    ['rezC has no meta.baselineSD -- rerun batch_DAPeth_globalAndMotifWeighted ' ...
     '(the field was added when this question came up).']);
tint  = r0.meta.tint;
respI = tint >= opt.respWin(1) & tint <= opt.respWin(2);

%% ---- gather ----
n = height(S);
sd = nan(n,1); pk = nan(n,1);
for i = 1:n
    r = rezC{S.a(i), S.s(i)};
    sd(i) = getSD(r);
    if isfield(r.peth, tt) && r.peth.(tt).n > 0
        y = getTrace(r.peth.(tt));
        pk(i) = max(y(respI));
    end
end
T = table(S.header, S.date, S.sessOrder, S.isDay4, sd, pk, pk .* sd, 1 ./ sd, ...
    'VariableNames', {'header','date','sessOrder','isDay4','baselineSD','peakZ','peakDff','invSD'});

if isfield(r0.meta, 'normMode') && ~startsWith(r0.meta.normMode, 'zscore')
    warning('BaselineSD:NotZscored', ...
        ['This run used normMode = %s, so the PETHs are not divided by baselineSD. ' ...
         'The SD is still shown, but peakZ is not a z-score.'], r0.meta.normMode);
end

%% ---- stats: does SD trend with session order? ----
ok = isfinite(T.sessOrder) & isfinite(T.baselineSD);
[rSD, pSD] = corr(T.sessOrder(ok), T.baselineSD(ok));
okP = ok & isfinite(T.peakZ);
[rZ,  pZ]  = corr(T.sessOrder(okP), T.peakZ(okP));
[rD,  pD]  = corr(T.sessOrder(okP), T.peakDff(okP));
[rI,  pI]  = corr(T.invSD(okP),     T.peakZ(okP));

fprintf('\n%s | %s | %s trials, peak in [%g %g] s  (n = %d sessions)\n', ...
    mouseId, sigLabel, upper(tt), opt.respWin(1), opt.respWin(2), sum(okP));
fprintf('  baselineSD vs session order : r = %+.2f, p = %.3g\n', rSD, pSD);
fprintf('  peak z     vs session order : r = %+.2f, p = %.3g\n', rZ,  pZ);
fprintf('  peak dF/F  vs session order : r = %+.2f, p = %.3g   <- the denominator-free version\n', rD, pD);
fprintf('  peak z     vs 1/baselineSD  : r = %+.2f, p = %.3g   <- a pure artifact would be near 1\n', rI, pI);

%% ---- plot ----
fig = figure('Color','w', 'Position', [100 100 1150 360]);
useD4 = opt.markDay4 && ismember('isDay4', S.Properties.VariableNames);
ax = gobjects(1,3);

ax(1) = subplot(1,3,1); hold on;
plotSplit(ax(1), T.sessOrder, T.baselineSD, T.isDay4, useD4, opt.color);
xlabel('session order (within animal)'); ylabel(sprintf('baseline SD of %s (dF/F)', sigLabel));
title(sprintf('baseline SD vs session\nr = %+.2f, p = %.3g', rSD, pSD), 'FontWeight','normal');

ax(2) = subplot(1,3,2); hold on;
plotSplit(ax(2), T.invSD, T.peakZ, T.isDay4, useD4, opt.color);
xlabel('1 / baseline SD'); ylabel(sprintf('peak %s (z)', sigLabel));
title(sprintf('peak z vs 1/SD\nr = %+.2f, p = %.3g', rI, pI), 'FontWeight','normal');

ax(3) = subplot(1,3,3); hold on;
plotSplit(ax(3), T.sessOrder, T.peakZ,   T.isDay4, useD4, [0.05 0.05 0.86]);
plotSplit(ax(3), T.sessOrder, T.peakDff, T.isDay4, useD4, [0.85 0.33 0.10]);
xlabel('session order (within animal)'); ylabel('peak response');
title(sprintf('z (blue) vs dF/F (orange)\nz: p = %.3g   dF/F: p = %.3g', pZ, pD), 'FontWeight','normal');

for i = 1:3
    set(ax(i), 'TickDir','out', 'Box','off');
    grid(ax(i), 'on'); hold(ax(i), 'off');
end

d4note = '';
if useD4, d4note = ' | open = pre-day 4, filled = day 4 onward'; end
sgtitle(fig, sprintf('%s | %s | %s trials%s', mouseId, sigLabel, upper(tt), d4note), ...
    'FontWeight','bold', 'FontSize', 11);

out = struct('T', T, 'fig', fig, 'ax', ax, ...
    'stats', struct('r_SD_vs_session', rSD, 'p_SD_vs_session', pSD, ...
                    'r_z_vs_session', rZ,  'p_z_vs_session', pZ, ...
                    'r_dff_vs_session', rD, 'p_dff_vs_session', pD, ...
                    'r_z_vs_invSD', rI,    'p_z_vs_invSD', pI));

if strlength(strtrim(string(opt.figSaveDir))) > 0
    d = char(string(opt.figSaveDir));
    if exist(d,'dir') ~= 7, mkdir(d); end
    f = fullfile(d, sprintf('baselineSD_%s_%s_%s.pdf', mouseId, sigTag, tt));
    set(fig, 'InvertHardcopy', 'off');
    print(fig, f, '-dpdf', '-painters', '-bestfit');
    fprintf('Saved:\n  %s\n', f);
end
end

%% ========================================================================
function plotSplit(ax, x, y, isD4, useD4, col)
% Filled = day-4 onward, open = before. Connecting line drawn through all
% points in x order so the trend is visible regardless of the split.
ok = isfinite(x) & isfinite(y);
[xs, si] = sort(x(ok)); ys = y(ok); ys = ys(si);
plot(ax, xs, ys, '-', 'Color', [col 0.35], 'LineWidth', 1, 'HandleVisibility','off');
if useD4
    d4 = logical(isD4(ok)); d4 = d4(si);
    scatter(ax, xs(~d4), ys(~d4), 45, col, 'LineWidth', 1.2);
    scatter(ax, xs(d4),  ys(d4),  45, col, 'filled', 'MarkerFaceAlpha', 0.85);
else
    scatter(ax, xs, ys, 45, col, 'filled', 'MarkerFaceAlpha', 0.85);
end
end