function out = plotDAPethLearningBins(rezC, sessT, mouseId, varargin)
%PLOTDAPETHLEARNINGBINS
%   One animal's trial-averaged DA PETH, split into three learning bins,
%   one panel per trial type (CR left, Hit right by default).
%
%   BINS (day4MarkC-anchored)
%     early        : every session BEFORE the animal's day-4 date
%     intermediate : first half of the sessions on or after day 4
%     late         : second half of the sessions on or after day 4
%   The day-4 sessions are ordered by recording date and split in half;
%   with an odd count the extra session goes to LATE (intermediate gets
%   floor(n/2)). Both choices are printed so the split is never implicit.
%
%   AVERAGING: each bin's trace is the mean of its sessions' mean PETHs --
%   sessions weighted equally, not trials, consistent with the two-stage
%   averaging used throughout this project. A session contributing 150
%   trials and one contributing 30 count the same.
%
%   out = plotDAPethLearningBins(rezC, sessT, mouseId, 'day4MarkC', day4MarkC, ...)
%
% INPUTS
%   rezC, sessT : from batch_DAPeth_globalAndMotifWeighted
%   mouseId     : e.g. 'm1045'
%
% NAME-VALUE
%   'day4MarkC'   : REQUIRED. {mouseId, "MMDDYY"; ...}. Only this animal's
%                   row is used, so extra rows (e.g. m1237) are harmless.
%   'signal'      : 'global' (default) or a motif number 1..K, which selects
%                   that motif's footprint-weighted DA (peth.motifMean(k,:))
%   'trialTypes'  : {'cr','hit'} -- one panel each, in this order
%   'showSEM'     : false. true shades +/- SEM ACROSS SESSIONS within a bin
%                   (not across trials -- sessions are the unit here).
%   'colors'      : [3 x 3] RGB for early/intermediate/late. Default matches
%                   the reference figure: pale lavender -> saturated blue.
%   'lineWidth'   : 2.5
%   'yLim'        : [] (shared auto) or [lo hi], applied to all panels
%   'xLim'        : [] -> full grid
%   'epochs'      : struct array with .name, .t (e.g. [0 2]), .filled.
%                   Default: tone [0 2] filled, response [2 4] outlined.
%   'showEpochs'  : true -- bar above each panel plus dashed boundaries
%   'figSaveDir'  : '' -> no save
%
% OUTPUT
%   .binMeans.<type>  [3 x T]   bin traces (rows early/intermediate/late)
%   .binSem.<type>    [3 x T]
%   .binSessions      {3 x 1}   headers in each bin
%   .fig, .ax

%% ---- parse ----
p = inputParser;
p.addParameter('day4MarkC', {}, @(c) iscell(c) || isstring(c));
p.addParameter('signal', 'global', @(x) (ischar(x) || isstring(x)) || (isnumeric(x) && isscalar(x)));
p.addParameter('trialTypes', {'cr','hit'}, @(c) iscellstr(c) || isstring(c));
p.addParameter('showSEM', false, @(x) islogical(x) && isscalar(x));
p.addParameter('colors', [0.80 0.80 0.91; 0.42 0.42 0.96; 0.05 0.05 0.86], ...
    @(x) isnumeric(x) && isequal(size(x), [3 3]));
p.addParameter('lineWidth', 2.5, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('yLim', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2 && x(2)>x(1)));
p.addParameter('xLim', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2 && x(2)>x(1)));
p.addParameter('epochs', struct('name', {'tone','response'}, 't', {[0 2], [2 4]}, ...
    'filled', {true, false}), @isstruct);
p.addParameter('showEpochs', true, @(x) islogical(x) && isscalar(x));
p.addParameter('figSaveDir', '', @(s) ischar(s) || isstring(s));
p.parse(varargin{:});
opt = p.Results;

mouseId = char(string(mouseId));
trialTypes = cellstr(opt.trialTypes);
binNames = {'early', 'intermediate', 'late'};

%% ---- this animal's sessions ----
rows = find(string(sessT.animal) == string(mouseId));
assert(~isempty(rows), 'mouseId %s not found in sessT.', mouseId);
S = sessT(rows, :);
S = sortrows(S, 'date');

%% ---- day-4 cutoff ----
assert(~isempty(opt.day4MarkC), ...
    'day4MarkC is required -- the early bin is defined as pre-day-4 sessions.');
d4 = opt.day4MarkC;
if isstring(d4), d4 = cellstr(d4); end
assert(size(d4, 2) == 2, 'day4MarkC must be Nx2: {mouseId, "MMDDYY"; ...}.');
hit = find(string(d4(:,1)) == string(mouseId), 1, 'first');
assert(~isempty(hit), 'mouseId %s has no row in day4MarkC.', mouseId);
cutoff = datetime(char(string(d4{hit, 2})), 'InputFormat', 'MMddyy');

% Compare on calendar day, so a same-day rerun header ('-1' suffix, parsed
% as a few seconds past midnight) is never pushed across the boundary.
isD4 = dateshift(S.date, 'start', 'day') >= cutoff;

iEarly = find(~isD4);
iD4    = find(isD4);
nInt   = floor(numel(iD4) / 2);
binIdx = {iEarly, iD4(1:nInt), iD4(nInt+1:end)};

fprintf('%s | day 4 = %s | early %d, intermediate %d, late %d sessions', ...
    mouseId, datestr(cutoff, 'mm/dd/yy'), numel(binIdx{1}), numel(binIdx{2}), numel(binIdx{3}));
if mod(numel(iD4), 2) == 1
    fprintf('  (odd day-4 count: extra session assigned to late)');
end
fprintf('\n');
for b = 1:3
    if isempty(binIdx{b})
        warning('DAPethBins:EmptyBin', '%s: the %s bin has no sessions and will not be drawn.', ...
            mouseId, binNames{b});
    end
end

%% ---- signal selector ----
if isnumeric(opt.signal)
    k = opt.signal;
    getTrace = @(pe) pe.motifMean(k, :);
    sigLabel = sprintf('motif %d DA', k);
    sigTag   = sprintf('motif%02d', k);
else
    assert(strcmpi(opt.signal, 'global'), 'signal must be ''global'' or a motif number.');
    getTrace = @(pe) pe.globalMean;
    sigLabel = 'global DA';
    sigTag   = 'global';
end

firstR = rezC{S.a(1), S.s(1)};
tint = firstR.meta.tint;
if isnumeric(opt.signal)
    K = size(firstR.peth.(trialTypes{1}).motifMean, 1);
    assert(opt.signal >= 1 && opt.signal <= K, 'signal must be in 1..%d.', K);
end
normTag = '';
if isfield(firstR.meta, 'normMode') && startsWith(firstR.meta.normMode, 'zscore')
    normTag = ' (z)';
end

%% ---- bin means ----
nT = numel(tint);
out = struct();
out.binSessions = cellfun(@(ix) S.header(ix), binIdx, 'UniformOutput', false)';

for i = 1:numel(trialTypes)
    tt = trialTypes{i};
    M = nan(3, nT); E = nan(3, nT);
    for b = 1:3
        ix = binIdx{b};
        if isempty(ix), continue; end
        X = nan(numel(ix), nT);
        for j = 1:numel(ix)
            r = rezC{S.a(ix(j)), S.s(ix(j))};
            if isfield(r.peth, tt) && r.peth.(tt).n > 0
                X(j, :) = getTrace(r.peth.(tt));
            end
        end
        ok = any(isfinite(X), 2);
        X = X(ok, :);
        if isempty(X), continue; end
        M(b, :) = mean(X, 1, 'omitnan');
        if size(X, 1) > 1
            E(b, :) = std(X, 0, 1, 'omitnan') ./ sqrt(size(X, 1));
        end
    end
    out.binMeans.(tt) = M;
    out.binSem.(tt)   = E;
end

%% ---- plot ----
nP = numel(trialTypes);
fig = figure('Color', 'w', 'Position', [100 100 440*nP 430]);
ax = gobjects(1, nP);

for i = 1:nP
    tt = trialTypes{i};
    ax(i) = subplot(1, nP, i);
    hold(ax(i), 'on');
    M = out.binMeans.(tt); E = out.binSem.(tt);

    hL = gobjects(0); lbl = {};
    for b = 1:3
        if all(~isfinite(M(b, :))), continue; end
        c = opt.colors(b, :);
        if opt.showSEM && any(isfinite(E(b, :)))
            ok = isfinite(M(b,:)) & isfinite(E(b,:));
            fill(ax(i), [tint(ok), fliplr(tint(ok))], ...
                [M(b,ok) - E(b,ok), fliplr(M(b,ok) + E(b,ok))], c, ...
                'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        hL(end+1) = plot(ax(i), tint, M(b, :), '-', 'Color', c, 'LineWidth', opt.lineWidth); %#ok<AGROW>
        lbl{end+1} = binNames{b}; %#ok<AGROW>
    end

    xlabel(ax(i), 'Time (s)', 'FontWeight', 'bold');
    ylabel(ax(i), sprintf('%s%s (%s)', sigLabel, normTag, upperFirst(tt)), 'FontWeight', 'bold');
    set(ax(i), 'TickDir', 'out', 'Box', 'on', 'LineWidth', 1, 'FontSize', 11);
    if ~isempty(opt.xLim), xlim(ax(i), opt.xLim); else, xlim(ax(i), [tint(1) tint(end)]); end

    if i == 1 && ~isempty(hL)
        lg = legend(ax(i), hL, lbl, 'Location', 'north', 'Box', 'off', ...
            'FontAngle', 'italic', 'FontWeight', 'bold');
        lg.ItemTokenSize = [18 18];
    end
end

% shared y-limits: the bins are only comparable across panels on one axis
if ~isempty(opt.yLim)
    set(ax, 'YLim', opt.yLim);
else
    yl = cell2mat(arrayfun(@(h) ylim(h), ax(:), 'UniformOutput', false));
    set(ax, 'YLim', [min(yl(:,1)), max(yl(:,2))]);
end

% dashed epoch boundaries, drawn after y-limits are final
if opt.showEpochs
    bounds = unique([opt.epochs.t]);
    for i = 1:nP
        for bx = bounds
            xline(ax(i), bx, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1, ...
                'HandleVisibility', 'off');
        end
        hold(ax(i), 'off');
    end
    drawnow;
    for i = 1:nP
        local_epochBar(fig, ax(i), opt.epochs);
    end
end

sgtitle(fig, sprintf('%s | %s | early / intermediate / late relative to day 4', mouseId, sigLabel), ...
    'FontWeight', 'bold', 'FontSize', 11);

out.fig = fig; out.ax = ax; out.tint = tint; out.cutoff = cutoff;

%% ---- save ----
if strlength(strtrim(string(opt.figSaveDir))) > 0
    d = char(string(opt.figSaveDir));
    if exist(d, 'dir') ~= 7, mkdir(d); end
    f = fullfile(d, sprintf('DAPethBins_%s_%s_%s.pdf', mouseId, sigTag, strjoin(trialTypes, '')));
    set(fig, 'InvertHardcopy', 'off');
    print(fig, f, '-dpdf', '-painters', '-bestfit');
    fprintf('Saved:\n  %s\n', f);
end
end

%% ========================================================================
function local_epochBar(fig, ax, epochs)
% Thin bar above the axes: filled epochs dark grey with white italic text,
% unfilled ones outlined with grey italic text -- as in the reference figure.
pos = ax.Position;
h = 0.05;
axB = axes(fig, 'Position', [pos(1), pos(2) + pos(4) + 0.012, pos(3), h]);
hold(axB, 'on');
xl = xlim(ax);
grey = [0.45 0.45 0.45];
for e = 1:numel(epochs)
    t = epochs(e).t;
    t = [max(t(1), xl(1)), min(t(2), xl(2))];
    if t(2) <= t(1), continue; end
    if epochs(e).filled
        rectangle(axB, 'Position', [t(1) 0 diff(t) 1], 'FaceColor', grey, 'EdgeColor', grey, 'LineWidth', 1);
        txtCol = [1 1 1];
    else
        rectangle(axB, 'Position', [t(1) 0 diff(t) 1], 'FaceColor', 'w', 'EdgeColor', grey, 'LineWidth', 1);
        txtCol = [0.6 0.6 0.6];
    end
    text(axB, mean(t), 0.5, epochs(e).name, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'Color', txtCol, 'FontAngle', 'italic', ...
        'FontWeight', 'bold', 'FontSize', 11);
end
xlim(axB, xl); ylim(axB, [0 1]);
axis(axB, 'off');
end

function s = upperFirst(s)
s = char(s);
if strcmpi(s, 'cr') || strcmpi(s, 'fa'), s = upper(s); return; end
s(1) = upper(s(1));
end