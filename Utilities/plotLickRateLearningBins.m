function out = plotLickRateLearningBins(rezC, sessT, mouseId, varargin)
%PLOTLICKRATELEARNINGBINS
%   One animal's mean lick-rate function, split into the same three learning
%   bins as plotDAPethLearningBins, one panel per trial type. The companion
%   behavioural figure to the DA PETHs: same animal, same sessions, same
%   bins, so a DA change can be read against the licking that accompanies it.
%
%   BINS come from sessT.learnBin when present (written by
%   batch_DAPeth_globalAndMotifWeighted), so this figure and the DA figure
%   cannot disagree about which session is early/intermediate/late. Only if
%   that column is missing are the bins re-derived from day4MarkC:
%     early        = sessions before the animal's day-4 date
%     intermediate = first half of the day-4-onward sessions
%     late         = second half (odd count -> extra to late)
%
%   THE TRACE is the combined lick rate built in the batch: video licks
%   before lickSplit (default 2.3 s), contact licks after. The detector
%   changes at that boundary, so a step there can reflect detector
%   sensitivity rather than behaviour -- it is marked with a dotted line
%   rather than left for the reader to rediscover.
%
%   AVERAGING: each bin's trace is the mean of its sessions' mean lick
%   rates -- sessions weighted equally, matching the DA figure and the
%   two-stage convention used throughout this project.
%
%   out = plotLickRateLearningBins(rezC, sessT, mouseId, ...)
%
% NAME-VALUE
%   'trialTypes'  : {'cr','hit'} -- one panel each, in this order. Use
%                   {'hit','cr','miss','fa'} for all four.
%   'day4MarkC'   : only needed if sessT has no learnBin column
%   'showSEM'     : false -- shade +/- SEM ACROSS SESSIONS within a bin
%   'colors'      : [3 x 3] early/intermediate/late; default pale -> saturated blue
%   'lineWidth'   : 2.5
%   'yLim'        : [] (shared auto) or [lo hi]
%   'xLim'        : []
%   'showEpochs'  : true -- tone/response bar and dashed boundaries
%   'showSplit'   : true -- dotted line at the video/contact boundary
%   'epochs'      : struct array (.name, .t, .filled)
%   'figSaveDir'  : ''
%
% OUTPUT
%   .binMeans.<type> [3 x T], .binSem.<type> [3 x T], .binSessions {3x1},
%   .tint, .fig, .ax

%% ---- parse ----
p = inputParser;
p.addParameter('trialTypes', {'cr','hit'}, @(c) iscellstr(c) || isstring(c));
p.addParameter('day4MarkC', {}, @(c) isempty(c) || iscell(c) || isstring(c));
p.addParameter('showSEM', false, @(x) islogical(x) && isscalar(x));
p.addParameter('colors', [0.80 0.80 0.91; 0.42 0.42 0.96; 0.05 0.05 0.86], ...
    @(x) isnumeric(x) && isequal(size(x), [3 3]));
p.addParameter('lineWidth', 2.5, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('yLim', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2 && x(2)>x(1)));
p.addParameter('xLim', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2 && x(2)>x(1)));
p.addParameter('showEpochs', true, @(x) islogical(x) && isscalar(x));
p.addParameter('showSplit', true, @(x) islogical(x) && isscalar(x));
p.addParameter('epochs', struct('name', {'tone','response'}, 't', {[0 2], [2 4]}, ...
    'filled', {true, false}), @isstruct);
p.addParameter('figSaveDir', '', @(s) ischar(s) || isstring(s));
p.parse(varargin{:});
opt = p.Results;

mouseId = char(string(mouseId));
trialTypes = cellstr(opt.trialTypes);
binNames = {'early','intermediate','late'};

%% ---- this animal's sessions, sorted by date ----
rows = find(string(sessT.animal) == string(mouseId));
assert(~isempty(rows), 'mouseId %s not found in sessT.', mouseId);
S = sortrows(sessT(rows, :), 'date');

%% ---- bins ----
if ismember('learnBin', S.Properties.VariableNames) && any(strlength(string(S.learnBin)) > 0)
    lb = string(S.learnBin);
    binIdx = arrayfun(@(b) find(lb == binNames{b}), 1:3, 'UniformOutput', false);
    binSrc = 'sessT.learnBin';
else
    assert(~isempty(opt.day4MarkC), ...
        ['sessT has no learnBin column, so day4MarkC is required. ' ...
         '(Rerunning the batch would store learnBin and remove the need.)']);
    d4 = opt.day4MarkC; if isstring(d4), d4 = cellstr(d4); end
    hit = find(string(d4(:,1)) == string(mouseId), 1, 'first');
    assert(~isempty(hit), 'mouseId %s has no row in day4MarkC.', mouseId);
    cutoff = datetime(char(string(d4{hit,2})), 'InputFormat', 'MMddyy');
    isD4 = dateshift(S.date, 'start', 'day') >= cutoff;
    iD4 = find(isD4); nInt = floor(numel(iD4)/2);
    binIdx = {find(~isD4), iD4(1:nInt), iD4(nInt+1:end)};
    binSrc = 'day4MarkC (recomputed)';
end

fprintf('%s | bins from %s | early %d, intermediate %d, late %d sessions\n', ...
    mouseId, binSrc, numel(binIdx{1}), numel(binIdx{2}), numel(binIdx{3}));
for b = 1:3
    if isempty(binIdx{b})
        warning('LickBins:EmptyBin', '%s: the %s bin has no sessions.', mouseId, binNames{b});
    end
end

%% ---- gather ----
r0 = rezC{S.a(1), S.s(1)};
assert(isfield(r0, 'lick'), ...
    ['rezC has no lick field -- rerun batch_DAPeth_globalAndMotifWeighted with ' ...
     'the lick source available.']);
tint = r0.meta.tint;
nT = numel(tint);
lickSplit = NaN;
if isfield(r0.meta, 'lickSplit'), lickSplit = r0.meta.lickSplit; end

out = struct();
out.binSessions = cellfun(@(ix) S.header(ix), binIdx, 'UniformOutput', false)';
nMissing = 0;

for i = 1:numel(trialTypes)
    tt = trialTypes{i};
    M = nan(3, nT); E = nan(3, nT);
    for b = 1:3
        ix = binIdx{b};
        if isempty(ix), continue; end
        X = nan(numel(ix), nT);
        for j = 1:numel(ix)
            r = rezC{S.a(ix(j)), S.s(ix(j))};
            if ~isfield(r, 'lick') || ~isfield(r.lick, tt)
                nMissing = nMissing + 1; continue;
            end
            if r.lick.(tt).n > 0
                X(j, :) = r.lick.(tt).rateMean;
            end
        end
        X = X(any(isfinite(X), 2), :);
        if isempty(X), continue; end
        M(b, :) = mean(X, 1, 'omitnan');
        if size(X,1) > 1
            E(b, :) = std(X, 0, 1, 'omitnan') ./ sqrt(size(X,1));
        end
    end
    out.binMeans.(tt) = M;
    out.binSem.(tt)   = E;
end
if nMissing > 0
    fprintf('  %d session x type combination(s) had no lick data and were skipped.\n', nMissing);
end

%% ---- plot ----
nP = numel(trialTypes);
fig = figure('Color','w', 'Position', [100 100 440*nP 430]);
ax = gobjects(1, nP);

for i = 1:nP
    tt = trialTypes{i};
    ax(i) = subplot(1, nP, i); hold(ax(i), 'on');
    M = out.binMeans.(tt); E = out.binSem.(tt);

    hL = gobjects(0); lbl = {};
    for b = 1:3
        if all(~isfinite(M(b,:))), continue; end
        c = opt.colors(b,:);
        if opt.showSEM && any(isfinite(E(b,:)))
            ok = isfinite(M(b,:)) & isfinite(E(b,:));
            fill(ax(i), [tint(ok), fliplr(tint(ok))], ...
                [M(b,ok)-E(b,ok), fliplr(M(b,ok)+E(b,ok))], c, ...
                'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        hL(end+1) = plot(ax(i), tint, M(b,:), '-', 'Color', c, 'LineWidth', opt.lineWidth); %#ok<AGROW>
        lbl{end+1} = binNames{b}; %#ok<AGROW>
    end

    xlabel(ax(i), 'Time (s)', 'FontWeight', 'bold');
    ylabel(ax(i), sprintf('lick rate (licks/s) (%s)', upperType(tt)), 'FontWeight', 'bold');
    set(ax(i), 'TickDir', 'out', 'Box', 'on', 'LineWidth', 1, 'FontSize', 11);
    if ~isempty(opt.xLim), xlim(ax(i), opt.xLim); else, xlim(ax(i), [tint(1) tint(end)]); end
    if i == 1 && ~isempty(hL)
        lg = legend(ax(i), hL, lbl, 'Location', 'northwest', 'Box', 'off', ...
            'FontAngle', 'italic', 'FontWeight', 'bold');
        lg.ItemTokenSize = [18 18];
    end
end

if ~isempty(opt.yLim)
    set(ax, 'YLim', opt.yLim);
else
    yl = cell2mat(arrayfun(@(h) ylim(h), ax(:), 'UniformOutput', false));
    set(ax, 'YLim', [min(yl(:,1)), max(yl(:,2))]);
end

for i = 1:nP
    if opt.showEpochs
        for bx = unique([opt.epochs.t])
            xline(ax(i), bx, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1, 'HandleVisibility', 'off');
        end
    end
    % detector boundary: a step here can be sensitivity, not behaviour
    if opt.showSplit && isfinite(lickSplit)
        xline(ax(i), lickSplit, ':', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.2, ...
            'HandleVisibility', 'off');
    end
    hold(ax(i), 'off');
end

if opt.showEpochs
    drawnow;
    for i = 1:nP, local_epochBar(fig, ax(i), opt.epochs); end
end

splitNote = '';
if opt.showSplit && isfinite(lickSplit)
    splitNote = sprintf(' | dotted: video/contact split at %.2f s', lickSplit);
end
sgtitle(fig, sprintf('%s | lick rate | early / intermediate / late relative to day 4%s', ...
    mouseId, splitNote), 'FontWeight', 'bold', 'FontSize', 11);

out.tint = tint; out.fig = fig; out.ax = ax;

%% ---- save ----
if strlength(strtrim(string(opt.figSaveDir))) > 0
    d = char(string(opt.figSaveDir));
    if exist(d, 'dir') ~= 7, mkdir(d); end
    f = fullfile(d, sprintf('lickRateBins_%s_%s.pdf', mouseId, strjoin(trialTypes, '')));
    set(fig, 'InvertHardcopy', 'off');
    print(fig, f, '-dpdf', '-painters', '-bestfit');
    fprintf('Saved:\n  %s\n', f);
end
end

%% ========================================================================
function local_epochBar(fig, ax, epochs)
pos = ax.Position;
axB = axes(fig, 'Position', [pos(1), pos(2) + pos(4) + 0.012, pos(3), 0.05]);
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
xlim(axB, xl); ylim(axB, [0 1]); axis(axB, 'off');
end

function s = upperType(tt)
switch lower(tt)
    case 'cr', s = 'CR';
    case 'fa', s = 'FA';
    otherwise, s = [upper(tt(1)) tt(2:end)];
end
end