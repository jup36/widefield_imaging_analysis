function out = plotTrialTypeGlobalDAPeths(tbytDat_DAglobalAligned, trI, varargin)
%PLOTTRIALTYPEGLOBALDAPETHS
%   Trial-type PETHs of the GLOBAL cortical dopamine transient, aligned to
%   tone onset. One trace per trial type (Hit / CR / Miss / FA), mean +/-
%   SEM across trials.
%
%   Unlike the motif-projected version this replaces, the global DA
%   transient is a single [1 x T] trace per trial, so there is no motif
%   index to select -- the function reads row 1 of every trial cell.
%
%   out = plotTrialTypeGlobalDAPeths(tbytDat_DAglobalAligned, trI, ...)
%
% INPUTS
%   tbytDat_DAglobalAligned : {2 x nTrials} cell. Row 1 = [1 x T] global
%                             DA transient for that trial; row 2 = [1 x T]
%                             timestamps relative to tone onset (sec).
%   trI                     : struct of logical trial masks (hitI, crI,
%                             missI, faI), each [1 x nTrials] or [nTrials x 1].
%
% NAME-VALUE ARGS
%   'trialTypes'  : which types to plot and in what order. Default
%                   {'Hit','CR','Miss','FA'}; types missing from trI are
%                   skipped with a note.
%   'showSEM'     : true (default). Shaded +/- SEM band per trace.
%   'minTrials'   : 3 (default). Skip a type with fewer usable trials.
%   'colors'      : [nTypes x 3] RGB in trialTypes order. Default: a fixed
%                   palette so Hit/CR/Miss/FA colours are stable across
%                   sessions and figures.
%   'ax'          : axes to draw into. Default: new figure.
%   'titleStr'    : override the default title.
%
% OUTPUT (out)
%   .timeX        : [1 x T] time axis
%   .peth         : struct, one field per plotted type, each [nTrials x T]
%   .nTrials      : struct, trial count per plotted type
%   .fig, .ax     : handles
%
% NOTES
%   Trials are taken at face value from trI. Any trial whose DA cell is
%   empty, or whose length disagrees with the shared time axis, is dropped
%   with a count reported -- a silent length mismatch would otherwise
%   corrupt the mean.

%% -------------------- parse --------------------
p = inputParser;
p.addParameter('trialTypes', {'Hit','CR','Miss','FA'}, @(c) iscellstr(c) || isstring(c));
p.addParameter('showSEM', true, @(x) islogical(x) && isscalar(x));
p.addParameter('minTrials', 3, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('colors', [], @(x) isempty(x) || (isnumeric(x) && size(x,2) == 3));
p.addParameter('ax', [], @(h) isempty(h) || isgraphics(h, 'axes'));
p.addParameter('titleStr', '', @(s) ischar(s) || isstring(s));
p.parse(varargin{:});
opt = p.Results;

trialTypes = cellstr(opt.trialTypes);

% Fixed palette keyed by type name so colours mean the same thing in
% every figure, regardless of which types happen to be present.
palette = containers.Map( ...
    {'Hit', 'CR', 'Miss', 'FA'}, ...
    {[0.85 0.33 0.10], [0.00 0.45 0.74], [0.93 0.69 0.13], [0.49 0.18 0.56]});
maskField = containers.Map( ...
    {'Hit', 'CR', 'Miss', 'FA'}, ...
    {'hitI', 'crI', 'missI', 'faI'});

%% -------------------- time axis --------------------
timeX = local_firstNonEmptyTime(tbytDat_DAglobalAligned);
assert(~isempty(timeX), 'No valid time vector found in row 2 of tbytDat_DAglobalAligned.');
timeX = double(timeX(:)');
T = numel(timeX);

%% -------------------- extract per type --------------------
nTrialsTotal = size(tbytDat_DAglobalAligned, 2);
pethS   = struct();
nTrS    = struct();
plotted = {};
nDroppedTotal = 0;

for i = 1:numel(trialTypes)
    tt = trialTypes{i};
    if ~isKey(maskField, tt)
        warning('plotTrialTypeGlobalDAPeths:unknownType', 'Unknown trial type "%s" -- skipping.', tt);
        continue;
    end
    fld = maskField(tt);
    if ~isfield(trI, fld)
        fprintf('  trI has no field "%s" -- skipping %s.\n', fld, tt);
        continue;
    end

    mask = logical(trI.(fld)(:)');
    assert(numel(mask) == nTrialsTotal, ...
        'trI.%s has %d entries but tbytDat_DAglobalAligned has %d trials.', ...
        fld, numel(mask), nTrialsTotal);

    [M, nDropped] = local_extractGlobalDA(tbytDat_DAglobalAligned, mask, T);
    nDroppedTotal = nDroppedTotal + nDropped;

    if size(M, 1) < opt.minTrials
        fprintf('  %s: only %d usable trial(s) (< minTrials=%d) -- skipping.\n', tt, size(M,1), opt.minTrials);
        continue;
    end

    pethS.(tt) = M;
    nTrS.(tt)  = size(M, 1);
    plotted{end+1} = tt; %#ok<AGROW>
end

if nDroppedTotal > 0
    fprintf('  %d trial(s) dropped across types (empty cell or length != %d).\n', nDroppedTotal, T);
end
assert(~isempty(plotted), 'No trial type had enough usable trials to plot.');

%% -------------------- colours --------------------
nP = numel(plotted);
if isempty(opt.colors)
    cols = zeros(nP, 3);
    for i = 1:nP, cols(i,:) = palette(plotted{i}); end
else
    assert(size(opt.colors, 1) >= numel(trialTypes), 'colors needs one row per entry of trialTypes.');
    % map user colours by position in the ORIGINAL trialTypes list
    cols = zeros(nP, 3);
    for i = 1:nP
        cols(i,:) = opt.colors(strcmp(trialTypes, plotted{i}), :);
    end
end

%% -------------------- plot --------------------
if isempty(opt.ax)
    fig = figure('Color', 'w');
    ax  = axes(fig);
else
    ax  = opt.ax;
    fig = ancestor(ax, 'figure');
end
hold(ax, 'on');

hLines = gobjects(1, nP);
for i = 1:nP
    M  = pethS.(plotted{i});
    mu = mean(M, 1, 'omitnan');
    if opt.showSEM
        se = std(M, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(M), 1));
        ok = isfinite(mu) & isfinite(se);
        fill(ax, [timeX(ok), fliplr(timeX(ok))], [mu(ok) - se(ok), fliplr(mu(ok) + se(ok))], ...
            cols(i,:), 'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    hLines(i) = plot(ax, timeX, mu, 'Color', cols(i,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('%s (n=%d)', plotted{i}, nTrS.(plotted{i})));
end

xline(ax, 0, '--', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
yline(ax, 0, '-',  'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off');

xlabel(ax, 'Time from tone onset (s)');
ylabel(ax, 'Global DA (dF/F)');
if strlength(string(opt.titleStr)) > 0
    title(ax, char(opt.titleStr));
else
    title(ax, sprintf('Trial-aligned global DA transient (mean %s SEM)', ...
        char(ternary(opt.showSEM, '\pm', 'only'))));
end
legend(ax, hLines, 'Location', 'best');
xlim(ax, [timeX(1), timeX(end)]);
set(ax, 'TickDir', 'out');
box(ax, 'off');
hold(ax, 'off');

%% -------------------- output --------------------
out = struct('timeX', timeX, 'peth', pethS, 'nTrials', nTrS, 'fig', fig, 'ax', ax);
end

%% ========================= local helpers =========================
function timeX = local_firstNonEmptyTime(alignedCell)
timeX = [];
if size(alignedCell, 1) < 2, return; end
for tr = 1:size(alignedCell, 2)
    if ~isempty(alignedCell{2, tr})
        timeX = alignedCell{2, tr};
        return;
    end
end
end


function [M, nDropped] = local_extractGlobalDA(alignedCell, mask, T)
% Stack row-1 traces for the masked trials into [nTrials x T]. Trials
% whose trace is empty or not length T are dropped and counted.
idx = find(mask);
M = nan(numel(idx), T);
keep = false(numel(idx), 1);
for ii = 1:numel(idx)
    tr = idx(ii);
    if tr > size(alignedCell, 2), continue; end
    x = alignedCell{1, tr};
    if isempty(x), continue; end
    x = double(x(:)');
    if numel(x) ~= T, continue; end
    M(ii, :) = x;
    keep(ii) = true;
end
nDropped = sum(~keep);
M = M(keep, :);
end


function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end