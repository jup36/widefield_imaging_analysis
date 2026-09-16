function out = plotTrialTypeHPeths(tbytDat_hAligned, trI, motifK, varargin)
%PLOTTRIALTYPEHPETHS
%   Trial-type PETHs of one motif's time-corrected H, aligned to tone
%   onset. One trace per trial type (Hit / CR / Miss / FA), mean +/- SEM
%   across trials. Companion to plotTrialTypeGlobalDAPeths, with the same
%   extraction rules, palette, and output structure so the two figures
%   are directly comparable.
%
%   out = plotTrialTypeHPeths(tbytDat_hAligned, trI, motifK, ...)
%
% INPUTS
%   tbytDat_hAligned : {2 x nTrials} cell. Row 1 = [K x T] motif H for
%                      that trial; row 2 = [1 x T] timestamps (sec).
%   trI              : struct of logical trial masks (hitI, crI, missI, faI).
%   motifK           : which motif row to plot.
%
% NAME-VALUE ARGS
%   'trialTypes', 'showSEM', 'minTrials', 'colors', 'ax', 'titleStr' --
%   identical to plotTrialTypeGlobalDAPeths.
%
% OUTPUT (out)
%   .timeX, .peth.<type> ([nTrials x T]), .nTrials.<type>, .fig, .ax
%
% LENGTH HANDLING
%   Trials whose H matrix is empty, has fewer than motifK rows, or whose
%   time dimension does not match the shared time axis are DROPPED and
%   counted. The previous version concatenated rows directly and errored
%   on the first length mismatch; dropping with a report is the safer
%   behaviour, but a large drop count means the alignment upstream needs
%   looking at, not that this function has handled it.

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

assert(isnumeric(motifK) && isscalar(motifK) && motifK >= 1 && motifK == round(motifK), ...
    'motifK must be a positive integer.');
trialTypes = cellstr(opt.trialTypes);

palette = containers.Map( ...
    {'Hit', 'CR', 'Miss', 'FA'}, ...
    {[0.85 0.33 0.10], [0.00 0.45 0.74], [0.93 0.69 0.13], [0.49 0.18 0.56]});
maskField = containers.Map( ...
    {'Hit', 'CR', 'Miss', 'FA'}, ...
    {'hitI', 'crI', 'missI', 'faI'});

%% -------------------- time axis --------------------
timeX = local_firstNonEmptyTime(tbytDat_hAligned);
assert(~isempty(timeX), 'No valid time vector found in row 2 of tbytDat_hAligned.');
timeX = double(timeX(:)');
T = numel(timeX);

% Report the length distribution once, so a mismatch problem is visible
% up front rather than inferred from the drop count.
lens = cellfun(@(x) size(x, 2), tbytDat_hAligned(1, :));
lens = lens(lens > 0);
uLens = unique(lens);
if numel(uLens) > 1
    fprintf('  NOTE: H trial lengths vary -- %s. Using T = %d (from the time axis); others are dropped.\n', ...
        strjoin(arrayfun(@(L) sprintf('%d (x%d)', L, sum(lens == L)), uLens, 'UniformOutput', false), ', '), T);
end

%% -------------------- extract per type --------------------
nTrialsTotal = size(tbytDat_hAligned, 2);
pethS   = struct();
nTrS    = struct();
plotted = {};
nDroppedTotal = 0;

for i = 1:numel(trialTypes)
    tt = trialTypes{i};
    if ~isKey(maskField, tt)
        warning('plotTrialTypeHPeths:unknownType', 'Unknown trial type "%s" -- skipping.', tt);
        continue;
    end
    fld = maskField(tt);
    if ~isfield(trI, fld)
        fprintf('  trI has no field "%s" -- skipping %s.\n', fld, tt);
        continue;
    end

    mask = logical(trI.(fld)(:)');
    assert(numel(mask) == nTrialsTotal, ...
        'trI.%s has %d entries but tbytDat_hAligned has %d trials.', fld, numel(mask), nTrialsTotal);

    [M, nDropped] = local_extractMotifRow(tbytDat_hAligned, mask, motifK, T);
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
    fprintf('  %d trial(s) dropped across types (empty, <%d motif rows, or length != %d).\n', ...
        nDroppedTotal, motifK, T);
end
assert(~isempty(plotted), 'No trial type had enough usable trials to plot.');

%% -------------------- colours --------------------
nP = numel(plotted);
if isempty(opt.colors)
    cols = zeros(nP, 3);
    for i = 1:nP, cols(i,:) = palette(plotted{i}); end
else
    assert(size(opt.colors, 1) >= numel(trialTypes), 'colors needs one row per entry of trialTypes.');
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
ylabel(ax, sprintf('Motif %d H', motifK));
if strlength(string(opt.titleStr)) > 0
    title(ax, char(opt.titleStr));
else
    title(ax, sprintf('Trial-aligned H, motif %d (mean %s SEM)', motifK, ...
        char(ternary(opt.showSEM, '\pm', 'only'))));
end
legend(ax, hLines, 'Location', 'best');
xlim(ax, [timeX(1), timeX(end)]);
set(ax, 'TickDir', 'out');
box(ax, 'off');
hold(ax, 'off');

%% -------------------- output --------------------
out = struct('timeX', timeX, 'motifK', motifK, 'peth', pethS, 'nTrials', nTrS, 'fig', fig, 'ax', ax);
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


function [M, nDropped] = local_extractMotifRow(alignedCell, mask, motifK, T)
% Stack row motifK of each masked trial's [K x T] matrix into [nTrials x T].
% Trials that are empty, have too few rows, or the wrong length are dropped.
idx = find(mask);
M = nan(numel(idx), T);
keep = false(numel(idx), 1);
for ii = 1:numel(idx)
    tr = idx(ii);
    if tr > size(alignedCell, 2), continue; end
    A = alignedCell{1, tr};
    if isempty(A) || size(A, 1) < motifK || size(A, 2) ~= T, continue; end
    M(ii, :) = double(A(motifK, :));
    keep(ii) = true;
end
nDropped = sum(~keep);
M = M(keep, :);
end


function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

% 
% function plotTrialTypeHPeths(tbytDat_hAligned, trI, motifK)
% 
% timeX = getFirstNonEmptyTime(tbytDat_hAligned);
% 
% if isempty(timeX)
%     warning('No valid H aligned time vector found.');
%     return;
% end
% 
% trialTypes = {};
% trialMasks = {};
% 
% if isfield(trI, 'crI')
%     trialTypes{end+1} = 'CR'; %#ok<AGROW>
%     trialMasks{end+1} = trI.crI; %#ok<AGROW>
% end
% 
% if isfield(trI, 'hitI')
%     trialTypes{end+1} = 'Hit'; %#ok<AGROW>
%     trialMasks{end+1} = trI.hitI; %#ok<AGROW>
% end
% 
% if isfield(trI, 'missI')
%     trialTypes{end+1} = 'Miss'; %#ok<AGROW>
%     trialMasks{end+1} = trI.missI; %#ok<AGROW>
% end
% 
% if isfield(trI, 'faI')
%     trialTypes{end+1} = 'FA'; %#ok<AGROW>
%     trialMasks{end+1} = trI.faI; %#ok<AGROW>
% end
% 
% figure;
% hold on;
% 
% for i = 1:numel(trialTypes)
% 
%     trialMask = logical(trialMasks{i});
%     pethMat = extractMotifAlignedMatrix(tbytDat_hAligned, trialMask, motifK);
% 
%     if isempty(pethMat)
%         continue;
%     end
% 
%     plot(timeX, mean(pethMat, 1, 'omitnan'), 'LineWidth', 2);
% end
% 
% xline(0, '--');
% xlabel('Time from tone onset, sec');
% ylabel(sprintf('Motif %d H', motifK));
% title(sprintf('Trial-aligned H, motif %d', motifK));
% legend(trialTypes, 'Location', 'best');
% box off;
% 
% end