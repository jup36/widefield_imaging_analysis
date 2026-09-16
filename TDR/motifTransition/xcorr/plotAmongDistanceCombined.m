function out = plotAmongDistanceCombined(amongOuts, streamLabels, varargin)
%PLOTAMONGDISTANCECOMBINED
%   Combines two or more runAmongDistancePermutationTest outputs (e.g.
%   Hit and CR) into ONE trajectory figure and ONE difference-vs-null
%   figure. Draws entirely from the stored amongOut structs -- no
%   permutations are rerun, so this is instant and always consistent
%   with the stats you already have.
%
%   out = plotAmongDistanceCombined({amongOut_hit, amongOut_cr}, {'Hit','CR'}, ...)
%
% LAYOUT ('layout')
%   'panels'  (default) -- one panel per stream, side by side, with a
%              SHARED y-axis so the streams are directly comparable.
%              Recommended: each panel keeps its own null band and its
%              own significance markers, so nothing is lost.
%   'overlay' -- all streams on one axis. Streams are distinguished by
%              line style (solid / dashed / dotted) and marker shape;
%              groups keep their colours. For the difference figure the
%              null bands are drawn as light outlines rather than fills,
%              because two filled bands over each other are unreadable.
%              Use this when the point is the streams' relative
%              magnitude, not the per-stream null.
%
% NAME-VALUE ARGS
%   'layout'        : 'panels' (default) | 'overlay'
%   'yLimTraj'      : [] (auto, shared across panels) | [lo hi]
%   'yLimDiff'      : [] (auto, shared across panels) | [lo hi]
%   'alpha'         : 0.05 -- marker threshold; should match the runs
%   'groupColors'   : [2 x 3] RGB for groups 1 and 2. Default: MATLAB
%                     lines(2), i.e. the same blue/orange as the
%                     single-stream figures.
%   'printFig'      : false
%   'figSaveDir'    : required if printFig
%   'figureNameBase': 'amongDistanceCombined'
%
% OUTPUT (out)
%   .figTraj, .figDiff : figure handles
%   .axTraj,  .axDiff  : axes handles (one per panel, or one for overlay)
%
% REQUIREMENTS on each amongOut: fields amongX, amongY_group1,
%   amongY_group2, groupNames, obsPerPos, nullPerPos, pPerPos_uncorr,
%   pPerPos_maxStat, pValue, contrastStat, amongMetric, day4Filter. All
%   runs must have the same nSessAnalyze -- the window is what makes the
%   x-axes comparable -- and the same groupNames order.

%% -------------------- parse --------------------
p = inputParser;
p.addParameter('layout', 'panels', @(s) any(strcmpi(string(s), ["panels","overlay"])));
p.addParameter('yLimTraj', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2 && x(2)>x(1)));
p.addParameter('yLimDiff', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2 && x(2)>x(1)));
p.addParameter('alpha', 0.05, @(x) isnumeric(x) && isscalar(x) && x>0 && x<1);
p.addParameter('groupColors', flipud(lines(2)), @(x) isnumeric(x) && isequal(size(x), [2 3]));
p.addParameter('printFig', false, @(x) islogical(x) && isscalar(x));
p.addParameter('figSaveDir', '', @(s) ischar(s) || isstring(s));
p.addParameter('figureNameBase', 'amongDistanceCombined', @(s) ischar(s) || isstring(s));
p.parse(varargin{:});
opt = p.Results;
layout = lower(string(opt.layout));

%% -------------------- validate --------------------
assert(iscell(amongOuts) && numel(amongOuts) >= 2, 'Pass a cell array of at least two amongOut structs.');
assert(numel(streamLabels) == numel(amongOuts), 'One label per amongOut.');
nS = numel(amongOuts);

nPos = amongOuts{1}.nSessAnalyze;
gNames = amongOuts{1}.groupNames;
for s = 2:nS
    assert(amongOuts{s}.nSessAnalyze == nPos, ...
        'Stream "%s" has nSessAnalyze=%d but "%s" has %d -- windows must match to share an x-axis.', ...
        streamLabels{s}, amongOuts{s}.nSessAnalyze, streamLabels{1}, nPos);
    assert(isequal(amongOuts{s}.groupNames, gNames), ...
        'groupNames differ between streams -- runs must use the same groupDefs.');
end
xPos = (1:nPos)';

% Consistent day-4 state across streams, for the title
d4 = cellfun(@(a) a.day4Filter.applied, amongOuts);
if all(d4),      d4Str = 'day-4 sessions only';
elseif ~any(d4), d4Str = 'all sessions';
else,            d4Str = 'MIXED day-4 filtering (check inputs)';
end

styleLine   = {'-', '--', ':', '-.'};
styleMarker = {'o', 's', '^', 'd'};
styleLine   = styleLine(mod(0:nS-1, 4) + 1);
styleMarker = styleMarker(mod(0:nS-1, 4) + 1);

%% ======================================================================
%  FIGURE 1 -- trajectories
%% ======================================================================
figTraj = figure('Color','w');
if layout == "panels"
    figTraj.Position(3) = figTraj.Position(3) * (0.55 * nS + 0.45);
    tl = tiledlayout(figTraj, 1, nS, 'TileSpacing','compact', 'Padding','compact');
    axTraj = gobjects(1, nS);

    for s = 1:nS
        A = amongOuts{s};
        axTraj(s) = nexttile(tl, s);
        hold(axTraj(s), 'on');
        plot(axTraj(s), xPos, A.amongY_group1, '-o', 'Color', opt.groupColors(1,:), ...
            'LineWidth', 2, 'MarkerFaceColor', opt.groupColors(1,:), 'DisplayName', gNames{1});
        plot(axTraj(s), xPos, A.amongY_group2, '-o', 'Color', opt.groupColors(2,:), ...
            'LineWidth', 2, 'MarkerFaceColor', opt.groupColors(2,:), 'DisplayName', gNames{2});
        local_markSig(axTraj(s), xPos, [A.amongY_group1; A.amongY_group2], ...
            A.pPerPos_maxStat, A.pPerPos_uncorr, opt.alpha, 1);
        title(axTraj(s), sprintf('%s   (%s p = %.3g)', streamLabels{s}, A.contrastStat, A.pValue));
        xlabel(axTraj(s), 'Session position (right-aligned)');
        if s == 1
            ylabel(axTraj(s), sprintf('Mean pairwise distance among group (%s)', A.amongMetric));
            legend(axTraj(s), 'Location', 'best');
        end
        grid(axTraj(s), 'on'); box(axTraj(s), 'off');
        set(gca, 'TickDir', 'out')
        xticks(axTraj(s), xPos); xlim(axTraj(s), [0.8, nPos + 0.2]);
        hold(axTraj(s), 'off');
    end
    local_shareY(axTraj, opt.yLimTraj);
    sgtitle(tl, sprintf('Among-group dispersion -- %s   |   * FWER p<%.2f, o uncorrected only', ...
        d4Str, opt.alpha));

else  % overlay
    axTraj = axes(figTraj);
    hold(axTraj, 'on');
    for s = 1:nS
        A = amongOuts{s};
        plot(axTraj, xPos, A.amongY_group1, [styleLine{s} styleMarker{s}], ...
            'Color', opt.groupColors(1,:), 'LineWidth', 2, ...
            'MarkerFaceColor', opt.groupColors(1,:), ...
            'DisplayName', sprintf('%s -- %s', gNames{1}, streamLabels{s}));
        plot(axTraj, xPos, A.amongY_group2, [styleLine{s} styleMarker{s}], ...
            'Color', opt.groupColors(2,:), 'LineWidth', 2, ...
            'MarkerFaceColor', opt.groupColors(2,:), ...
            'DisplayName', sprintf('%s -- %s', gNames{2}, streamLabels{s}));
    end
    % markers stacked by stream so they don't collide
    yAll = cell2mat(cellfun(@(a) [a.amongY_group1; a.amongY_group2], amongOuts, 'UniformOutput', false));
    for s = 1:nS
        A = amongOuts{s};
        local_markSig(axTraj, xPos, yAll, A.pPerPos_maxStat, A.pPerPos_uncorr, opt.alpha, s, styleMarker{s});
    end
    pStr = strjoin(arrayfun(@(s) sprintf('%s p=%.3g', streamLabels{s}, amongOuts{s}.pValue), 1:nS, ...
        'UniformOutput', false), ', ');
    title(axTraj, sprintf('Among-group dispersion -- %s\n%s   |   * FWER p<%.2f, o uncorrected only', ...
        d4Str, pStr, opt.alpha));
    xlabel(axTraj, 'Session position (right-aligned; last = each animal''s final session)');
    ylabel(axTraj, sprintf('Mean pairwise distance among group (%s)', amongOuts{1}.amongMetric));
    legend(axTraj, 'Location', 'best');
    set(gca, 'TickDir', 'out')
    grid(axTraj, 'on'); box(axTraj, 'off');
    xticks(axTraj, xPos); xlim(axTraj, [0.8, nPos + 0.2]);
    if ~isempty(opt.yLimTraj), ylim(axTraj, opt.yLimTraj); end
    hold(axTraj, 'off');
end

%% ======================================================================
%  FIGURE 2 -- difference vs. null
%% ======================================================================
obsColor = [0.85 0.20 0.10];
figDiff = figure('Color','w');
if layout == "panels"
    figDiff.Position(3) = figDiff.Position(3) * (0.55 * nS + 0.45);
    tl2 = tiledlayout(figDiff, 1, nS, 'TileSpacing','compact', 'Padding','compact');
    axDiff = gobjects(1, nS);

    for s = 1:nS
        A = amongOuts{s};
        [lo, hi, med, thrPos, thrFwer] = local_nullBands(A, opt.alpha);
        axDiff(s) = nexttile(tl2, s);
        hold(axDiff(s), 'on');
        fill(axDiff(s), [xPos; flipud(xPos)], [lo; flipud(hi)], [0.6 0.6 0.6], ...
            'FaceAlpha', 0.25, 'EdgeColor', 'none', 'DisplayName', 'null 2.5-97.5%');
        plot(axDiff(s), xPos, med, '-', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.2, 'DisplayName', 'null median');
        plot(axDiff(s), xPos, thrPos, ':', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.4, ...
            'DisplayName', sprintf('null %g%% (per position)', 100*(1-opt.alpha)));
        yline(axDiff(s), thrFwer, '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.4, ...
            'DisplayName', sprintf('max-stat %g%% (FWER)', 100*(1-opt.alpha)));
        yline(axDiff(s), 0, '-', 'Color', [0.75 0.75 0.75], 'HandleVisibility', 'off');
        plot(axDiff(s), xPos, A.obsPerPos, '-o', 'Color', obsColor, 'LineWidth', 2.2, ...
            'MarkerFaceColor', obsColor, 'MarkerSize', 7, ...
            'DisplayName', sprintf('observed (%s - %s)', gNames{2}, gNames{1}));
        local_markSig(axDiff(s), xPos, [A.obsPerPos; lo; hi], A.pPerPos_maxStat, A.pPerPos_uncorr, opt.alpha, 1);
        title(axDiff(s), sprintf('%s   (%s p = %.3g)', streamLabels{s}, A.contrastStat, A.pValue));
        xlabel(axDiff(s), 'Session position (right-aligned)');
        if s == 1
            ylabel(axDiff(s), sprintf('Among-distance difference, %s - %s', gNames{2}, gNames{1}));
            legend(axDiff(s), 'Location', 'best');
        end
        grid(axDiff(s), 'on'); box(axDiff(s), 'off');
        xticks(axDiff(s), xPos); xlim(axDiff(s), [0.8, nPos + 0.2]);
        set(gca, 'TickDir', 'out')
        hold(axDiff(s), 'off');
    end
    local_shareY(axDiff, opt.yLimDiff);
    sgtitle(tl2, sprintf('Observed group difference vs. permutation null -- %s   |   * FWER p<%.2f, o uncorrected only', ...
        d4Str, opt.alpha));

else  % overlay: outlines instead of fills, one observed curve per stream
    axDiff = axes(figDiff);
    hold(axDiff, 'on');
    streamColors = lines(nS + 2); streamColors = streamColors(3:end, :);   % avoid clashing with group colours
    yAllD = [];
    for s = 1:nS
        A = amongOuts{s};
        [lo, hi, ~, ~, thrFwer] = local_nullBands(A, opt.alpha);
        c = streamColors(s, :);
        plot(axDiff, xPos, lo, styleLine{s}, 'Color', [c 0.45], 'LineWidth', 0.9, 'HandleVisibility', 'off');
        plot(axDiff, xPos, hi, styleLine{s}, 'Color', [c 0.45], 'LineWidth', 0.9, ...
            'DisplayName', sprintf('%s null 2.5-97.5%%', streamLabels{s}));
        yline(axDiff, thrFwer, styleLine{s}, 'Color', c, 'LineWidth', 1.3, ...
            'DisplayName', sprintf('%s max-stat %g%%', streamLabels{s}, 100*(1-opt.alpha)));
        plot(axDiff, xPos, A.obsPerPos, [styleLine{s} styleMarker{s}], 'Color', c, 'LineWidth', 2.2, ...
            'MarkerFaceColor', c, 'MarkerSize', 7, ...
            'DisplayName', sprintf('%s observed (%s - %s)', streamLabels{s}, gNames{2}, gNames{1}));
        yAllD = [yAllD; A.obsPerPos; lo; hi]; %#ok<AGROW>
    end
    yline(axDiff, 0, '-', 'Color', [0.75 0.75 0.75], 'HandleVisibility', 'off');
    for s = 1:nS
        A = amongOuts{s};
        local_markSig(axDiff, xPos, yAllD, A.pPerPos_maxStat, A.pPerPos_uncorr, opt.alpha, s, styleMarker{s});
    end
    pStr = strjoin(arrayfun(@(s) sprintf('%s p=%.3g', streamLabels{s}, amongOuts{s}.pValue), 1:nS, ...
        'UniformOutput', false), ', ');
    title(axDiff, sprintf('Observed group difference vs. permutation null -- %s\n%s   |   * FWER p<%.2f, o uncorrected only', ...
        d4Str, pStr, opt.alpha));
    xlabel(axDiff, 'Session position (right-aligned; last = each animal''s final session)');
    ylabel(axDiff, sprintf('Among-distance difference, %s - %s', gNames{2}, gNames{1}));
    legend(axDiff, 'Location', 'best');
    grid(axDiff, 'on'); box(axDiff, 'off');
    xticks(axDiff, xPos); xlim(axDiff, [0.8, nPos + 0.2]);
    if ~isempty(opt.yLimDiff), ylim(axDiff, opt.yLimDiff); end
    set(gca, 'TickDir', 'out')
    hold(axDiff, 'off');
end

%% -------------------- save --------------------
if opt.printFig
    assert(strlength(strtrim(string(opt.figSaveDir))) > 0, 'printFig is true but no figSaveDir given.');
    outDir = char(string(opt.figSaveDir));
    if exist(outDir, 'dir') ~= 7, mkdir(outDir); end
    tag = sprintf('%s_%s_%s_%s', char(string(opt.figureNameBase)), ...
        strjoin(cellfun(@(x) char(string(x)), streamLabels, 'UniformOutput', false), ''), ...
        char(layout), char(datetime('today','Format','MMddyy')));
    for pair = {{figTraj, 'trajectory'}, {figDiff, 'diffVsNull'}}
        fH = pair{1}{1};
        f  = fullfile(outDir, sprintf('%s_%s.pdf', tag, pair{1}{2}));
        set(fH, 'InvertHardcopy', 'off');
        print(fH, f, '-dpdf', '-painters', '-bestfit');
        fprintf('Saved figure:\n  %s\n', f);
    end
end

out = struct('figTraj', figTraj, 'figDiff', figDiff, 'axTraj', axTraj, 'axDiff', axDiff);
end

%% ========================= local helpers =========================
function [lo, hi, med, thrPos, thrFwer] = local_nullBands(A, alpha)
lo      = prctile(A.nullPerPos, 2.5, 1)';
hi      = prctile(A.nullPerPos, 97.5, 1)';
med     = median(A.nullPerPos, 1, 'omitnan')';
thrPos  = prctile(A.nullPerPos, 100*(1-alpha), 1)';
nullMax = max(A.nullPerPos, [], 2, 'omitnan');
thrFwer = prctile(nullMax, 100*(1-alpha));
end


function local_markSig(ax, xPos, yRef, pFwer, pUnc, alpha, row, mk)
% Significance markers above the data. 'row' stacks markers for multiple
% streams on one axis; 'mk' overrides the uncorrected-only marker shape.
if nargin < 8, mk = 'o'; end
yTop = max(yRef, [], 'omitnan');
yPad = 0.05 * range(yRef) * row;
for pp = 1:numel(xPos)
    if pFwer(pp) < alpha
        plot(ax, xPos(pp), yTop + yPad, 'k*', 'MarkerSize', 9, 'HandleVisibility', 'off');
    elseif pUnc(pp) < alpha
        plot(ax, xPos(pp), yTop + yPad, ['k' mk], 'MarkerSize', 6, 'HandleVisibility', 'off');
    end
end
end


function local_shareY(axArr, yLimFixed)
% Common y-limits across panels: either the user's, or the union of the
% auto limits so no panel is clipped.
if ~isempty(yLimFixed)
    set(axArr, 'YLim', yLimFixed);
else
    yl = cell2mat(arrayfun(@(h) ylim(h), axArr(:), 'UniformOutput', false));
    set(axArr, 'YLim', [min(yl(:,1)), max(yl(:,2))]);
end
end