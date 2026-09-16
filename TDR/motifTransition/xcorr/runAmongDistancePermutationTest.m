function amongOut = runAmongDistancePermutationTest(Y, sessInfo, groupDefs, varargin)
%RUNAMONGDISTANCEPERMUTATIONTEST
%   Compares fast vs. slow learners on MEAN PAIRWISE DISTANCE AMONG A
%   GROUP'S OWN MEMBERS, as a right-aligned trajectory across sessions
%   (each animal's final session = position nSessAnalyze, working
%   backward). Unlike the centroid/ellipsoid convergence analyses
%   elsewhere in this project, this metric has NO mechanical self-
%   selection bias: a group's own MEAN is, by construction, the closest
%   possible point to its own members (guaranteed circularity), but
%   pairwise distance BETWEEN actual data points has no such property --
%   picking any random 4-of-9 animals and measuring how far apart they
%   are from each other isn't artificially deflated just because they
%   were selected. So a plain group-label permutation test is valid here
%   with no deeper circularity concern to control for.
%
%   WHY A COMMON RIGHT-ALIGNED WINDOW IS REQUIRED (not just tidiness)
%   ----------------------------------------------------------------
%   The contrast averages each group's among-distance over session
%   positions, so it is tempting to think uneven session counts wash out
%   in that average. They do not. At a deep position (far from the end),
%   only animals with long enough series have a point at all; a group
%   with fewer than 2 members there yields NaN, and 'omitnan' then
%   averages the two groups over DIFFERENT sets of positions. Under
%   permutation this is worse: every relabeling changes which positions
%   are usable, so the null contrasts are not computed on the same
%   footing as the observed one and exchangeability fails. Clipping to
%   nSessAnalyze = the shallowest animal's depth guarantees all members
%   of every possible partition contribute at every position.
%
%   DAY-4 FILTERING ('day4MarkC', optional)
%   ---------------------------------------
%   Restricts the analysis to sessions on or after each animal's own
%   "day 4" date -- the date the full Go/NoGo task with a 1:1 trial ratio
%   began. Pass an {mouseId, 'MMDDYY'} cell array. Sessions before that
%   date are dropped BEFORE any depth/position bookkeeping, so
%   nSessAnalyze is recomputed from the retained sessions only. Because
%   the filter removes an early contiguous block, each animal's retained
%   sessions stay contiguous and end at its true final session, so
%   right-alignment is unaffected.
%
%   PER-POSITION INFERENCE
%   ----------------------
%   The same partitions that build the whole-window null also give a null
%   AT EVERY SESSION POSITION at no extra cost, so the function now
%   reports three things per position:
%     .pPerPos_uncorr : the raw per-position permutation p. Read these
%                       descriptively -- with one test per position they
%                       carry no family-wise control.
%     .pPerPos_maxStat: family-wise corrected across positions by the
%                       max-statistic (a.k.a. tmax) method -- for each
%                       partition, take the LARGEST contrast across all
%                       positions, and compare each observed position
%                       against that single null distribution. This is
%                       the standard permutation FWER control and is
%                       exact here; it needs no independence assumption,
%                       which matters because positions share animals and
%                       are strongly correlated.
%     .pPerPos_BH     : Benjamini-Hochberg FDR across positions, as a
%                       less conservative alternative.
%
%   RESOLUTION LIMIT: with exact enumeration the smallest attainable
%   p is 1/(nPartitions+1) -- 1/127 ~= 0.0079 for a 4-vs-5 split. No
%   per-position p, corrected or not, can go below that. If several
%   positions sit at the floor they are tied, not ranked.
%
%   A NOTE ON WHAT TO PRE-SPECIFY: if the hypothesis is "the groups
%   diverge as animals learn," the sharper test is a single contrast on
%   the TREND rather than six separate positional tests. Pass
%   'contrastStat','slopeDiff' to test whether the between-group gap
%   grows linearly across the window -- one test, full power, no
%   correction needed. The per-position values are then the descriptive
%   follow-up, not the primary result.
%
%   amongOut = runAmongDistancePermutationTest(Y, sessInfo, groupDefs, ...)
%
% INPUTS
%   Y         : [S x dimFull] MDS coordinates for ONE stream.
%   sessInfo  : matching sessInfo table (needs mouseId, sessWithin; also a
%               header-like column if day4MarkC is used).
%   groupDefs : scalar struct, exactly 2 groups (e.g. fast/slow).
%
% NAME-VALUE ARGS
%   'day4MarkC'      : {} (default) or {N x 2} cell of {mouseId,'MMDDYY'}.
%   'sessionDateVar' : '' (auto-detect) or the sessInfo column name.
%   'dims'           : columns of Y to use. Default: 1:size(Y,2).
%   'amongMetric'    : 'euclidean' (default) or 'mahal'.
%   'nSessAnalyze'   : common window depth. Default [] = auto.
%   'contrastStat'   : 'meanDiff' (default) | 'finalDiff' | 'slopeDiff'.
%                      slopeDiff = OLS slope of the per-position group
%                      difference against position; positive means the
%                      gap widens toward the final session.
%   'direction'      : 'slow_minus_fast' (default) or 'fast_minus_slow'.
%   'alphaPerPos'    : 0.05 (default). Used only to mark the plot.
%   'useExact'       : true (default).
%   'maxExact'       : 5000 (default).
%   'nPerm'          : 1000 (default, only if exact is intractable).
%   'permSeed'       : [] (default) or scalar seed.
%   'doPlot'         : true (default). Trajectory figure (both groups'
%                      among-distance curves).
%   'yLimTraj'       : [] (default, automatic) or [lo hi] for the
%                      TRAJECTORY figure's y-axis only. The difference
%                      figure is on an unrelated scale and is unaffected.
%   'groupColors'    : [2 x 3] RGB for groups 1 and 2, in
%                      fieldnames(groupDefs) order. Default flipud(lines(2))
%                      = project convention: fast (group 1) orange,
%                      slow (group 2) blue.
%   'doPlotDiff'     : true (default). Difference figure -- the observed
%                      per-position gap drawn against the permutation
%                      null band, with the same significance markers.
%   'bandPct'        : [2.5 97.5] (default). Percentiles of the null used
%                      for the shaded band on the difference figure.
%   'printFig'       : false (default). Save both figures as vector PDFs.
%   'figSaveDir'     : required when printFig is true.
%   'verbose'        : true (default).
%   'doSave', 'saveDir', 'saveTag' : as in the other stats functions here.
%
% OUTPUT (amongOut)
%   .amongX, .amongY_group1, .amongY_group2, .groupNames
%   .obsContrast, .nullContrast, .pValue            (whole-window)
%   .obsPerPos, .nullPerPos                          (per position)
%   .pPerPos_uncorr, .pPerPos_maxStat, .pPerPos_BH
%   .perPosTable                                     (tidy summary)
%   .pFloor, .nPartitions, .isExact
%   .direction, .contrastStat, .amongMetric, .nSessAnalyze, .day4Filter

%% -------------------- parse options --------------------
p = inputParser;
p.addParameter('day4MarkC', {}, @(x) isempty(x) || (iscell(x) && size(x,2) == 2));
p.addParameter('sessionDateVar', '', @(s) ischar(s) || isstring(s));
p.addParameter('dims', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
p.addParameter('amongMetric', 'euclidean', @(s) any(strcmpi(string(s), ["euclidean","mahal"])));
p.addParameter('nSessAnalyze', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>=2 && x==round(x)));
p.addParameter('contrastStat', 'meanDiff', @(s) any(strcmpi(string(s), ["meanDiff","finalDiff","slopeDiff"])));
p.addParameter('direction', 'slow_minus_fast', @(s) any(strcmpi(string(s), ["slow_minus_fast","fast_minus_slow"])));
p.addParameter('alphaPerPos', 0.05, @(x) isnumeric(x) && isscalar(x) && x>0 && x<1);
p.addParameter('useExact', true, @(x) islogical(x) && isscalar(x));
p.addParameter('maxExact', 5000, @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('nPerm', 1000, @(x) isnumeric(x) && isscalar(x) && x>=1 && x==round(x));
p.addParameter('permSeed', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('doPlot', true, @(x) islogical(x) && isscalar(x));
p.addParameter('yLimTraj', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2 && x(2)>x(1)));
p.addParameter('groupColors', flipud(lines(2)), @(x) isnumeric(x) && isequal(size(x), [2 3]));
p.addParameter('verbose', true, @(x) islogical(x) && isscalar(x));
p.addParameter('doSave', true, @(x) islogical(x) && isscalar(x));
p.addParameter('saveDir', '', @(s) ischar(s) || isstring(s));
p.addParameter('saveTag', '', @(s) ischar(s) || isstring(s));
p.addParameter('doPlotDiff', true, @(x) islogical(x) && isscalar(x));
p.addParameter('bandPct', [2.5 97.5], @(x) isnumeric(x) && numel(x)==2 && x(1)<x(2));
p.addParameter('printFig', false, @(x) islogical(x) && isscalar(x));
p.addParameter('figSaveDir', '', @(s) ischar(s) || isstring(s));
p.parse(varargin{:});
opt = p.Results;

amongMetric  = lower(string(opt.amongMetric));
direction    = lower(string(opt.direction));
contrastStat = lower(string(opt.contrastStat));

if ~isempty(opt.permSeed)
    rng(opt.permSeed);
end

%% -------------------- validate + build fixed mouse universe --------------------
groupNames = fieldnames(groupDefs);
assert(numel(groupNames) == 2, 'runAmongDistancePermutationTest requires exactly 2 groups.');

group1IdC = groupDefs.(groupNames{1});
group2IdC = groupDefs.(groupNames{2});
overlap = intersect(group1IdC, group2IdC);
assert(isempty(overlap), 'Animal(s) %s appear in both groups.', strjoin(overlap, ', '));

allMiceIdC = [group1IdC(:)', group2IdC(:)'];
nAllMice   = numel(allMiceIdC);
n1         = numel(group1IdC);

assert(all(ismember({'mouseId','sessWithin'}, sessInfo.Properties.VariableNames)), ...
    'sessInfo must have mouseId and sessWithin columns.');

mouseId    = string(sessInfo.mouseId);
sessWithin = sessInfo.sessWithin;

missingMice = allMiceIdC(~ismember(string(allMiceIdC), unique(mouseId)));
if ~isempty(missingMice)
    error('The following animal(s) from groupDefs were not found in sessInfo.mouseId: %s', ...
        strjoin(missingMice, ', '));
end

if isempty(opt.dims)
    dims = 1:size(Y,2);
else
    dims = opt.dims(:)';
end
Yd = Y(:, dims);

%% -------------------- day-4 filtering (optional) --------------------
day4Filter = struct('applied', false, 'cutoffDates', [], 'nBefore', [], 'nAfter', [], 'mice', {allMiceIdC});
keepRow = true(numel(mouseId), 1);

if ~isempty(opt.day4MarkC)
    dateVar  = local_findDateVar(sessInfo, opt.sessionDateVar);
    sessDate = local_parseHeaderDate(sessInfo.(dateVar));

    nUnparsed = sum(isnat(sessDate));
    if nUnparsed > 0
        warning('runAmongDistancePermutationTest:unparsedDates', ...
            '%d session(s) in column "%s" could not be date-parsed and will be DROPPED by the day-4 filter.', ...
            nUnparsed, dateVar);
    end

    markIds   = string(opt.day4MarkC(:,1));
    markDates = local_parseMMDDYY(opt.day4MarkC(:,2));
    assert(~any(isnat(markDates)), 'One or more day4MarkC dates could not be parsed as MMDDYY.');

    missingMark = allMiceIdC(~ismember(string(allMiceIdC), markIds));
    assert(isempty(missingMark), ...
        'day4MarkC has no entry for: %s. Add them, or omit day4MarkC entirely.', ...
        strjoin(missingMark, ', '));

    cutoffPerMouse = NaT(nAllMice, 1);
    nBefore = zeros(nAllMice, 1);
    nAfter  = zeros(nAllMice, 1);

    keepRow = false(numel(mouseId), 1);
    for im = 1:nAllMice
        rowsIm = (mouseId == allMiceIdC{im});
        cutoffPerMouse(im) = markDates(find(markIds == string(allMiceIdC{im}), 1, 'first'));
        keepIm = rowsIm & ~isnat(sessDate) & (sessDate >= cutoffPerMouse(im));
        keepRow = keepRow | keepIm;
        nBefore(im) = sum(rowsIm);
        nAfter(im)  = sum(keepIm);
    end

    if opt.verbose
        fprintf('\nDay-4 filter applied (column "%s"):\n', dateVar);
        fprintf('  %-8s %-10s %8s %8s\n', 'animal', 'cutoff', 'before', 'after');
        for im = 1:nAllMice
            fprintf('  %-8s %-10s %8d %8d\n', allMiceIdC{im}, ...
                datestr(cutoffPerMouse(im), 'mm/dd/yy'), nBefore(im), nAfter(im));
        end
    end

    tooFew = allMiceIdC(nAfter < 2);
    assert(isempty(tooFew), ...
        'After day-4 filtering these animals have <2 sessions: %s. Analysis is not possible.', ...
        strjoin(tooFew, ', '));

    day4Filter.applied     = true;
    day4Filter.cutoffDates = cutoffPerMouse;
    day4Filter.nBefore     = nBefore;
    day4Filter.nAfter      = nAfter;
    day4Filter.dateVar     = dateVar;
end

mouseId    = mouseId(keepRow);
sessWithin = sessWithin(keepRow);
Yd         = Yd(keepRow, :);

%% -------------------- fixed right-aligned bookkeeping --------------------
perMouseMax   = zeros(nAllMice, 1);
perMouseCount = zeros(nAllMice, 1);
for im = 1:nAllMice
    sw = sessWithin(mouseId == allMiceIdC{im});
    perMouseMax(im)   = max(sw);
    perMouseCount(im) = numel(sw);
end
minDepth = min(perMouseCount);

if isempty(opt.nSessAnalyze)
    nSessAnalyze = minDepth;
else
    nSessAnalyze = opt.nSessAnalyze;
    if nSessAnalyze > minDepth
        error('Requested nSessAnalyze=%d exceeds the minimum available depth (%d).', nSessAnalyze, minDepth);
    end
end
assert(nSessAnalyze >= 2, ...
    'nSessAnalyze resolved to %d -- need at least 2 session positions. Check the day-4 filter.', nSessAnalyze);
if contrastStat == "slopediff"
    assert(nSessAnalyze >= 3, ...
        'contrastStat=''slopeDiff'' needs at least 3 positions to be meaningful (have %d).', nSessAnalyze);
end

sessFromEnd = nan(numel(mouseId), 1);
for im = 1:nAllMice
    rowsIm = (mouseId == allMiceIdC{im});
    sessFromEnd(rowsIm) = perMouseMax(im) - sessWithin(rowsIm);
end
sessPos = nSessAnalyze - sessFromEnd;

for im = 1:nAllMice
    posIm = sessPos(mouseId == allMiceIdC{im});
    missingPos = setdiff(1:nSessAnalyze, posIm(:)');
    assert(isempty(missingPos), ...
        'Animal %s has no session at right-aligned position(s) %s -- the common window is not complete (gap in sessWithin?).', ...
        allMiceIdC{im}, mat2str(missingPos));
end

if opt.verbose
    fprintf('\nRight-aligned window: nSessAnalyze = %d (min retained depth across %d mice).\n', ...
        nSessAnalyze, nAllMice);
end

%% -------------------- shared, non-circular Mahalanobis ruler --------------------
W_mahal = [];
if amongMetric == "mahal"
    SigmaGlobal = cov(Yd, 'omitrows');
    SigmaGlobal = (SigmaGlobal + SigmaGlobal')/2 + 1e-10*eye(size(SigmaGlobal,1));
    W_mahal = local_cholStable(SigmaGlobal);
end

%% -------------------- core reusable subroutines --------------------
    function amongTraj = local_amongTrajectory(memberIdC)
        amongTraj = nan(nSessAnalyze, 1);
        for pp = 1:nSessAnalyze
            rowsAtP = ismember(mouseId, string(memberIdC)) & (sessPos == pp);
            pts = Yd(rowsAtP, :);
            if size(pts, 1) < 2
                continue;   % need at least 2 points to have a pairwise distance
            end
            switch amongMetric
                case "euclidean"
                    dists = pdist(pts, 'euclidean');
                case "mahal"
                    ptsW = pts / W_mahal;
                    dists = pdist(ptsW, 'euclidean');
            end
            amongTraj(pp) = mean(dists);
        end
    end

    % Per-position signed difference -- the quantity every summary below
    % is built from, so the whole-window statistic and the per-position
    % tests are guaranteed to describe the same underlying contrast.
    function d = local_diffPerPos(traj1, traj2)
        switch direction
            case "slow_minus_fast", d = traj2 - traj1;
            case "fast_minus_slow", d = traj1 - traj2;
        end
    end

    function c = local_summarize(d)
        switch contrastStat
            case "meandiff"
                c = mean(d, 'omitnan');
            case "finaldiff"
                c = d(end);
            case "slopediff"
                % OLS slope of the gap against position; positive means the
                % groups separate further toward each animal's final session.
                x = (1:numel(d))';
                ok = isfinite(d);
                if sum(ok) < 3
                    c = NaN;
                else
                    b = polyfit(x(ok), d(ok), 1);
                    c = b(1);
                end
        end
    end

%% -------------------- observed --------------------
amongY_group1_obs = local_amongTrajectory(group1IdC);
amongY_group2_obs = local_amongTrajectory(group2IdC);
obsPerPos   = local_diffPerPos(amongY_group1_obs, amongY_group2_obs);
obsContrast = local_summarize(obsPerPos);

fprintf('\n============================================================\n');
fprintf('Among-distance permutation test (metric: %s, contrast: %s, direction: %s)\n', ...
    amongMetric, contrastStat, direction);
fprintf('Day-4 filter: %s | nSessAnalyze = %d (right-aligned common window)\n', ...
    local_onOff(day4Filter.applied), nSessAnalyze);
fprintf('Observed whole-window contrast: %.4f\n', obsContrast);

%% -------------------- null: exact enumeration or Monte Carlo --------------------
nPartitions = nchoosek(nAllMice, n1);
isExact = opt.useExact && (nPartitions <= opt.maxExact);

if isExact
    fprintf('Using EXACT enumeration: all C(%d,%d) = %d partitions.\n', nAllMice, n1, nPartitions);
    combos = nchoosek(1:nAllMice, n1);
    nDraws = size(combos, 1);
else
    fprintf('Using Monte Carlo: %d random partitions.\n', opt.nPerm);
    nDraws = opt.nPerm;
end

nullContrast = nan(nDraws, 1);
nullPerPos   = nan(nDraws, nSessAnalyze);   % the per-position null, free of extra cost
tPerm = tic;
for pI = 1:nDraws
    if isExact
        idx1 = combos(pI, :);
    else
        permOrder = randperm(nAllMice);
        idx1 = permOrder(1:n1);
    end
    idx2 = setdiff(1:nAllMice, idx1);

    traj1 = local_amongTrajectory(allMiceIdC(idx1));
    traj2 = local_amongTrajectory(allMiceIdC(idx2));

    dPerm = local_diffPerPos(traj1, traj2);
    nullPerPos(pI, :)  = dPerm(:)';
    nullContrast(pI)   = local_summarize(dPerm);

    if opt.verbose && (mod(pI, 500) == 0 || pI == nDraws)
        fprintf('  %d/%d (%.1fs elapsed)\n', pI, nDraws, toc(tPerm));
    end
end

pFloor = 1 / (nDraws + 1);
pValue = (1 + sum(nullContrast >= obsContrast, 'omitnan')) / (nDraws + 1);

%% -------------------- per-position inference --------------------
% Uncorrected: each position against its own null column.
pPerPos_uncorr = nan(nSessAnalyze, 1);
for pp = 1:nSessAnalyze
    pPerPos_uncorr(pp) = (1 + sum(nullPerPos(:,pp) >= obsPerPos(pp), 'omitnan')) / (nDraws + 1);
end

% Max-statistic FWER control: one null built from the LARGEST contrast
% across positions within each partition. Because the same relabeling
% generates the whole row, this automatically accounts for the
% correlation between positions -- no independence assumption, and no
% need to know how strongly neighbouring positions covary.
nullMax = max(nullPerPos, [], 2, 'omitnan');
pPerPos_maxStat = nan(nSessAnalyze, 1);
for pp = 1:nSessAnalyze
    pPerPos_maxStat(pp) = (1 + sum(nullMax >= obsPerPos(pp), 'omitnan')) / (nDraws + 1);
end

% BH-FDR across positions, as the less conservative alternative.
pPerPos_BH = local_bhFDR(pPerPos_uncorr);

perPosTable = table((1:nSessAnalyze)', amongY_group1_obs, amongY_group2_obs, obsPerPos, ...
    pPerPos_uncorr, pPerPos_maxStat, pPerPos_BH, ...
    'VariableNames', {'position', ['among_' groupNames{1}], ['among_' groupNames{2}], ...
                      'diff', 'p_uncorr', 'p_maxStat', 'p_BH'});

fprintf('Completed in %.1f sec.\n', toc(tPerm));
fprintf('Whole-window: observed = %.4f | null mean=%.4f SD=%.4f | p = %.4f', ...
    obsContrast, mean(nullContrast,'omitnan'), std(nullContrast,'omitnan'), pValue);
if pValue <= pFloor + eps
    fprintf('  <- at the floor (beat all %d partitions)', nDraws);
end
fprintf('\n');
fprintf('Exact p floor = 1/%d = %.4f -- no p below this is attainable.\n', nDraws + 1, pFloor);
fprintf('\nPer-position (position %d = each animal''s final session):\n', nSessAnalyze);
disp(perPosTable);

%% -------------------- optional plot --------------------
% Handles initialised here (not inside the if-blocks) so the save section
% below can test them whether or not either figure was drawn.
figTraj = [];
figDiff = [];

if opt.doPlot
    figTraj = figure('Color','w');
    hold on;
    gc = opt.groupColors;
    h1 = plot(1:nSessAnalyze, amongY_group1_obs, '-o', 'Color', gc(1,:), ...
        'MarkerFaceColor', gc(1,:), 'LineWidth', 2, 'DisplayName', groupNames{1});
    h2 = plot(1:nSessAnalyze, amongY_group2_obs, '-o', 'Color', gc(2,:), ...
        'MarkerFaceColor', gc(2,:), 'LineWidth', 2, 'DisplayName', groupNames{2});

    % Mark positions passing the FWER-corrected threshold, and (open
    % symbol) those passing only uncorrected -- so the figure never
    % implies more evidence than the correction supports.
    yTop = max([amongY_group1_obs; amongY_group2_obs], [], 'omitnan');
    yPad = 0.04 * range([amongY_group1_obs; amongY_group2_obs]);
    for pp = 1:nSessAnalyze
        if pPerPos_maxStat(pp) < opt.alphaPerPos
            plot(pp, yTop + yPad, 'k*', 'MarkerSize', 9, 'HandleVisibility', 'off');
        elseif pPerPos_uncorr(pp) < opt.alphaPerPos
            plot(pp, yTop + yPad, 'ko', 'MarkerSize', 6, 'HandleVisibility', 'off');
        end
    end

    xlabel('Session position (right-aligned; last = each animal''s final session)');
    ylabel(sprintf('Mean pairwise distance among group (%s)', amongMetric));
    title(sprintf(['Among-group dispersion (day-4 filter: %s)\n' ...
                   'whole-window %s p = %.4g   |   * FWER p<%.2f, o uncorrected only'], ...
        local_onOff(day4Filter.applied), char(contrastStat), pValue, opt.alphaPerPos));
    legend([h1 h2], 'Location','best');
    grid on; box off;
    xticks(1:nSessAnalyze);
    xlim([0.8, nSessAnalyze + 0.2]);
    if ~isempty(opt.yLimTraj)
        % An explicit range can hide points -- say so rather than silently
        % clipping, since a clipped trajectory still looks like a complete
        % curve. Applied AFTER the markers so yMark stays in view when the
        % range has headroom; if it doesn't, the markers are clipped too.
        yAllTraj = [amongY_group1_obs; amongY_group2_obs];
        nOut = sum(yAllTraj < opt.yLimTraj(1) | yAllTraj > opt.yLimTraj(2));
        if nOut > 0
            warning('runAmongDistancePermutationTest:yLimClipping', ...
                '%d trajectory point(s) fall outside yLimTraj [%g %g] and are clipped.', ...
                nOut, opt.yLimTraj(1), opt.yLimTraj(2));
        end
        ylim(opt.yLimTraj);
    end
    set(gca, 'TickDir', 'out');
    hold off;
end

%% -------------------- difference figure: observed vs. null band --------------------
%  The same information the test uses, drawn directly: the observed
%  per-position gap against the distribution of gaps the relabelings
%  produce. Three reference levels are drawn --
%    * the null median, i.e. where a typical relabeling lands;
%    * the per-position (1-alpha) null quantile, which is exactly what
%      p_uncorr < alpha means; and
%    * the max-statistic (1-alpha) quantile, drawn as a SINGLE horizontal
%      line because that null is taken over the max across positions and
%      so has no positional structure. Points above it are the
%      FWER-significant ones, which is why the '*' markers sit precisely
%      where the observed curve crosses it.
if opt.doPlotDiff
    nullLo    = prctile(nullPerPos, opt.bandPct(1), 1)';
    nullHi    = prctile(nullPerPos, opt.bandPct(2), 1)';
    nullMed   = median(nullPerPos, 1, 'omitnan')';
    uncorrThr = prctile(nullPerPos, 100*(1 - opt.alphaPerPos), 1)';
    fwerThr   = prctile(nullMax,    100*(1 - opt.alphaPerPos));

    xPos = (1:nSessAnalyze)';

    figDiff = figure('Color','w');
    hold on;

    fill([xPos; flipud(xPos)], [nullLo; flipud(nullHi)], [0.6 0.6 0.6], ...
        'FaceAlpha', 0.25, 'EdgeColor', 'none', ...
        'DisplayName', sprintf('null %g-%g%%', opt.bandPct(1), opt.bandPct(2)));
    plot(xPos, nullMed, '-', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.2, ...
        'DisplayName', 'null median');
    plot(xPos, uncorrThr, ':', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.4, ...
        'DisplayName', sprintf('null %g%% (per position)', 100*(1-opt.alphaPerPos)));
    yline(fwerThr, '--', 'Color', [0.20 0.20 0.20], 'LineWidth', 1.4, ...
        'DisplayName', sprintf('max-stat %g%% (FWER)', 100*(1-opt.alphaPerPos)));
    yline(0, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8, 'HandleVisibility', 'off');

    % Neutral dark for the observed DIFFERENCE: it is not a group, so it
    % must not wear either group's colour (fast is now orange, which the
    % old red would have been mistaken for).
    obsColor = [0.15 0.15 0.15];
    plot(xPos, obsPerPos, '-o', 'Color', obsColor, 'LineWidth', 2.2, ...
        'MarkerFaceColor', obsColor, 'MarkerSize', 7, ...
        'DisplayName', sprintf('observed (%s - %s)', groupNames{2}, groupNames{1}));

    % same significance markers as the trajectory figure
    yAll  = [obsPerPos; nullLo; nullHi];
    yMark = max(yAll, [], 'omitnan') + 0.05*range(yAll);
    for pp = 1:nSessAnalyze
        if pPerPos_maxStat(pp) < opt.alphaPerPos
            plot(pp, yMark, 'k*', 'MarkerSize', 9, 'HandleVisibility', 'off');
        elseif pPerPos_uncorr(pp) < opt.alphaPerPos
            plot(pp, yMark, 'ko', 'MarkerSize', 6, 'HandleVisibility', 'off');
        end
    end

    xlabel('Session position (right-aligned; last = each animal''s final session)');
    ylabel(sprintf('Among-distance difference, %s - %s (%s)', ...
        groupNames{2}, groupNames{1}, amongMetric));
    title(sprintf(['Observed group difference vs. permutation null (day-4 filter: %s)\n' ...
                   'whole-window %s p = %.4g   |   * FWER p<%.2f, o uncorrected only'], ...
        local_onOff(day4Filter.applied), char(contrastStat), pValue, opt.alphaPerPos));
    legend('Location', 'best');
    grid on; box off;
    xticks(1:nSessAnalyze);
    xlim([0.8, nSessAnalyze + 0.2]);
    set(gca, 'TickDir', 'out');
    hold off;
end

%% -------------------- optional figure save --------------------
if opt.printFig
    if strlength(strtrim(string(opt.figSaveDir))) == 0
        error('printFig is true but no ''figSaveDir'' was provided.');
    end
    outDirFig = char(string(opt.figSaveDir));
    if exist(outDirFig, 'dir') ~= 7
        mkdir(outDirFig);
    end

    % Same tag scheme as the .mat save, so a run's figures and its stats
    % file are matched by filename.
    figTag = sprintf('%s_%s_%s', local_tagOrDefault(opt.saveTag), ...
        local_onOffTag(day4Filter.applied), char(contrastStat));
    dStr = char(datetime('today','Format','MMddyy'));

    figsToSave = {figTraj, 'trajectory'; figDiff, 'diffVsNull'};
    for k = 1:size(figsToSave, 1)
        fH = figsToSave{k, 1};
        if isempty(fH) || ~isgraphics(fH, 'figure'), continue; end
        outFile = fullfile(outDirFig, sprintf('amongDistance_%s_%s_%s.pdf', ...
            figsToSave{k, 2}, figTag, dStr));
        set(fH, 'InvertHardcopy', 'off');
        print(fH, outFile, '-dpdf', '-painters', '-bestfit');
        fprintf('Saved figure:\n  %s\n', outFile);
    end
end

%% -------------------- package output --------------------
amongOut = struct();
amongOut.amongX        = (1:nSessAnalyze)';
amongOut.amongY_group1 = amongY_group1_obs;
amongOut.amongY_group2 = amongY_group2_obs;
amongOut.groupNames    = groupNames;

amongOut.obsContrast   = obsContrast;
amongOut.nullContrast  = nullContrast;
amongOut.pValue        = pValue;

amongOut.obsPerPos       = obsPerPos;
amongOut.nullPerPos      = nullPerPos;
amongOut.pPerPos_uncorr  = pPerPos_uncorr;
amongOut.pPerPos_maxStat = pPerPos_maxStat;
amongOut.pPerPos_BH      = pPerPos_BH;
amongOut.perPosTable     = perPosTable;

amongOut.pFloor        = pFloor;
amongOut.nPartitions   = nPartitions;
amongOut.isExact       = isExact;
amongOut.direction     = char(direction);
amongOut.contrastStat  = char(contrastStat);
amongOut.amongMetric   = char(amongMetric);
amongOut.nSessAnalyze  = nSessAnalyze;
amongOut.day4Filter    = day4Filter;
amongOut.figs          = struct('trajectory', figTraj, 'difference', figDiff);

%% -------------------- optional save --------------------
if opt.doSave
    if strlength(strtrim(string(opt.saveDir))) == 0
        error('doSave is true but no ''saveDir'' was provided.');
    end
    outDir = char(string(opt.saveDir));
    if exist(outDir, 'dir') ~= 7
        mkdir(outDir);
    end
    tag = strtrim(char(string(opt.saveTag)));
    if isempty(tag)
        tag = 'unlabeled';
        warning('runAmongDistancePermutationTest:noSaveTag', 'No ''saveTag'' provided -- using "unlabeled".');
    end
    % Filter state and contrast both go in the filename so variants from
    % the same stream and date cannot overwrite each other.
    tag = sprintf('%s_%s_%s', tag, local_onOffTag(day4Filter.applied), char(contrastStat));
    dateStr  = char(datetime('today','Format','MMddyy'));
    saveName = sprintf('amongDistancePermTest_%s_%s.mat', tag, dateStr);
    saveFullPath = fullfile(outDir, saveName);
    save(saveFullPath, 'amongOut', 'groupDefs');
    amongOut.save = struct('didSave', true, 'file', saveFullPath);
    fprintf('\nSaved to:\n%s\n', saveFullPath);
else
    amongOut.save = struct('didSave', false, 'file', '');
end

end

%% ========================= local helpers =========================
function s = local_onOff(tf)
if tf, s = 'ON'; else, s = 'off'; end
end


function t = local_tagOrDefault(s)
t = strtrim(char(string(s)));
if isempty(t), t = 'unlabeled'; end
end


function s = local_onOffTag(tf)
if tf, s = 'day4on'; else, s = 'allSess'; end
end


function q = local_bhFDR(pvals)
% Benjamini-Hochberg step-up FDR. Returns monotonic q-values capped at 1.
pvals = pvals(:);
n = numel(pvals);
[pSorted, sortIdx] = sort(pvals);
qSorted = pSorted .* n ./ (1:n)';
qSorted = flipud(cummin(flipud(qSorted)));
qSorted = min(qSorted, 1);
q = nan(n, 1);
q(sortIdx) = qSorted;
end


function varName = local_findDateVar(sessInfo, requested)
if strlength(strtrim(string(requested))) > 0
    varName = char(string(requested));
    assert(ismember(varName, sessInfo.Properties.VariableNames), ...
        'Requested sessionDateVar "%s" is not a column of sessInfo. Available: %s', ...
        varName, strjoin(sessInfo.Properties.VariableNames, ', '));
    return;
end

candidates = {'header', 'sessionHeader', 'sessHeader', 'sessionId', 'sessName', 'sessionName'};
hit = candidates(ismember(candidates, sessInfo.Properties.VariableNames));
assert(~isempty(hit), ...
    ['day4MarkC was supplied but no session-date column was found in sessInfo ' ...
     '(looked for: %s). Pass ''sessionDateVar'' explicitly. Available columns: %s'], ...
    strjoin(candidates, ', '), strjoin(sessInfo.Properties.VariableNames, ', '));
varName = hit{1};
end


function dt = local_parseHeaderDate(headerCol)
% Parses headers like 'm1613_050725' or 'm1613_050725-1' into datetimes.
% MMddyy convention; the optional '-N' same-day-rerun suffix is ignored
% for date comparison. Returns NaT for anything that doesn't match.
headerCol = string(headerCol(:));
n = numel(headerCol);
dt = NaT(n, 1);
for i = 1:n
    tok = regexp(char(headerCol(i)), '_(\d{6})(?:-\d+)?$', 'tokens', 'once');
    if isempty(tok), continue; end
    try
        dt(i) = datetime(tok{1}, 'InputFormat', 'MMddyy');
    catch
        % leave as NaT
    end
end
end


function dt = local_parseMMDDYY(dateCol)
dateCol = string(dateCol(:));
n = numel(dateCol);
dt = NaT(n, 1);
for i = 1:n
    s = strtrim(char(dateCol(i)));
    if numel(s) ~= 6, continue; end
    try
        dt(i) = datetime(s, 'InputFormat', 'MMddyy');
    catch
        % leave as NaT
    end
end
end


function R = local_cholStable(Sigma)
Sigma = (Sigma + Sigma')/2;
reg = 0;
for it = 1:10
    try
        R = chol(Sigma + reg*eye(size(Sigma,1)));
        return;
    catch
        if reg == 0
            reg = 1e-12;
        else
            reg = reg * 10;
        end
    end
end
error('Cholesky failed even after regularization. Sigma may be badly conditioned.');
end