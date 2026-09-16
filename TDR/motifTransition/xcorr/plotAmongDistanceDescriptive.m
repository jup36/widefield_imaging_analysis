function amongDescOut = plotAmongDistanceDescriptive(Y, sessInfo, groupDefs, varargin)
%PLOTAMONGDISTANCEDESCRIPTIVE
%   PURELY DESCRIPTIVE visualization (NOT a statistical test) of mean
%   pairwise distance among each group's own members, at EVERY session
%   position where at least 2 group members have a valid session --
%   shown for BOTH left-alignment (raw session order) and right-
%   alignment (aligned to each group's own final session), side by side
%   in one figure, with all groups overlaid on each panel.
%
%   Unlike runAmongDistancePermutationTest / _blocked (which restrict to
%   a common window so every animal contributes a complete row for a
%   valid permutation-test design), this function makes NO such
%   restriction -- it uses every session position where the data allow a
%   pairwise distance to be computed at all (>=2 members present), and
%   leaves NaN (a visible gap in the plotted line) everywhere else. No
%   p-value, no null distribution -- this is for looking at the full
%   shape of the data, nothing more.
%
%   DAY-4 FILTERING ('day4MarkC', optional)
%   ---------------------------------------
%   Restricts the plot to sessions on or after each animal's own "day 4"
%   date -- the date the full Go/NoGo task with a 1:1 trial ratio began.
%   Pass an {mouseId, 'MMDDYY'} cell array, the same one used for the
%   permutation test, so the descriptive figure and the formal test
%   describe the same session set.
%
%   NOTE ON LEFT ALIGNMENT UNDER FILTERING: each animal's day-4 date
%   falls at a different ORIGINAL sessWithin index, so plotting the raw
%   index would start each group at a different x even though both are at
%   the same task stage. 'reindexLeftAfterFilter' (default TRUE) therefore
%   renumbers each animal's retained sessions from 1 for the LEFT panel,
%   giving every animal a common origin. Set it false to see the true
%   original session numbers instead, at the cost of staggered starts.
%   The x-label states which convention is in force.
%
%   RIGHT ALIGNMENT is computed from a single GLOBAL shift across every
%   animal in groupDefs, so all groups share one right edge and each
%   animal's final session lands on the same x. (An earlier version
%   shifted per group, aligning each group to its own longest-recorded
%   animal, which placed the groups on different x ranges and made the
%   overlay misleading.)
%
%   amongDescOut = plotAmongDistanceDescriptive(Y, sessInfo, groupDefs, ...)
%
% INPUTS
%   Y         : [S x dimFull] MDS coordinates for ONE stream.
%   sessInfo  : matching sessInfo table (needs mouseId, sessWithin; also a
%               header-like column if day4MarkC is used).
%   groupDefs : scalar struct, ANY number of groups (not restricted to
%               exactly 2, since this is descriptive, not a permutation
%               test), e.g.
%                 groupDefs.fast = {'m1044','m1045','m1092','m1094'};
%                 groupDefs.slow = {'m1048','m1049','m1613','m1859','m1873'};
%
% NAME-VALUE ARGS
%   'day4MarkC'    : {} (default, no filtering) or an {N x 2} cell array
%                    of {mouseId, 'MMDDYY'}. Every animal appearing in
%                    groupDefs must have an entry.
%   'sessionDateVar' : '' (default = auto-detect) or the name of the
%                    sessInfo column holding the session header/date.
%   'reindexLeftAfterFilter' : true (default). When the day-4 filter is
%                    on, renumbers each animal's retained sessions 1..n
%                    for the LEFT panel so all groups share an origin.
%                    Set false to keep the true original session index.
%                    Has no effect when the filter is off.
%   'dims'         : which columns of Y to use. Default: 1:size(Y,2).
%   'amongMetric'  : 'euclidean' (default) or 'mahal'. If 'mahal', the
%                    whitening covariance is computed ONCE from the
%                    RETAINED pooled dataset (all groups, post-filter) --
%                    a shared ruler, not built from any one group.
%   'groupColors'  : [nGroups x 3] RGB, one row per group in
%                    fieldnames(groupDefs) order. Default: MATLAB's
%                    standard color order.
%   'lineWidth'    : default 2.0
%   'markerSize'   : default 6
%   'gridOn'       : default true
%   'figWidthScale': default 1.8 (two panels side by side need more width)
%   'figSaveDir'   : directory to save a PDF. Default: '' (no save).
%   'printFig'     : true/false, whether to save. Default: false.
%   'reprint'      : overwrite existing file. Default: false.
%   'figureNameBase' : filename prefix. Default: 'amongDistanceDescriptive_'.
%
% OUTPUT (amongDescOut)
%   .left, .right        : struct arrays, one entry per group, each with
%                          .groupName, .x, .y
%   .fig, .ax            : figure and [1x2] axes handles
%   .amongMetric, .day4Filter, .save

%% -------------------- parse options --------------------
p = inputParser;
p.addParameter('day4MarkC', {}, @(x) isempty(x) || (iscell(x) && size(x,2) == 2));
p.addParameter('sessionDateVar', '', @(s) ischar(s) || isstring(s));
p.addParameter('reindexLeftAfterFilter', true, @(x) islogical(x) && isscalar(x));
p.addParameter('dims', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
p.addParameter('amongMetric', 'euclidean', @(s) any(strcmpi(string(s), ["euclidean","mahal"])));
p.addParameter('groupColors', [], @(x) isempty(x) || (isnumeric(x) && size(x,2)==3));
p.addParameter('lineWidth', 2.0, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('markerSize', 6, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('gridOn', true, @(x) islogical(x) && isscalar(x));
p.addParameter('figWidthScale', 1.8, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('figSaveDir', '', @(s) ischar(s) || isstring(s));
p.addParameter('printFig', false, @(x) islogical(x) && isscalar(x));
p.addParameter('reprint', false, @(x) islogical(x) && isscalar(x));
p.addParameter('figureNameBase', 'amongDistanceDescriptive_', @(s) ischar(s) || isstring(s));
p.parse(varargin{:});
opt = p.Results;

amongMetric = lower(string(opt.amongMetric));

%% -------------------- validate + setup --------------------
groupNames = fieldnames(groupDefs);
nGroups = numel(groupNames);
    
assert(all(ismember({'mouseId','sessWithin'}, sessInfo.Properties.VariableNames)), ...
    'sessInfo must have mouseId and sessWithin columns.');

mouseId    = string(sessInfo.mouseId);
sessWithin = sessInfo.sessWithin;

if isempty(opt.dims)
    dims = 1:size(Y,2);
else
    dims = opt.dims(:)';
end
Yd = Y(:, dims);

if isempty(opt.groupColors)
    groupColors = flipud(lines(nGroups));
else
    assert(size(opt.groupColors,1) == nGroups, 'groupColors must have one row per group.');
    groupColors = opt.groupColors;
end

allMembers = {};
for g = 1:nGroups
    allMembers = [allMembers, groupDefs.(groupNames{g})(:)']; %#ok<AGROW>
end
allMembers = unique(allMembers, 'stable');

%% -------------------- day-4 filtering (optional) --------------------
% Applied before anything else, so both alignments, the Mahalanobis ruler,
% and every plotted point refer to retained sessions only.
day4Filter = struct('applied', false, 'cutoffDates', [], 'nBefore', [], 'nAfter', [], ...
    'mice', {allMembers}, 'reindexedLeft', false);
keepRow = true(numel(mouseId), 1);

if ~isempty(opt.day4MarkC)
    dateVar  = local_findDateVar(sessInfo, opt.sessionDateVar);
    sessDate = local_parseHeaderDate(sessInfo.(dateVar));

    nUnparsed = sum(isnat(sessDate));
    if nUnparsed > 0
        warning('plotAmongDistanceDescriptive:unparsedDates', ...
            '%d session(s) in column "%s" could not be date-parsed and will be DROPPED by the day-4 filter.', ...
            nUnparsed, dateVar);
    end

    markIds   = string(opt.day4MarkC(:,1));
    markDates = local_parseMMDDYY(opt.day4MarkC(:,2));
    assert(~any(isnat(markDates)), 'One or more day4MarkC dates could not be parsed as MMDDYY.');

    missingMark = allMembers(~ismember(string(allMembers), markIds));
    assert(isempty(missingMark), ...
        'day4MarkC has no entry for: %s. Add them, or omit day4MarkC entirely.', ...
        strjoin(missingMark, ', '));

    nMice = numel(allMembers);
    cutoffPerMouse = NaT(nMice, 1);
    nBefore = zeros(nMice, 1);
    nAfter  = zeros(nMice, 1);

    keepRow = false(numel(mouseId), 1);
    for im = 1:nMice
        rowsIm = (mouseId == allMembers{im});
        cutoffPerMouse(im) = markDates(find(markIds == string(allMembers{im}), 1, 'first'));
        keepIm = rowsIm & ~isnat(sessDate) & (sessDate >= cutoffPerMouse(im));
        keepRow = keepRow | keepIm;
        nBefore(im) = sum(rowsIm);
        nAfter(im)  = sum(keepIm);
    end

    fprintf('\nDay-4 filter applied (column "%s"):\n', dateVar);
    fprintf('  %-8s %-10s %8s %8s\n', 'animal', 'cutoff', 'before', 'after');
    for im = 1:nMice
        fprintf('  %-8s %-10s %8d %8d\n', allMembers{im}, ...
            datestr(cutoffPerMouse(im), 'mm/dd/yy'), nBefore(im), nAfter(im));
    end

    tooFew = allMembers(nAfter < 1);
    if ~isempty(tooFew)
        warning('plotAmongDistanceDescriptive:noSessionsAfterFilter', ...
            'These animals have NO sessions after the day-4 cutoff and drop out entirely: %s', ...
            strjoin(tooFew, ', '));
    end

    day4Filter.applied     = true;
    day4Filter.cutoffDates = cutoffPerMouse;
    day4Filter.nBefore     = nBefore;
    day4Filter.nAfter      = nAfter;
    day4Filter.dateVar     = dateVar;
end

mouseId    = mouseId(keepRow);
sessWithin = sessWithin(keepRow);
Yd         = Yd(keepRow, :);

% Optional renumbering for the LEFT panel only. Off by default because it
% replaces the true session index with a within-filter index; on, it gives
% every animal a common origin so the left panel is interpretable again.
sessLeft = sessWithin;
if day4Filter.applied && opt.reindexLeftAfterFilter
    for im = 1:numel(allMembers)
        rowsIm = (mouseId == allMembers{im});
        if ~any(rowsIm), continue; end
        [~, ord] = sort(sessWithin(rowsIm), 'ascend');
        newIdx = nan(sum(rowsIm), 1);
        newIdx(ord) = 1:sum(rowsIm);
        sessLeft(rowsIm) = newIdx;
    end
    day4Filter.reindexedLeft = true;
end

%% -------------------- shared, non-circular Mahalanobis ruler (if requested) --------------------
% Computed ONCE from the RETAINED pooled dataset (every group's members,
% every retained session) -- a shared metric, not built from any single
% group. NOTE: with the day-4 filter on, this ruler differs from the
% unfiltered one, so absolute distance values are not comparable between
% a filtered and an unfiltered run -- only the shapes are.
W_mahal = [];
if amongMetric == "mahal"
    poolRows = ismember(mouseId, string(allMembers));
    SigmaGlobal = cov(Yd(poolRows, :), 'omitrows');
    SigmaGlobal = (SigmaGlobal + SigmaGlobal')/2 + 1e-10*eye(numel(dims));
    W_mahal = local_cholStable(SigmaGlobal);
end

%% -------------------- core: among-distance trajectory over a given x-vector --------------------
    function yVec = local_amongOverX(memberIdC, xForMember, xGrid)
        % xForMember: same length as mouseId, giving each row's ALIGNED x
        %             position (only meaningful for rows belonging to
        %             memberIdC; other rows are ignored below).
        yVec = nan(numel(xGrid), 1);
        rowsMember = ismember(mouseId, string(memberIdC));
        for xi = 1:numel(xGrid)
            rowsAtX = rowsMember & (xForMember == xGrid(xi));
            pts = Yd(rowsAtX, :);
            if size(pts,1) < 2
                continue;   % descriptive gap -- fewer than 2 members present at this position
            end
            switch amongMetric
                case "euclidean"
                    dists = pdist(pts, 'euclidean');
                case "mahal"
                    ptsW = pts / W_mahal;
                    dists = pdist(ptsW, 'euclidean');
            end
            yVec(xi) = mean(dists);
        end
    end

%% -------------------- global right-alignment shift --------------------
% Computed ONCE across every animal in groupDefs, not per group, so all
% groups share a common right edge: each animal's final session lands on
% the same x. Aligning each group to its own longest animal would put the
% groups on different x ranges, which makes an overlay unreadable.
uAll = unique(mouseId, 'stable');
uAll = uAll(ismember(uAll, string(allMembers)));
perMouseMaxAll = zeros(numel(uAll), 1);
for im = 1:numel(uAll)
    perMouseMaxAll(im) = max(sessWithin(mouseId == uAll(im)));
end
maxFinalAll = max(perMouseMaxAll);

xRightAll = sessWithin;
for im = 1:numel(uAll)
    rowsIm = (mouseId == uAll(im));
    xRightAll(rowsIm) = sessWithin(rowsIm) + (maxFinalAll - perMouseMaxAll(im));
end

%% -------------------- compute LEFT- and RIGHT-aligned trajectories, per group --------------------
leftResults  = struct('groupName', {}, 'x', {}, 'y', {});
rightResults = struct('groupName', {}, 'x', {}, 'y', {});

for g = 1:nGroups
    members = groupDefs.(groupNames{g});
    memberRows = ismember(mouseId, string(members));

    if ~any(memberRows)
        warning('plotAmongDistanceDescriptive:emptyGroup', ...
            'Group "%s" has no retained sessions -- skipping.', groupNames{g});
        leftResults(g)  = struct('groupName', groupNames{g}, 'x', [], 'y', []);
        rightResults(g) = struct('groupName', groupNames{g}, 'x', [], 'y', []);
        continue;
    end

    % ---- LEFT alignment: raw (or reindexed) session number, no shift ----
    xGridLeft = unique(sessLeft(memberRows), 'sorted');
    yLeft = local_amongOverX(members, sessLeft, xGridLeft);

    leftResults(g).groupName = groupNames{g};
    leftResults(g).x = xGridLeft;
    leftResults(g).y = yLeft;

    % ---- RIGHT alignment: uses the GLOBAL shift computed once above, so
    % every group's final session lands on the same rightmost x. (A
    % per-group shift would align each group to its own longest animal,
    % putting the two lines on different x ranges and defeating the
    % point of overlaying them.) ----
    xGridRight = unique(xRightAll(memberRows), 'sorted');
    yRight = local_amongOverX(members, xRightAll, xGridRight);

    rightResults(g).groupName = groupNames{g};
    rightResults(g).x = xGridRight;
    rightResults(g).y = yRight;
end

%% -------------------- plot --------------------
fig = figure('Color','w');
pos = fig.Position;
pos(3) = pos(3) * opt.figWidthScale;
fig.Position = pos;

tl = tiledlayout(fig, 1, 2, 'TileSpacing','compact', 'Padding','compact');

% Left-panel x-label depends on whether filtering shifted the origin.
if day4Filter.applied && day4Filter.reindexedLeft
    leftXLabel = 'Session (left-aligned; renumbered from each animal''s day-4 session)';
elseif day4Filter.applied
    leftXLabel = 'Session (left-aligned; ORIGINAL index -- animals start at different x)';
else
    leftXLabel = 'Session (left-aligned; raw session order)';
end

axLeft = nexttile(tl, 1);
hold(axLeft, 'on');
for g = 1:nGroups
    if isempty(leftResults(g).x), continue; end
    plot(axLeft, leftResults(g).x, leftResults(g).y, '-o', ...
        'LineWidth', opt.lineWidth, 'MarkerSize', opt.markerSize, ...
        'Color', groupColors(g,:), 'MarkerFaceColor', groupColors(g,:), ...
        'DisplayName', groupNames{g});
end
xlabel(axLeft, leftXLabel);
ylabel(axLeft, sprintf('Mean pairwise distance among group (%s)', amongMetric));
title(axLeft, 'Left-aligned');
if opt.gridOn, grid(axLeft, 'on'); end
box(axLeft, 'off');
legend(axLeft, 'Location', 'best');
hold(axLeft, 'off');

axRight = nexttile(tl, 2);
hold(axRight, 'on');
for g = 1:nGroups
    if isempty(rightResults(g).x), continue; end
    plot(axRight, rightResults(g).x, rightResults(g).y, '-o', ...
        'LineWidth', opt.lineWidth, 'MarkerSize', opt.markerSize, ...
        'Color', groupColors(g,:), 'MarkerFaceColor', groupColors(g,:), ...
        'DisplayName', groupNames{g});
end
xlabel(axRight, 'Session (right-aligned; last = each group''s own final session)');
ylabel(axRight, sprintf('Mean pairwise distance among group (%s)', amongMetric));
title(axRight, 'Right-aligned');
if opt.gridOn, grid(axRight, 'on'); end
box(axRight, 'off');
legend(axRight, 'Location', 'best');
hold(axRight, 'off');

if day4Filter.applied
    sgtitle(tl, 'Among-group pairwise dispersion (descriptive; no statistical test) -- DAY-4 SESSIONS ONLY');
else
    sgtitle(tl, 'Among-group pairwise dispersion (descriptive; no statistical test)');
end

%% -------------------- optional save --------------------
saveInfo = struct('didSave', false, 'file', '');
if opt.printFig
    if strlength(strtrim(string(opt.figSaveDir))) == 0
        error('printFig is true but no ''figSaveDir'' was provided.');
    end
    outDir = char(string(opt.figSaveDir));
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end
    dstr = datestr(now, 'mmddyy');
    % Filter state in the filename so a filtered and an unfiltered figure
    % from the same day cannot overwrite each other.
    fnBase = sprintf('%s%s_%s', char(string(opt.figureNameBase)), ...
        local_onOffTag(day4Filter.applied), dstr);
    outFile = fullfile(outDir, [fnBase '.pdf']);
    if ~(exist(outFile, 'file') && ~opt.reprint)
        set(fig, 'InvertHardcopy', 'off');
        print(fig, outFile, '-dpdf', '-painters', '-bestfit');
        saveInfo.didSave = true;
        saveInfo.file = outFile;
    end
end

%% -------------------- package output --------------------
amongDescOut = struct();
amongDescOut.left        = leftResults;
amongDescOut.right       = rightResults;
amongDescOut.fig         = fig;
amongDescOut.ax          = [axLeft, axRight];
amongDescOut.amongMetric = char(amongMetric);
amongDescOut.day4Filter  = day4Filter;
amongDescOut.save        = saveInfo;

end

%% ========================= local helpers =========================
function s = local_onOffTag(tf)
if tf, s = 'day4on'; else, s = 'allSess'; end
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