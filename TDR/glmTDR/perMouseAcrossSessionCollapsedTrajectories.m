function out = perMouseAcrossSessionCollapsedTrajectories(rezProjStats, timestamps, varargin)
%PERMOUSEACROSSSESSIONCOLLAPSEDTRAJECTORIES
% Plot 1D projected score trajectories across sessions for one mouse,
% using rezProjStats.perMouse{row,1}.collapsed.<field>, where each field is
% [nSessions x nTime x nAxes].
%
% out = perMouseAcrossSessionCollapsedTrajectories(rezProjStats, timestamps, 'Name', value, ...)
%
% REQUIRED INPUTS
%   rezProjStats : struct with field .perMouse
%                  Expected usage:
%                  rezProjStats.perMouse{rowIdx,1}.collapsed.muGo   [S x T x D]
%                  rezProjStats.perMouse{rowIdx,1}.collapsed.muNoGo [S x T x D]
%                  rezProjStats.perMouse{rowIdx,1}.collapsed.diff   [S x T x D]
%                  rezProjStats.perMouse{rowIdx,1}.collapsed.energy [S x T x D]
%
%   timestamps   : [1 x T] or [T x 1] time vector aligned to 2nd dim
%                  of collapsed arrays
%
% NAME-VALUE (core)
%   'mouseId'         : "m1045" (required)
%   'targetName'      : "NoGoToneOn_1" (required; must exist in axes names)
%   'trialField'      : "go" | "nogo" | "diff" | "energy"
%                       case-insensitive
%                       default = "go"
%
%   'FlipDiffForNoGoAxes'
%                     : true (default)
%                       If true and trialField="diff" and targetName contains
%                       "nogo" case-insensitively, flip sign of diff.
%
%                       Stored diff convention:
%                           diff = Go - NoGo
%
%                       Plotted convention after flip for NoGo axes:
%                           -diff = NoGo - Go
%
%                       This is useful for visualizing positive emergence of
%                       NoGo-axis selectivity across learning.
%
% NAME-VALUE (time + filtering)
%   'tBounds'         : [] (default) or [tMin tMax]
%   'smoothingFactor' : 0 (default) or positive integer
%   'dateLaterThan'   : [] (default) or "MMDDYY" (include sessions >= this date)
%   'dateEarlierThan' : [] (default) or "MMDDYY" (include sessions <= this date)
%   'day4MarkC'       : [] (default) or Nx2 cell/string array {mouseId, "MMDDYY"; ...}
%                      If provided, overrides dateLaterThan by using the mouse-specific Day4 date.
%
% NAME-VALUE (visual)
%   'FadeToWhite'     : 0.85 (default)
%   'LineWidth'       : 2 (default)
%   'MakeFigure'      : true (default)
%   'lineColor'       : [0 0 1] (default)
%   'ylim'            : [] (default). If provided as [yMin yMax], sets the
%                       y-axis limits of the plot accordingly. Ignored if
%                       MakeFigure=false.
%
% NAME-VALUE (save)
%   'figSaveDir'      : '' (default). If non-empty, saves PDF to this directory.
%
% OUTPUT
%   out : struct with fields
%       .mouseId
%       .targetName
%       .targetCol
%       .collapsedField
%       .plotFieldLabel
%       .doFlipDiff
%       .diffOriginalConvention
%       .diffPlotConvention
%       .headersUsed
%       .sessDtUsed
%       .t
%       .trajC
%       .figHandle
%       .rowIdx
%

% -------------------- parse --------------------
p = inputParser;
p.FunctionName = mfilename;

p.addRequired('rezProjStats', @(s)isstruct(s) && isfield(s,'perMouse'));
p.addRequired('timestamps', @(t) isnumeric(t) && isvector(t));

p.addParameter('mouseId', "", @(s)ischar(s)||isstring(s));
p.addParameter('targetName', "", @(s)ischar(s)||isstring(s));
p.addParameter('trialField', "go", @(s)ischar(s)||isstring(s));

p.addParameter('FlipDiffForNoGoAxes', true, @(x) islogical(x) && isscalar(x));

p.addParameter('tBounds', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2 && x(1)<x(2)));
p.addParameter('smoothingFactor', 0, @(x) isnumeric(x) && isscalar(x) && x>=0);

p.addParameter('dateLaterThan', [], @(s) isempty(s) || ischar(s) || isstring(s));
p.addParameter('dateEarlierThan', [], @(s) isempty(s) || ischar(s) || isstring(s));
p.addParameter('day4MarkC', [], @(c) isempty(c) || iscell(c) || isstring(c));

p.addParameter('FadeToWhite', 0.85, @(x) isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('LineWidth', 2, @(x) isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('MakeFigure', true, @(x) islogical(x)&&isscalar(x));
p.addParameter('lineColor', [0 0 1], @(x) isnumeric(x) && numel(x)==3 && all(x>=0 & x<=1));
p.addParameter('ylim', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2 && x(1)<x(2)));

p.addParameter('figSaveDir', '', @(s) isempty(s) || ischar(s) || isstring(s));

p.parse(rezProjStats, timestamps, varargin{:});
opt = p.Results;

mouseId    = string(opt.mouseId);
targetName = string(opt.targetName);

assert(strlength(mouseId)>0, 'mouseId is required, e.g. "m1045".');
assert(strlength(targetName)>0, 'targetName is required, e.g. "NoGoToneOn_1".');

% -------------------- timestamps sanity --------------------
tvec = timestamps(:)'; % 1 x T

% -------------------- resolve effective dateLaterThan --------------------
effectiveDateLaterThan = opt.dateLaterThan;

if ~isempty(opt.day4MarkC)
    d4 = opt.day4MarkC;
    if isstring(d4)
        d4 = cellstr(d4);
    end

    if iscell(d4)
        assert(size(d4,2) == 2, ...
            'day4MarkC must be an Nx2 cell/string array: {mouseId, "MMDDYY"; ...}.');

        mouseCol = string(d4(:,1));
        dateCol  = string(d4(:,2));

        hit = find(mouseCol == mouseId, 1, 'first');
        if isempty(hit)
            error('day4MarkC provided, but mouseId=%s was not found in day4MarkC.', mouseId);
        end

        effectiveDateLaterThan = dateCol(hit);
    else
        error('day4MarkC must be empty or an Nx2 cell/string array.');
    end
end

% -------------------- locate row for this mouse USING rezProjStats.perMouse --------------------
assert(iscell(rezProjStats.perMouse), 'rezProjStats.perMouse must be a cell array.');

nRows = size(rezProjStats.perMouse, 1);
mouseIdPerRow = strings(nRows,1);

for r = 1:nRows
    if ~isempty(rezProjStats.perMouse{r,1}) && isstruct(rezProjStats.perMouse{r,1}) ...
            && isfield(rezProjStats.perMouse{r,1}, 'mouseId') ...
            && ~isempty(rezProjStats.perMouse{r,1}.mouseId)
        mouseIdPerRow(r) = string(rezProjStats.perMouse{r,1}.mouseId);
    end
end

rowIdx = find(strcmpi(mouseIdPerRow, mouseId), 1, 'first');
assert(~isempty(rowIdx), ...
    'Could not find mouseId=%s in rezProjStats.perMouse{:,1}.mouseId.', mouseId);

% -------------------- access per-mouse stats --------------------
pm = rezProjStats.perMouse{rowIdx,1};

assert(~isempty(pm) && isstruct(pm), ...
    'rezProjStats.perMouse{%d,1} is empty or invalid.', rowIdx);

assert(isfield(pm,'mouseId') && strcmpi(string(pm.mouseId), mouseId), ...
    'Resolved rowIdx=%d does not match requested mouseId=%s.', rowIdx, mouseId);

assert(isfield(pm,'collapsed') && isstruct(pm.collapsed), ...
    'rezProjStats.perMouse{%d,1} must contain a struct field .collapsed.', rowIdx);

assert(isfield(pm,'sessions') && isstruct(pm.sessions) && isfield(pm.sessions,'headers'), ...
    'rezProjStats.perMouse{%d,1}.sessions.headers is required.', rowIdx);

% -------------------- map trialField -> collapsed field --------------------
collapsedField = local_resolve_collapsed_field(opt.trialField);

assert(isfield(pm.collapsed, collapsedField), ...
    'collapsed field "%s" was not found in rezProjStats.perMouse{%d,1}.collapsed.', ...
    collapsedField, rowIdx);

Xall = pm.collapsed.(collapsedField);   % [S x T x D]
assert(ndims(Xall) == 3, ...
    'collapsed.%s must be [nSessions x nTime x nAxes].', collapsedField);

% -------------------- resolve axis names --------------------
nameList = local_resolve_axis_names(pm);

targetCol = glmNameToColumns(nameList, cellstr(targetName));
assert(~isempty(targetCol), ...
    'targetName=%s not found in available axis names.', targetName);
targetCol = targetCol(1);

% -------------------- optional sign flip for NoGo-axis diff --------------------
% Stored convention:
%   collapsed.diff = muGo - muNoGo
%
% For NoGo axes, the sign flip gives:
%   -collapsed.diff = muNoGo - muGo
%
% This makes positive values correspond to stronger NoGo-vs-Go expression.
isDiffField = strcmpi(collapsedField, 'diff');
isNoGoAxis  = contains(lower(char(targetName)), 'nogo');

doFlipDiff = opt.FlipDiffForNoGoAxes && isDiffField && isNoGoAxis;

if doFlipDiff
    plotFieldLabel = 'diff flipped: NoGo - Go';
    diffPlotConvention = 'NoGo - Go';
    fprintf('[%s] Sign-flipping diff for NoGo axis "%s": plotting NoGo - Go.\n', ...
        mfilename, char(targetName));
else
    plotFieldLabel = collapsedField;
    if isDiffField
        diffPlotConvention = 'Go - NoGo';
    else
        diffPlotConvention = '';
    end
end

diffOriginalConvention = 'Go - NoGo';

% -------------------- gather session headers from pm.sessions.headers --------------------
headersRaw = pm.sessions.headers;

if isstring(headersRaw)
    headersUsedAll = cellstr(headersRaw(:)');
elseif iscell(headersRaw)
    headersUsedAll = cell(size(headersRaw));
    for i = 1:numel(headersRaw)
        headersUsedAll{i} = char(string(headersRaw{i}));
    end
    headersUsedAll = reshape(headersUsedAll, 1, []);
else
    error('rezProjStats.perMouse{%d,1}.sessions.headers must be a cell array or string array.', rowIdx);
end

nSessFromHeaders = numel(headersUsedAll);
nSessFromData = size(Xall,1);

assert(nSessFromHeaders == nSessFromData, ...
    ['Mouse %s row %d: sessions.headers has %d sessions, ', ...
     'but collapsed data has %d sessions. These must match exactly.'], ...
     mouseId, rowIdx, nSessFromHeaders, nSessFromData);

nSess = nSessFromData;

% -------------------- session date parsing + filtering --------------------
sessDt = NaT(1, nSess);

for jj = 1:nSess
    h = headersUsedAll{jj};
    if isempty(h)
        continue;
    end

    [dtSess, ok] = parse_header_mmddyy_ignore_suffix(h);
    if ok
        sessDt(jj) = dtSess;
    end
end

if ~isempty(effectiveDateLaterThan)
    dt0 = datetime(char(string(effectiveDateLaterThan)), 'InputFormat','MMddyy');
    keepLater = (sessDt >= dt0);
else
    keepLater = true(size(sessDt));
end

if ~isempty(opt.dateEarlierThan)
    dt1 = datetime(char(string(opt.dateEarlierThan)), 'InputFormat','MMddyy');
    keepEarlier = (sessDt <= dt1);
else
    keepEarlier = true(size(sessDt));
end

validDataMask = false(1, nSess);
for jj = 1:nSess
    xj = squeeze(Xall(jj,:,targetCol));
    validDataMask(jj) = ~isempty(xj) && any(~isnan(xj(:)));
end

keepMask = validDataMask & keepLater & keepEarlier;
idxKeep = find(keepMask);

assert(~isempty(idxKeep), ...
    'No sessions survived selection for mouseId=%s.', mouseId);

dtKeep = sessDt(idxKeep);

if all(~isnat(dtKeep))
    [~, ord] = sort(dtKeep, 'ascend');
    idxKeep = idxKeep(ord);
    dtKeep = dtKeep(ord);
else
    dtKeep = sessDt(idxKeep);
end

% -------------------- time bounds --------------------
T = size(Xall,2);

assert(numel(tvec)==T, ...
    'timestamps length (%d) must match data time dimension T (%d).', ...
    numel(tvec), T);

tMask = true(1, T);

if ~isempty(opt.tBounds)
    tMask = (tvec >= opt.tBounds(1)) & (tvec <= opt.tBounds(2));
end

tIdx = find(tMask);
assert(~isempty(tIdx), 'tBounds excluded all timestamps.');

tPlot = tvec(tIdx);

% -------------------- plot --------------------
if opt.MakeFigure
    hFig = figure('Color','w');
    hold on;
else
    hFig = [];
end

Ns = numel(idxKeep);
trajC = cell(1, Ns);
headersUsed = cell(1, Ns);

wVec = linspace(opt.FadeToWhite, 0, Ns);
baseCol = opt.lineColor;

hLeg = gobjects(1, Ns);
legC = cell(1, Ns);

for k = 1:Ns
    jj = idxKeep(k);

    x = squeeze(Xall(jj, :, targetCol));
    x = x(:);
    x = x(tIdx);

    % Stored diff = Go - NoGo.
    % For NoGo axes, optionally plot NoGo - Go.
    if doFlipDiff
        x = -x;
    end

    if opt.smoothingFactor > 0
        x = smooth2a(double(x), opt.smoothingFactor, 0);
    else
        x = double(x);
    end

    trajC{k} = x;
    headersUsed{k} = headersUsedAll{jj};

    col = blend_to_white(baseCol, wVec(k));

    if opt.MakeFigure
        hLeg(k) = plot(tPlot, x, ...
            'LineWidth', opt.LineWidth, ...
            'Color', col);
    end

    legC{k} = headersUsedAll{jj};
end

if opt.MakeFigure
    xlabel('time', 'Interpreter', 'none');

    ylabel(sprintf('%s (%s)', plotFieldLabel, targetName), ...
        'Interpreter', 'none');

    title(sprintf('%s | %s | %s', mouseId, plotFieldLabel, char(targetName)), ...
        'Interpreter','none');

    legend(hLeg(isgraphics(hLeg)), legC(isgraphics(hLeg)), ...
        'Location','eastoutside', ...
        'Interpreter','none');

    box off;

    if ~isempty(opt.ylim)
        ylim(opt.ylim);
    end
end

% -------------------- save figure --------------------
figSaveDir = char(string(opt.figSaveDir));

if opt.MakeFigure && ~isempty(figSaveDir)
    if ~exist(figSaveDir, 'dir')
        mkdir(figSaveDir);
    end

    todayStr = datestr(now, 'mmddyy');

    if doFlipDiff
        fieldForFile = 'diff_NoGoMinusGo';
    else
        fieldForFile = collapsedField;
    end

    fbase = sprintf('%s_%s_%s_acrossSession_%s.pdf', ...
        char(mouseId), fieldForFile, char(targetName), todayStr);

    fbase = sanitize_filename(fbase);
    fpath = fullfile(figSaveDir, fbase);

    set(gcf, 'PaperPositionMode','auto');
    print(gcf, fpath, '-dpdf', '-painters', '-bestfit');

    fprintf('[perMouseAcrossSessionCollapsedTrajectories] saved: %s\n', fpath);
end

% -------------------- pack output --------------------
out = struct();

out.mouseId        = char(mouseId);
out.targetName     = char(targetName);
out.targetCol      = targetCol;

out.collapsedField = collapsedField;
out.plotFieldLabel = plotFieldLabel;

out.doFlipDiff     = doFlipDiff;
out.diffOriginalConvention = diffOriginalConvention;
out.diffPlotConvention     = diffPlotConvention;

out.headersUsed    = headersUsed;
out.sessDtUsed     = dtKeep(:);
out.t              = tPlot(:);
out.trajC          = trajC;
out.figHandle      = hFig;
out.rowIdx         = rowIdx;

end

%% ===================== HELPERS =====================

function collapsedField = local_resolve_collapsed_field(trialField)

tf = lower(strtrim(char(string(trialField))));
tf = regexprep(tf, '\s+', '');

switch tf
    case 'go'
        collapsedField = 'muGo';

    case {'nogo','no-go','cr','correctrejection'}
        collapsedField = 'muNoGo';

    case 'diff'
        collapsedField = 'diff';

    case 'energy'
        collapsedField = 'energy';

    otherwise
        error(['Unrecognized trialField="%s". Supported values are: ', ...
               '"go", "nogo", "diff", "energy".'], char(string(trialField)));
end

end

function nameList = local_resolve_axis_names(pm)

nameList = [];

if isfield(pm, 'axes') && isstruct(pm.axes) && ...
        isfield(pm.axes, 'names') && ~isempty(pm.axes.names)
    nameList = pm.axes.names;
end

assert(~isempty(nameList), ...
    'Could not resolve axis names from pm.axes.names.');

nameList = cellstr(string(nameList(:)'));

end

function cols = glmNameToColumns(glmNameListC, targetNameC)

glmNameS = string(glmNameListC(:));
targS    = string(targetNameC(:));

cols = [];

for i = 1:numel(targS)
    hit = find(glmNameS == targS(i), 1, 'first');
    if ~isempty(hit)
        cols(end+1) = hit; %#ok<AGROW>
    end
end

end

function [dtSess, ok] = parse_header_mmddyy_ignore_suffix(header)

h = char(string(header));
tok = regexp(h, '_(\d{6})(?:-\d+)?$', 'tokens', 'once');

if isempty(tok)
    dtSess = NaT;
    ok = false;
    return;
end

mmddyy = tok{1};

try
    dtSess = datetime(mmddyy, 'InputFormat','MMddyy');
    ok = true;
catch
    dtSess = NaT;
    ok = false;
end

end

function colOut = blend_to_white(colIn, w)

colIn = colIn(:)';

if numel(colIn) ~= 3
    colIn = [0 0 0];
end

w = max(min(w,1),0);

colOut = (1-w)*colIn + w*[1 1 1];
colOut = max(min(colOut,1),0);

end

function fn = sanitize_filename(fn)

fn = char(fn);
bad = '<>:"/\|?*';

for k = 1:numel(bad)
    fn(fn==bad(k)) = '_';
end

fn = regexprep(fn, '\s+', '_');

end