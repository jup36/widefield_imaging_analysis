function out = plotMotifHAcrossSessions(HstatC, headerC, timeC, mId, motifId, trialType, varargin)
%PLOTMOTIFHACROSSSESSIONS
% Plot mean H (raw/z-scored motif activity) across sessions for one mouse
% and one motif -- with an early(light blue) -> late(dark blue) session
% color gradient, matching the style of the GLM-TDR projection-score
% trajectory panels elsewhere in this pipeline.
%
% This is the intuitive, single-motif companion to the GLM-TDR projection
% analyses: it plots the motif's raw temporal activation (from
% descriptiveH_fourWay) directly, not a GLM-fitted or axis-projected
% quantity.
%
% out = plotMotifHAcrossSessions(HstatC, headerC, timeC, mId, motifId, trialType, 'Name', value, ...)
%
% REQUIRED INPUTS
%   HstatC    : [J x S] cell array (from batchComputeMotifH_fourWayStats).
%               HstatC{j,s} is itself a 1 x K cell array (one struct per
%               motif), each with fields .mean.(hit|miss|cr|fa) and
%               .sem.(hit|miss|cr|fa), each [1 x T].
%   headerC   : [J x S] cell array of session headers (e.g. 'm1045_122424'),
%               same size as HstatC.
%   timeC     : [J x S] cell array of per-session time vectors (e.g. Hs.winCtrs
%               from the batch script), same size as HstatC.
%   mId       : char/string, mouse ID, e.g. "m1045"
%   motifId   : scalar, motif index (matches HstatC{j,s}{motifId}.motifIdx)
%   trialType : char/string, one of "hit","miss","cr","fa" (case-insensitive;
%               also accepts common synonyms: "correctrejection" -> cr,
%               "falsealarm" -> fa).
%               MAY ALSO be a cell array of multiple trial types, e.g.
%               {'hit','cr'} -- one subplot panel is created per trial
%               type, arranged side by side in a single figure, and the
%               figure width is scaled by the number of panels.
%
% NAME-VALUE
%   'tBounds'          : [] (default) or [tMin tMax] to restrict x-axis/plotted range
%   'smoothingFactor'  : 0 (default) or positive integer, passed to smooth2a
%   'FadeToWhite'      : 0.85 (default) -- earliest session/subgroup's blend toward white
%   'LineWidth'        : 2 (default)
%   'ShowSEM'          : false (default) -- if true, adds a light shaded
%                        SEM band around each line (same color, low alpha)
%   'SEMFaceAlpha'     : 0.10 (default), used only if ShowSEM=true
%   'eventLines'       : [0 2 4] (default) -- x-locations for dashed vertical
%                        reference lines (e.g. tone onset / response window
%                        bounds / tone offset). Pass [] to omit.
%   'subGroupSession'  : [] (default, numeric scalar otherwise). If empty,
%                        every session is plotted as its own line (original
%                        behavior), colored on a continuous early->late
%                        gradient, with a small gradient legend inset.
%                        If a positive integer G is given, the mouse's
%                        chronologically-sorted sessions (the ones with
%                        usable data for a given trial type) are split into
%                        G non-overlapping, contiguous, chronological
%                        subgroups (G=2 -> early/late; G=3 ->
%                        early/intermediate/late; G>3 -> early/mid.../late),
%                        sized as evenly as possible with any remainder
%                        distributed to the EARLIEST groups first (e.g.
%                        N=8,G=3 -> [3,3,2]; N=7,G=3 -> [3,2,2]). Within
%                        each subgroup, sessions' mean traces are averaged
%                        (mean +/- SEM across sessions in that subgroup,
%                        via meanstdsem) into a single representative line;
%                        a discrete legend (with each subgroup's label and
%                        session count) is shown instead of the continuous
%                        gradient inset.
%   'matchYLim'        : true (default). Only relevant when multiple trial
%                        types are given (subplot mode): if true, all panels'
%                        y-axis limits are set to the same shared range (the
%                        min/max across all panels' auto-scaled ranges) after
%                        plotting, so magnitudes are directly comparable
%                        across subplots. Event lines/labels are drawn AFTER
%                        this shared range is applied, so the dashed
%                        tone/response markers span the full shared range in
%                        every panel. Ignored (no-op) when only one trial
%                        type is given.
%   'MakeFigure'       : true (default)
%   'figSaveLogic'     : false (default)
%   'figSaveDir'       : '' (default; required if figSaveLogic=true)
%   'figName'          : "" (default) -- if empty, auto-built from
%                        mId/motifId/trialType(s)
%
% OUTPUT
%   out : struct.
%     If trialType was a single (non-cell) value, out has the SAME flat
%     fields as before (.mId, .motifId, .trialType, .headersUsed,
%     .sessDtUsed, .t, .trajC, .semC, .nTrialC, .figHandle), for backward
%     compatibility, PLUS out.byType.(canonicalTrialType) mirroring the
%     same content.
%
%     If trialType was a cell array of multiple types, the flat top-level
%     fields are omitted; use out.byType.(canonicalTrialType) for each
%     type's results, plus out.figHandle (shared figure) and out.mId /
%     out.motifId.
%
%     When 'subGroupSession' is used, each out.byType.(type) entry also
%     has fields .groupLabels, .groupSizes, .groupTrajC, .groupSemC
%     (one row per subgroup) instead of/alongside the per-session .trajC.
%
% EXAMPLE
%   % Single trial type, per-session gradient (original behavior)
%   plotMotifHAcrossSessions(HstatC, headerC, timeC, "m1045", 6, "cr", ...
%       'tBounds', [-1 5], 'smoothingFactor', 3);
%
%   % Two trial types side by side
%   plotMotifHAcrossSessions(HstatC, headerC, timeC, "m1045", 6, {'hit','cr'}, ...
%       'tBounds', [-1 5], 'smoothingFactor', 3);
%
%   % Early vs late subgroups (2 groups) for a single trial type
%   plotMotifHAcrossSessions(HstatC, headerC, timeC, "m1045", 6, "cr", ...
%       'subGroupSession', 2);
%
%   % Terciles, two trial types
%   plotMotifHAcrossSessions(HstatC, headerC, timeC, "m1045", 6, {'hit','cr'}, ...
%       'subGroupSession', 3);
%
% See also: descriptiveH_fourWay, batchComputeMotifH_fourWayStats

%% -------------------- parse --------------------
p = inputParser;
p.FunctionName = mfilename;

p.addRequired('HstatC', @(x) iscell(x));
p.addRequired('headerC', @(x) iscell(x) && isequal(size(x), size(HstatC)));
p.addRequired('timeC', @(x) iscell(x) && isequal(size(x), size(HstatC)));
p.addRequired('mId', @(s) ischar(s) || isstring(s));
p.addRequired('motifId', @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addRequired('trialType', @(s) ischar(s) || isstring(s) || iscell(s));

p.addParameter('tBounds', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2 && x(1)<x(2)));
p.addParameter('smoothingFactor', 0, @(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('FadeToWhite', 0.85, @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
p.addParameter('LineWidth', 2, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('ShowSEM', false, @(x) islogical(x) && isscalar(x));
p.addParameter('SEMFaceAlpha', 0.10, @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
p.addParameter('eventLines', [0 2 4], @(x) isnumeric(x));
p.addParameter('subGroupSession', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>=1 && x==round(x)));
p.addParameter('matchYLim', true, @(x) islogical(x) && isscalar(x));
p.addParameter('MakeFigure', true, @(x) islogical(x) && isscalar(x));
p.addParameter('figSaveLogic', false, @(x) islogical(x) && isscalar(x));
p.addParameter('figSaveDir', '', @(s) ischar(s) || isstring(s));
p.addParameter('figName', "", @(s) ischar(s) || isstring(s));

p.parse(HstatC, headerC, timeC, mId, motifId, trialType, varargin{:});
opt = p.Results;

mId = string(mId);

% -------------------- normalize trialType to a cellstr list --------------------
wasCell = iscell(trialType);
if wasCell
    ttList = trialType;
else
    ttList = {trialType};
end
nTypes = numel(ttList);
ttCanonList = cell(1, nTypes);
for i = 1:nTypes
    ttCanonList{i} = local_resolve_trial_type(ttList{i});
end

%% -------------------- locate sessions for this mouse (shared across trial types) --------------------
hdrFlat = headerC(:);
isHdr = ~cellfun(@isempty, hdrFlat);

linIdxAll = find(isHdr);
hdrFlatS  = string(hdrFlat(isHdr));

matchMask = contains(hdrFlatS, mId);
matchLinIdx = linIdxAll(matchMask);
matchHdrS   = hdrFlatS(matchMask);

assert(~isempty(matchLinIdx), 'No sessions found for mId=%s in headerC.', mId);

[dtKey, ok] = local_parse_header_datetime(matchHdrS);
if ok
    [~, ord] = sort(dtKey, 'ascend');
else
    ord = (1:numel(matchHdrS))';
    warning('plotMotifHAcrossSessions:DateParseFail', ...
        'Could not parse dates from one or more headers for mId=%s; using original order.', mId);
end

matchLinIdx = matchLinIdx(ord);
matchHdrS   = matchHdrS(ord);
if ok
    dtKey = dtKey(ord);
else
    dtKey = NaT(numel(matchHdrS),1);
end

[J, S] = size(HstatC);

%% -------------------- figure / subplot layout --------------------
if opt.MakeFigure
    if nTypes > 1
        panelW = 4.5; panelH = 4.0;
        hFig = figure('Color', 'w', 'Units', 'inches', ...
            'Position', [1, 1, panelW*nTypes, panelH]);
    else
        hFig = figure('Color', 'w');
    end
else
    hFig = [];
end

%% -------------------- per-trial-type processing + plotting --------------------
byType = struct();
axList = gobjects(1, nTypes);

for ti = 1:nTypes
    ttCanon = ttCanonList{ti};

    if opt.MakeFigure
        if nTypes > 1
            ax = subplot(1, nTypes, ti);
        else
            ax = gca;
        end
        hold(ax, 'on');
        axList(ti) = ax;
    else
        ax = [];
    end

    res = local_plot_data_one_type( ...
        ax, HstatC, timeC, matchLinIdx, matchHdrS, dtKey, J, S, ...
        mId, motifId, ttCanon, opt);

    byType.(ttCanon) = res;
end

% -------------------- match y-limits across subplots (before drawing
% event lines / legends, so dashed markers span the shared range) --------------------
if opt.MakeFigure && nTypes > 1 && opt.matchYLim
    ylAll = nan(nTypes, 2);
    for ti = 1:nTypes
        ylAll(ti,:) = ylim(axList(ti));
    end
    sharedYLim = [min(ylAll(:,1)), max(ylAll(:,2))];
    for ti = 1:nTypes
        ylim(axList(ti), sharedYLim);
    end
end

% -------------------- finalize each axes (event lines, labels, legend/gradient) --------------------
for ti = 1:nTypes
    ttCanon = ttCanonList{ti};
    if opt.MakeFigure
        local_finalize_one_axes(axList(ti), byType.(ttCanon), opt, mId, motifId, ttCanon);
    end
end

if opt.MakeFigure && nTypes > 1
    sgtitle(sprintf('%s | motif %d', mId, motifId), 'Interpreter', 'none');
end

%% -------------------- save --------------------
if opt.MakeFigure && opt.figSaveLogic
    figSaveDir = char(string(opt.figSaveDir));
    assert(~isempty(figSaveDir), 'figSaveDir must be provided when figSaveLogic=true.');
    if ~exist(figSaveDir, 'dir')
        mkdir(figSaveDir);
    end

    if strlength(string(opt.figName)) > 0
        fbase = char(string(opt.figName));
    else
        ttTag = strjoin(ttCanonList, '_');
        fbase = sprintf('%s_motif%d_%s_acrossSessionH', char(mId), motifId, ttTag);
    end
    fbase = local_sanitize_filename([fbase '.pdf']);
    fpath = fullfile(figSaveDir, fbase);

    set(hFig, 'PaperPositionMode', 'auto');
    print(hFig, fpath, '-dpdf', '-painters', '-bestfit');
    fprintf('[plotMotifHAcrossSessions] saved: %s\n', fpath);
end

%% -------------------- pack output --------------------
out = struct();
out.mId = char(mId);
out.motifId = motifId;
out.byType = byType;
out.figHandle = hFig;

if ~wasCell
    % Backward-compatible flat fields, mirroring the single result
    flat = byType.(ttCanonList{1});
    fns = fieldnames(flat);
    for i = 1:numel(fns)
        out.(fns{i}) = flat.(fns{i});
    end
end

end

%% ========================================================================
% Core: process + plot DATA ONLY for one trial type into one axes
% (event lines / labels / legend / gradient inset are added later, in
%  local_finalize_one_axes, AFTER y-limits are matched across subplots)
% ========================================================================
function res = local_plot_data_one_type( ...
    ax, HstatC, timeC, matchLinIdx, matchHdrS, dtKey, J, S, mId, motifId, trialTypeCanon, opt)

Ns = numel(matchLinIdx);

trajC   = cell(1, Ns);
semC    = cell(1, Ns);
nTrialC = nan(1, Ns);
headersUsed = cell(1, Ns);
tPlotFinal = [];
keepSess = false(1, Ns);

for k = 1:Ns
    linIdx = matchLinIdx(k);
    [j, s] = ind2sub([J, S], linIdx);

    hstatSess = HstatC{j, s};
    tvecSess  = timeC{j, s};

    if isempty(hstatSess) || isempty(tvecSess)
        continue;
    end

    if motifId > numel(hstatSess)
        warning('plotMotifHAcrossSessions:MotifOutOfRange', ...
            'Session %s: motifId=%d exceeds available motifs (%d). Skipping.', ...
            matchHdrS(k), motifId, numel(hstatSess));
        continue;
    end

    mstat = hstatSess{motifId};
    if isempty(mstat) || ~isfield(mstat,'mean') || ~isfield(mstat.mean, trialTypeCanon)
        warning('plotMotifHAcrossSessions:MissingField', ...
            'Session %s: motif %d missing field mean.%s. Skipping.', ...
            matchHdrS(k), motifId, trialTypeCanon);
        continue;
    end

    tvec = tvecSess(:)';
    xMean = mstat.mean.(trialTypeCanon);
    xMean = xMean(:)';
    xSem  = mstat.sem.(trialTypeCanon);
    xSem  = xSem(:)';

    if isfield(mstat,'nTrial') && isfield(mstat.nTrial, trialTypeCanon)
        nTrialC(k) = mstat.nTrial.(trialTypeCanon);
    end

    assert(numel(tvec) == numel(xMean), ...
        'Session %s: time vector length (%d) does not match mean trace length (%d).', ...
        matchHdrS(k), numel(tvec), numel(xMean));

    tMask = true(size(tvec));
    if ~isempty(opt.tBounds)
        tMask = (tvec >= opt.tBounds(1)) & (tvec <= opt.tBounds(2));
    end
    assert(any(tMask), 'tBounds excluded all timestamps for session %s.', matchHdrS(k));

    tPlot = tvec(tMask);
    xMeanPlot = xMean(tMask);
    xSemPlot  = xSem(tMask);

    if opt.smoothingFactor > 0
        xMeanPlot = smooth2a(double(xMeanPlot), 0, opt.smoothingFactor);
        xSemPlot  = smooth2a(double(xSemPlot),  0, opt.smoothingFactor);
    else
        xMeanPlot = double(xMeanPlot);
        xSemPlot  = double(xSemPlot);
    end

    trajC{k} = xMeanPlot;
    semC{k}  = xSemPlot;
    headersUsed{k} = char(matchHdrS(k));
    keepSess(k) = true;

    if isempty(tPlotFinal)
        tPlotFinal = tPlot;
    end
end

assert(any(keepSess), 'No sessions had usable data for mId=%s, motif=%d, trialType=%s.', ...
    mId, motifId, trialTypeCanon);

trajC       = trajC(keepSess);
semC        = semC(keepSess);
nTrialC     = nTrialC(keepSess);
headersUsed = headersUsed(keepSess);
dtKeyUsed   = dtKey(keepSess);
NsUsed      = sum(keepSess);

baseCol = [0 0 1]; % blue

res = struct();
res.mId = char(mId);
res.motifId = motifId;
res.trialType = trialTypeCanon;
res.headersUsed = headersUsed(:)';
res.sessDtUsed = dtKeyUsed(:);
res.t = tPlotFinal(:)';
res.trajC = trajC;
res.semC = semC;
res.nTrialC = nTrialC;

%% -------------------- either per-session OR subgrouped plotting --------------------
res.plotMode = 'perSession';
res.baseCol = baseCol;

if isempty(opt.subGroupSession)
    % ---- original per-session behavior ----
    wVec = linspace(opt.FadeToWhite, 0, NsUsed);

    if ~isempty(ax)
        for k = 1:NsUsed
            col = local_blend_to_white(baseCol, wVec(k));

            if opt.ShowSEM
                yU = trajC{k} + semC{k};
                yL = trajC{k} - semC{k};
                fill(ax, [tPlotFinal, fliplr(tPlotFinal)], [yU, fliplr(yL)], col, ...
                    'FaceAlpha', opt.SEMFaceAlpha, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            end

            plot(ax, tPlotFinal, trajC{k}, 'Color', col, 'LineWidth', opt.LineWidth);
        end
    end

else
    % ---- subgrouped (early/late or terciles etc.) behavior ----
    G = opt.subGroupSession;
    assert(G <= NsUsed, ...
        'subGroupSession=%d exceeds the number of usable sessions (%d) for trialType=%s.', ...
        G, NsUsed, trialTypeCanon);

    [groupIdxC, groupLabels] = local_partition_sessions(NsUsed, G);
    groupSizes = cellfun(@numel, groupIdxC);

    groupTrajC = cell(G, 1);
    groupSemC  = cell(G, 1);

    wVec = linspace(opt.FadeToWhite, 0, G);
    hLine = gobjects(1, G);

    for g = 1:G
        idx = groupIdxC{g};
        Mstack = cell2mat(trajC(idx)'); % [nSessInGroup x T]
        [mG, ~, seG] = meanstdsem(Mstack);
        groupTrajC{g} = mG;
        groupSemC{g}  = seG;

        if ~isempty(ax)
            col = local_blend_to_white(baseCol, wVec(g));

            if opt.ShowSEM
                yU = mG + seG;
                yL = mG - seG;
                fill(ax, [tPlotFinal, fliplr(tPlotFinal)], [yU, fliplr(yL)], col, ...
                    'FaceAlpha', opt.SEMFaceAlpha, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            end

            hLine(g) = plot(ax, tPlotFinal, mG, 'Color', col, 'LineWidth', opt.LineWidth);
        end
    end

    res.groupLabels = groupLabels;
    res.groupSizes  = groupSizes;
    res.groupTrajC  = groupTrajC;
    res.groupSemC   = groupSemC;
    res.plotMode    = 'subGroup';

    if ~isempty(ax)
        res.hLine = hLine; % handles for legend, built later in local_finalize_one_axes
    end
end

end

%% ========================================================================
% Local helpers
% ========================================================================

function local_finalize_one_axes(ax, res, opt, mId, motifId, trialTypeCanon)
% Draws event lines (using the CURRENT ylim -- already matched across
% subplots by this point if matchYLim was requested), axis labels/title,
% and the appropriate legend (discrete, subgroup mode) or gradient inset
% (continuous, per-session mode).

if ~isempty(opt.eventLines)
    yl = ylim(ax);
    for ev = opt.eventLines(:)'
        plot(ax, [ev ev], yl, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');
    end
    ylim(ax, yl);
end

xlabel(ax, 'Time (s)');
ylabel(ax, sprintf('Motif %d activity (%s, a.u.)', motifId, upper(trialTypeCanon)), 'Interpreter', 'none');
title(ax, sprintf('%s trials', upper(trialTypeCanon)), 'Interpreter', 'none');
box(ax, 'off');

switch res.plotMode
    case 'perSession'
        local_add_session_gradient_legend(ax, res.baseCol, opt.FadeToWhite);
    case 'subGroup'
        legStr = arrayfun(@(g) sprintf('%s (n=%d)', res.groupLabels{g}, res.groupSizes(g)), ...
            1:numel(res.groupLabels), 'UniformOutput', false);
        legend(ax, res.hLine, legStr, 'Location', 'best', 'Interpreter', 'none', 'Box', 'off');
end

end

function canon = local_resolve_trial_type(trialType)
tt = lower(strtrim(char(string(trialType))));
tt = regexprep(tt, '[\s_-]+', '');

switch tt
    case 'hit'
        canon = 'hit';
    case 'miss'
        canon = 'miss';
    case {'cr','correctrejection'}
        canon = 'cr';
    case {'fa','falsealarm'}
        canon = 'fa';
    otherwise
        error(['Unrecognized trialType="%s". Supported values (case-insensitive): ' ...
               '"hit", "miss", "cr"/"correctrejection", "fa"/"falsealarm".'], char(string(trialType)));
end
end

function [dtKey, ok] = local_parse_header_datetime(hdrS)
hdrS = string(hdrS(:));
n = numel(hdrS);
dtKey = NaT(n,1);
ok = true;

for i = 1:n
    h = char(hdrS(i));
    tok = regexp(h, '_(\d{6})(?:-(\d+))?$', 'tokens', 'once');
    if isempty(tok)
        ok = false;
        return;
    end
    mmddyy = tok{1};
    suf = 0;
    if numel(tok) >= 2 && ~isempty(tok{2})
        suf = str2double(tok{2});
        if ~isfinite(suf), suf = 0; end
    end
    try
        d0 = datetime(mmddyy, 'InputFormat', 'MMddyy');
    catch
        ok = false;
        return;
    end
    dtKey(i) = d0 + seconds(suf);
end
end

function colOut = local_blend_to_white(colIn, w)
colIn = colIn(:)';
if numel(colIn) ~= 3
    colIn = [0 0 1];
end
w = max(min(w,1),0);
colOut = (1-w)*colIn + w*[1 1 1];
colOut = max(min(colOut,1),0);
end

function [groupIdxC, groupLabels] = local_partition_sessions(Ns, G)
% Partition 1:Ns (chronologically ordered, earliest first) into G
% contiguous, non-overlapping groups. Sizes are as even as possible; any
% remainder is distributed to the EARLIEST groups first.
%   e.g. Ns=8, G=3 -> sizes [3,3,2]
%        Ns=7, G=3 -> sizes [3,2,2]

assert(G >= 1 && G <= Ns, 'subGroupSession (%d) must be between 1 and the number of usable sessions (%d).', G, Ns);

baseSize = floor(Ns / G);
remainder = mod(Ns, G);

sizes = repmat(baseSize, 1, G);
sizes(1:remainder) = sizes(1:remainder) + 1;

edges = [0, cumsum(sizes)];
groupIdxC = cell(1, G);
for g = 1:G
    groupIdxC{g} = (edges(g)+1):edges(g+1);
end

groupLabels = local_make_group_labels(G);
end

function labels = local_make_group_labels(G)
if G == 1
    labels = {'all'};
elseif G == 2
    labels = {'early', 'late'};
elseif G == 3
    labels = {'early', 'intermediate', 'late'};
else
    labels = cell(1, G);
    labels{1} = 'early';
    labels{G} = 'late';
    for g = 2:G-1
        labels{g} = sprintf('mid%d', g-1);
    end
end
end

function local_add_session_gradient_legend(ax, baseCol, fadeToWhite)
% Small inset gradient bar + "early"/"late" labels + "sessions" title,
% matching the style of existing GLM-TDR projection-score panels.

fig = ancestor(ax, 'figure');
axPos = get(ax, 'Position'); % normalized axes-in-figure units

insetW = axPos(3) * 0.32;
insetH = axPos(4) * 0.05;
insetX = axPos(1) + axPos(3) * 0.55;
insetY = axPos(2) + axPos(4) * 0.90;

axInset = axes('Parent', fig, 'Position', [insetX, insetY, insetW, insetH]);

nGrad = 256;
wVec = linspace(fadeToWhite, 0, nGrad);
gradImg = zeros(1, nGrad, 3);
for i = 1:nGrad
    gradImg(1,i,:) = local_blend_to_white(baseCol, wVec(i));
end

image(axInset, gradImg);
set(axInset, 'XTick', [], 'YTick', [], 'Box', 'off', 'Visible', 'off');

text(axInset, 0, 1.8, 'sessions', 'Units', 'normalized', ...
    'HorizontalAlignment', 'center', 'FontAngle', 'italic', 'FontSize', 9);
text(axInset, 0.02, -0.6, 'early', 'Units', 'normalized', ...
    'HorizontalAlignment', 'left', 'FontAngle', 'italic', 'FontSize', 8);
text(axInset, 0.98, -0.6, 'late', 'Units', 'normalized', ...
    'HorizontalAlignment', 'right', 'FontAngle', 'italic', 'FontSize', 8);

axes(ax); %#ok<LAXES> % restore focus to the main plotting axes
end

function fn = local_sanitize_filename(fn)
fn = char(fn);
bad = '<>:"/\|?*';
for k = 1:numel(bad)
    fn(fn==bad(k)) = '_';
end
fn = regexprep(fn, '\s+', '_');
end
