function anchor = collectPerMouseAndGlobalGlmPrjAxes(headerC, glmRezC, trIdC, varargin)
%collectPerMouseAndGlobalGlmPrjAxes  Build per-mouse + global GLM-TDR axes, plus per-session axes.
%
% PRINCIPLED POOLING (Purpose A):
%   For pooled anchors (global / allSessions pooled), we POOL DATA (Yhat_g) then SVD ONCE per group.
%   We do NOT average/merge axes across sessions.
%
% ASSUMPTIONS (as you specified):
%   - All sessions share the same predictor group names/order and motif order (K) is identical.
%   - Group definitions (cols) are consistent across sessions.
%
% REQUIRED INPUTS
%   headerC : cell (J x S) session headers
%   glmRezC : cell (J x S) glmRez structs (same size)
%   trIdC   : cell (J x S) trial-index structs (same size), fields:
%             .goI   : [N x 1] logical or index vector of Go trials
%             .nogoI : [N x 1] logical or index vector of NoGo trials
%
% NAME-VALUE
%   'expertHeaders_perMouse' : cell matrix (nMouse x nSess) headers (empties allowed)
%   'expertHeaders_global'   : cell matrix (nMouse x nSess) headers (empties allowed)
%
%   'ComputeAllSessions'       : true (default)  -> compute per-session axes for all sessions
%   'ComputeAllSessionsPooled' : true (default)  -> pooled across ALL sessions via concatYhat->SVD
%
%   'OrthMode'        : "GS" (default) | "none"
%       - applied to FINAL combined axis set per output (perMouse/global/pooled)
%       - within-group SVD axes are already orthogonal; GS only affects cross-group orthogonality
%
%   'SignFix'         : "projmean" (default) | "maxabs" | "none"
%   'SignFixWinSec'   : [0 4] (default)  % time window for projmean sign anchor
%   'SignFixTrialPolicy' : "GoNoGoByGroup" (default)
%       - groupName containing "nogo" -> NoGo trials
%       - otherwise -> Go trials
%
%   'Eps'             : 1e-10
%   'MultiDimGroups'  : {'GoToneOn','NoGoToneOn','ToneOffGo','ToneOffNoGo'}
%   'MultiDimK'       : 3
%   'PriorityNames'   : {'GoToneOn','NoGoToneOn','ToneOffGo','ToneOffNoGo','Lick'}
%   'Verbose'         : true
%   'PerSessionGS'    : true (default) -> store per-session GS axes in .A (else .A = .Araw_ord)
%
% OUTPUT (anchor struct)
%   .perMouse.axesByMouse{iMouse}   : axes from LAST provided expert session (no pooling per mouse)
%   .global.axes                    : pooled axes from LAST session per mouse (data pooling)
%   .allSessions.perSessionMat{j,s} : per-session axes struct (or [])
%   .allSessions.perSessionList(k)  : flat list of per-session structs
%   .allSessions.pooled             : pooled across ALL sessions (data pooling)
%
% NOTE:
%   - Per-session axes are computed per session (no pooling).
%   - Pooled anchors pool Yhat_g across sessions then SVD once per group.

% --- PRE-TRIM: allow headerC/glmRezC/trIdC to have trailing all-empty columns ---
[headerC, glmRezC, trIdC] = trimTrailingAllEmptyCols_(headerC, glmRezC, trIdC);

% -------------------- parse --------------------
p = inputParser;
p.addRequired('headerC', @(c) iscell(c));
p.addRequired('glmRezC', @(c) iscell(c) && isequal(size(c), size(headerC)));
p.addRequired('trIdC',   @(c) iscell(c) && isequal(size(c), size(headerC)));

p.addParameter('expertHeaders_perMouse', {}, @(c) iscell(c) || isstring(c));
p.addParameter('expertHeaders_global',   {}, @(c) iscell(c) || isstring(c));

p.addParameter('ComputeAllSessions', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('ComputeAllSessionsPooled', true, @(x)islogical(x)&&isscalar(x));

p.addParameter('OrthMode',"GS",@(s)ischar(s)||isstring(s));
p.addParameter('SignFix',"projmean",@(s)ischar(s)||isstring(s));
p.addParameter('SignFixWinSec',[0 4],@(x)isnumeric(x)&&numel(x)==2);
p.addParameter('SignFixTrialPolicy',"GoNoGoByGroup",@(s)ischar(s)||isstring(s));

p.addParameter('Eps',1e-10,@(x)isscalar(x)&&x>0);
p.addParameter('MultiDimGroups', {'GoToneOn','NoGoToneOn','ToneOffGo','ToneOffNoGo'}, @(c)iscell(c)||isstring(c));
p.addParameter('MultiDimK', 3, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('PriorityNames', {'GoToneOn','NoGoToneOn','ToneOffGo','ToneOffNoGo','Lick'}, @(c)iscell(c)||isstring(c));
p.addParameter('Verbose', true, @(x)islogical(x)&&isscalar(x));

p.addParameter('PerSessionGS', true, @(x)islogical(x)&&isscalar(x));

p.parse(headerC, glmRezC, trIdC, varargin{:});
opt = p.Results;
opt.MultiDimGroups = cellstr(string(opt.MultiDimGroups));
opt.PriorityNames  = cellstr(string(opt.PriorityNames));

% -------------------- normalize expert header matrices --------------------
hdrPM_mat = normalizeHeaderMatrix_(opt.expertHeaders_perMouse);
hdrG_mat  = normalizeHeaderMatrix_(opt.expertHeaders_global);

% -------------------- build lookup map from headerC -> glmRezC (+ trIdC) --------------------
[uniqHdr, glmFlatFirst, trFlatFirst] = buildHeaderLookup_(headerC, glmRezC, trIdC);

% ============================================================
% per-mouse (experts only): LAST provided session ONLY (no pooling)
% ============================================================
[nMousePM, nSessPM] = size(hdrPM_mat);

pm_foundMat     = false(nMousePM, nSessPM);
pm_glmMat       = cell(nMousePM, nSessPM);
pm_trIdMat      = cell(nMousePM, nSessPM);
pm_mouseIds     = cell(nMousePM, 1);
pm_axesByMouse  = cell(nMousePM, 1);
pm_headersUsed  = cell(nMousePM, 1);

for iM = 1:nMousePM
    hdrRow = hdrPM_mat(iM,:);
    pm_mouseIds{iM} = inferMouseIdFromRow_(hdrRow);

    usedHdr = strings(0,1);
    for iS = 1:nSessPM
        h = hdrRow{iS};
        if isempty(h), continue; end
        [g, t, ok] = fetchByHeaderFromLookup_(h, uniqHdr, glmFlatFirst, trFlatFirst);
        pm_foundMat(iM,iS) = ok;
        pm_glmMat{iM,iS}   = g;
        pm_trIdMat{iM,iS}  = t;
        if ok, usedHdr(end+1,1) = string(h); end %#ok<AGROW>
    end
    pm_headersUsed{iM} = cellstr(usedHdr);

    % pick last non-empty header in row (last-provided)
    hLast = "";
    for iS = nSessPM:-1:1
        if ~isempty(hdrRow{iS}), hLast = string(hdrRow{iS}); break; end
    end
    if strlength(hLast)==0, pm_axesByMouse{iM} = []; continue; end

    [gr, tr, ok] = fetchByHeaderFromLookup_(hLast, uniqHdr, glmFlatFirst, trFlatFirst);
    if ~ok, pm_axesByMouse{iM} = []; continue; end

    tdr = targetedDimRed_from_glmRez(gr, ...
        'trId', tr, ...
        'OrthMode', "none", ...   % per-session raw ordering; GS optional below
        'SignFix',  opt.SignFix, ...
        'SignFixWinSec', opt.SignFixWinSec, ...
        'SignFixTrialPolicy', opt.SignFixTrialPolicy, ...
        'Eps', opt.Eps, ...
        'ProjectWhichY', "Yz", ...
        'MultiDimGroups', opt.MultiDimGroups, ...
        'MultiDimK', opt.MultiDimK, ...
        'PriorityNames', opt.PriorityNames);

    pm_axesByMouse{iM} = packAxesOut_(tdr, opt.OrthMode, opt.Eps, char(hLast));
end

% ============================================================
% global pooling (experts only): LAST session per mouse then concatYhat->SVD
% ============================================================
[nMouseG, nSessG] = size(hdrG_mat);
g_foundMat = false(nMouseG, nSessG);
g_glmMat   = cell(nMouseG, nSessG);
g_trIdMat  = cell(nMouseG, nSessG);
hdrUsedGlobal = strings(0,1);

for iM = 1:nMouseG
    for iS = 1:nSessG
        h = hdrG_mat{iM,iS};
        if isempty(h), continue; end
        [g, t, ok] = fetchByHeaderFromLookup_(h, uniqHdr, glmFlatFirst, trFlatFirst);
        g_foundMat(iM,iS) = ok;
        g_glmMat{iM,iS}   = g;
        g_trIdMat{iM,iS}  = t;
        if ok, hdrUsedGlobal(end+1,1) = string(h); end %#ok<AGROW>
    end
end

if isempty(hdrUsedGlobal)
    warning('collectPerMouseAndGlobalGlmPrjAxes:NoGlobalExpertsFound', ...
        'None of expertHeaders_global matched headerC/glmRezC. anchor.global.axes will be empty.');
    globalAxes = [];
    hdrEndList = {};
else
    [glmEndList, trEndList, hdrEndList] = pickLastSessionPerMouse_(hdrG_mat, uniqHdr, glmFlatFirst, trFlatFirst);
    if isempty(glmEndList)
        warning('collectPerMouseAndGlobalGlmPrjAxes:NoEndSessions', ...
            'No global end sessions found. anchor.global.axes will be empty.');
        globalAxes = [];
    else
        globalAxes = buildPooledAxes_concatYhatSVD(glmEndList, trEndList, ...
            'headers', hdrEndList, ...
            'OrthMode', opt.OrthMode, ...
            'SignFix',  opt.SignFix, ...
            'SignFixWinSec', opt.SignFixWinSec, ...
            'SignFixTrialPolicy', opt.SignFixTrialPolicy, ...
            'Eps', opt.Eps, ...
            'MultiDimGroups', opt.MultiDimGroups, ...
            'MultiDimK', opt.MultiDimK, ...
            'PriorityNames', opt.PriorityNames, ...
            'Verbose', opt.Verbose);
    end
end

% ============================================================
% all sessions: per-session axes (no pooling), plus pooled across ALL sessions (optional)
% ============================================================
allSess = struct();
allSess.perSessionMat  = cell(size(headerC));
allSess.perSessionList = struct('j',{},'s',{},'header',{},'tdr',{},'trId',{}, ...
    'Araw_ord',{},'names_ord',{},'axisMeta_ord',{}, ...
    'A',{},'names',{},'keepAfterGS',{});
allSess.pooled = [];

if opt.ComputeAllSessions
    [perMat, perList, allHdrFlat, allGlmFlat, allTrFlat] = buildAllSessionAxes_(headerC, glmRezC, trIdC, opt);
    allSess.perSessionMat  = perMat;
    allSess.perSessionList = perList;

    if opt.ComputeAllSessionsPooled && ~isempty(allGlmFlat)
        allSess.pooled = buildPooledAxes_concatYhatSVD(allGlmFlat, allTrFlat, ...
            'headers', allHdrFlat, ...
            'OrthMode', opt.OrthMode, ...
            'SignFix',  opt.SignFix, ...
            'SignFixWinSec', opt.SignFixWinSec, ...
            'SignFixTrialPolicy', opt.SignFixTrialPolicy, ...
            'Eps', opt.Eps, ...
            'MultiDimGroups', opt.MultiDimGroups, ...
            'MultiDimK', opt.MultiDimK, ...
            'PriorityNames', opt.PriorityNames, ...
            'Verbose', opt.Verbose);
    end
end

% -------------------- pack output --------------------
anchor = struct();

anchor.perMouse = struct();
anchor.perMouse.headersMat    = hdrPM_mat;
anchor.perMouse.mouseIds      = pm_mouseIds;
anchor.perMouse.foundMat      = pm_foundMat;
anchor.perMouse.glmRezMat     = pm_glmMat;
anchor.perMouse.trIdMat       = pm_trIdMat;
anchor.perMouse.axesByMouse   = pm_axesByMouse;
anchor.perMouse.headersUsed   = pm_headersUsed;

anchor.global = struct();
anchor.global.headersMat      = hdrG_mat;
anchor.global.foundMat        = g_foundMat;
anchor.global.glmRezMat       = g_glmMat;
anchor.global.trIdMat         = g_trIdMat;
anchor.global.axes            = globalAxes;
anchor.global.headersUsed     = cellstr(hdrUsedGlobal);
anchor.global.endHeadersUsed  = hdrEndList;

anchor.allSessions = allSess;
anchor.opt = opt;

if opt.Verbose
    fprintf('[collectPerMouseAndGlobalGlmPrjAxes] perMouse=%d/%d | globalEnd=%d | allSessions=%d | allPooled=%d\n', ...
        sum(~cellfun(@isempty, pm_axesByMouse)), nMousePM, numel(hdrEndList), numel(allSess.perSessionList), ~isempty(allSess.pooled));
end

end

%% ======================================================================
% Helpers
% ======================================================================

function axesOut = packAxesOut_(tdr, orthMode, epsVal, headerStr)
axesOut = struct();
axesOut.Araw_ord      = tdr.Araw_ord;
axesOut.names_ord     = tdr.names_ord;
axesOut.axisMeta_ord  = tdr.axisMeta_ord;

if strcmpi(string(orthMode), "none")
    axesOut.A = tdr.Araw_ord;
    axesOut.keepAfterGS = true(1,size(tdr.Araw_ord,1));
    axesOut.names = tdr.names_ord;
    axesOut.axisMeta_final = tdr.axisMeta_ord;
else
    [A_gs, keep_gs] = gs_orth_rows(tdr.Araw_ord, epsVal);
    axesOut.A = A_gs;
    axesOut.keepAfterGS = keep_gs;
    axesOut.names = tdr.names_ord(keep_gs);
    axesOut.axisMeta_final = tdr.axisMeta_ord(keep_gs);
end

axesOut.session = struct('header', headerStr, 'tdr', tdr);
end

function [perMat, perList, hdrFlatOut, glmFlatOut, trFlatOut] = buildAllSessionAxes_(headerC, glmRezC, trIdC, opt)
[J,S] = size(headerC);
perMat = cell(J,S);

hdrFlatOut = {};
glmFlatOut = {};
trFlatOut  = {};

k = 0;
perList = struct('j',{},'s',{},'header',{},'tdr',{},'trId',{}, ...
    'Araw_ord',{},'names_ord',{},'axisMeta_ord',{}, ...
    'A',{},'names',{},'keepAfterGS',{});

for j = 1:J
    for s = 1:S
        gr = glmRezC{j,s};
        h  = headerC{j,s};
        tr = trIdC{j,s};

        if isempty(gr) || ~isstruct(gr) || isempty(h) || isempty(tr) || ~isstruct(tr)
            perMat{j,s} = [];
            continue;
        end

        try
            tdr = targetedDimRed_from_glmRez(gr, ...
                'trId', tr, ...
                'OrthMode', "none", ... % compute ordered raw axes; GS optional below
                'SignFix',  opt.SignFix, ...
                'SignFixWinSec', opt.SignFixWinSec, ...
                'SignFixTrialPolicy', opt.SignFixTrialPolicy, ...
                'Eps', opt.Eps, ...
                'ProjectWhichY', "Yz", ...
                'MultiDimGroups', opt.MultiDimGroups, ...
                'MultiDimK', opt.MultiDimK, ...
                'PriorityNames', opt.PriorityNames);

            sess = struct();
            sess.header       = h;
            sess.trId         = tr;
            sess.tdr          = tdr;
            sess.Araw_ord     = tdr.Araw_ord;
            sess.names_ord    = tdr.names_ord;
            sess.axisMeta_ord = tdr.axisMeta_ord;

            if opt.PerSessionGS
                [A_gs, keep_gs] = gs_orth_rows(tdr.Araw_ord, opt.Eps);
                sess.A            = A_gs;
                sess.keepAfterGS  = keep_gs;
                sess.names        = tdr.names_ord(keep_gs);
            else
                sess.A            = tdr.Araw_ord;
                sess.keepAfterGS  = true(1,size(tdr.Araw_ord,1));
                sess.names        = tdr.names_ord;
            end

            perMat{j,s} = sess;

            k = k + 1;
            perList(k).j = j;
            perList(k).s = s;
            perList(k).header = h;
            perList(k).trId = tr;
            perList(k).tdr = tdr;
            perList(k).Araw_ord = sess.Araw_ord;
            perList(k).names_ord = sess.names_ord;
            perList(k).axisMeta_ord = sess.axisMeta_ord;
            perList(k).A = sess.A;
            perList(k).names = sess.names;
            perList(k).keepAfterGS = sess.keepAfterGS;

            hdrFlatOut{end+1,1} = h;  %#ok<AGROW>
            glmFlatOut{end+1,1} = gr; %#ok<AGROW>
            trFlatOut{end+1,1}  = tr; %#ok<AGROW>

        catch ME
            perMat{j,s} = [];
            if opt.Verbose
                warning('buildAllSessionAxes:Fail', 'Skipping %d,%d (%s): %s', j, s, string(h), ME.message);
            end
        end
    end
end
end

function tdr = targetedDimRed_from_glmRez(glmRez, varargin)
% Per-session GLM-TDR:
%   For each group g: Yhat_g = Xz_g * beta_g; axis = top right singular vectors (in motif space).
% SignFix="projmean": flip sign so mean projection on (Go or NoGo) trials within winSec is positive.

p = inputParser;
p.addParameter('trId', struct(), @(x)isstruct(x) || isempty(x));

p.addParameter('OrthMode',"none",@(s)ischar(s)||isstring(s)); % (ignored here; kept for API consistency)
p.addParameter('SignFix',"projmean",@(s)ischar(s)||isstring(s));
p.addParameter('SignFixWinSec',[0 4],@(x)isnumeric(x)&&numel(x)==2);
p.addParameter('SignFixTrialPolicy',"GoNoGoByGroup",@(s)ischar(s)||isstring(s));

p.addParameter('Eps',1e-10,@(x)isscalar(x) && x>0);
p.addParameter('ProjectWhichY',"Yz",@(s)ischar(s)||isstring(s));

p.addParameter('MultiDimGroups', {'GoToneOn','NoGoToneOn','ToneOffGo','ToneOffNoGo'}, @(c)iscell(c)||isstring(c));
p.addParameter('MultiDimK', 3, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('PriorityNames', {'GoToneOn','NoGoToneOn','ToneOffGo','ToneOffNoGo','Lick'}, @(c)iscell(c)||isstring(c));
p.parse(varargin{:});
opt = p.Results;
opt.MultiDimGroups = cellstr(string(opt.MultiDimGroups));
opt.PriorityNames  = cellstr(string(opt.PriorityNames));

trId = opt.trId;

needFields = {'beta','X_design','muX','sdX','Yz','group','decBins'};
for f = needFields
    assert(isfield(glmRez,f{1}), 'glmRez missing required field: %s', f{1});
end
assert(isfield(glmRez.decBins,'time') && ~isempty(glmRez.decBins.time), 'glmRez.decBins.time required');

B    = glmRez.beta;        % [P x K]
Xraw = glmRez.X_design;    % [M x P]
muX  = glmRez.muX;         % [1 x P]
sdX  = glmRez.sdX;         % [1 x P]
Yz   = glmRez.Yz;          % [M x K]
grp  = glmRez.group;

[M, P]  = size(Xraw);
[Pb, K] = size(B);
assert(P==Pb, 'X_design cols (%d) != beta rows (%d)', P, Pb);
assert(all(size(Yz)==[M K]), 'Yz must be [M x K]');

Xz = (Xraw - muX) ./ sdX;
Xz(~isfinite(Xz)) = 0;

nW = numel(glmRez.decBins.time);
assert(rem(M,nW)==0, 'M=%d not divisible by nW=%d', M, nW);
N = M / nW;

G = local_num_groups_(grp);
groupNames = strings(1,G);
for g = 1:G
    gg = local_get_group_(grp, g);
    groupNames(g) = string(gg.name);
end

Araw_rows = [];
axisNames = {};
axisMeta  = struct('groupIdx',{},'groupName',{},'pcIdx',{},'sv',{},'expl',{});

for g = 1:G
    gg    = local_get_group_(grp, g);
    gName = string(gg.name);
    cols  = unique(round(gg.cols(:)'),'stable');
    cols  = cols(isfinite(cols) & cols>=1 & cols<=P);
    if isempty(cols), continue; end

    Yhat_g = Xz(:,cols) * B(cols,:); % [M x K]
    if norm(Yhat_g,'fro') < opt.Eps, continue; end

    nKeep = 1;
    if any(strcmpi(char(gName), opt.MultiDimGroups)), nKeep = opt.MultiDimK; end

    [~, Ssvd, V] = svd(Yhat_g, 'econ');
    r = size(V,2);
    nKeepEff = min(nKeep, r);

    svals = diag(Ssvd);
    denom = max(sum(svals.^2), opt.Eps);

    for pc = 1:nKeepEff
        v = V(:,pc)'; % 1 x K
        if norm(v) < opt.Eps, continue; end
        v = v / max(norm(v), opt.Eps);

        % if not projmean, can fix sign geometrically
        if ~strcmpi(string(opt.SignFix), "projmean")
            v = sign_fix_axis(v, opt.SignFix);
        end

        Araw_rows(end+1,:) = v; %#ok<AGROW>
        axisNames{end+1}   = char(gName + "_" + string(pc)); %#ok<AGROW>

        axisMeta(end+1).groupIdx  = g; %#ok<AGROW>
        axisMeta(end).groupName   = char(gName);
        axisMeta(end).pcIdx       = pc;
        axisMeta(end).sv          = svals(pc);
        axisMeta(end).expl        = (svals(pc)^2) / denom;
    end
end

assert(~isempty(Araw_rows), 'No valid axes constructed');

% order by PriorityNames then pcIdx
rankPerGroup = 1:G;
for i = 1:numel(opt.PriorityNames)
    gIdx = find(strcmpi(cellstr(groupNames), opt.PriorityNames{i}), 1, 'first');
    if ~isempty(gIdx), rankPerGroup(gIdx) = -1000 + i; end
end

nAxis = size(Araw_rows,1);
key1 = zeros(nAxis,1);
key2 = zeros(nAxis,1);
for a = 1:nAxis
    key1(a) = rankPerGroup(axisMeta(a).groupIdx);
    key2(a) = axisMeta(a).pcIdx;
end
[~, orderAxes] = sortrows([key1 key2],[1 2]);

Araw_ord     = Araw_rows(orderAxes,:);
names_ord    = axisNames(orderAxes);
axisMeta_ord = axisMeta(orderAxes);

% projmean sign fix (trial-type specific)
if strcmpi(string(opt.SignFix), "projmean")
    timeVec = glmRez.decBins.time(:);
    tMask = (timeVec >= opt.SignFixWinSec(1)) & (timeVec <= opt.SignFixWinSec(2));
    if any(tMask)
        switch lower(string(opt.ProjectWhichY))
            case "yz"
                Yscore = glmRez.Yz;
            case "ybig"
                assert(isfield(glmRez,'Ybig') && ~isempty(glmRez.Ybig), 'glmRez.Ybig missing/empty');
                Yscore = glmRez.Ybig;
            otherwise
                error('Unknown ProjectWhichY: %s', string(opt.ProjectWhichY));
        end

        for a = 1:size(Araw_ord,1)
            gName = string(axisMeta_ord(a).groupName);
            trialMask = local_trialMask_fromGroup_(trId, gName, N, opt.SignFixTrialPolicy);

            z = Yscore * Araw_ord(a,:)'; % [M x 1]
            Z = reshape(z, [N, nW]);     % [trial x time]
            score = mean(Z(trialMask, tMask), 'all', 'omitnan');
            if score < 0
                Araw_ord(a,:) = -Araw_ord(a,:);
            end
        end
    else
        warning('SignFixWinSec [%g %g] selects no timepoints; skipping projmean sign fix', opt.SignFixWinSec(1), opt.SignFixWinSec(2));
    end
end

% pack
tdr = struct();
tdr.M = M; tdr.P = P; tdr.K = K; tdr.nW = nW; tdr.N = N;
tdr.Xz = Xz;

tdr.Araw_ord     = Araw_ord;
tdr.names_ord    = names_ord;
tdr.axisMeta_ord = axisMeta_ord;
tdr.orderAxes    = orderAxes;

% optional projection (kept, because you’ve used it downstream)
switch lower(string(opt.ProjectWhichY))
    case "yz",   Yproj = glmRez.Yz;
    case "ybig", Yproj = glmRez.Ybig;
end
Zrows  = Yproj * Araw_ord';              % [M x nAxis]
Z      = reshape(Zrows, [N, nW, nAxis]); % [trial x time x axis]
tdr.Zrows = Zrows;
tdr.Z     = Z;

tdr.opt = opt;
end

function axesOut = buildPooledAxes_concatYhatSVD(glmRezList, trIdList, varargin)
% PRINCIPLED pooled anchor:
%   For each group g: pool Yhat_g rows across sessions, then SVD once -> shared motif axes.
%   Then optional SignFix (projmean) using each session’s Yz projections on those axes.
%   Then optional GS across groups (OrthMode).

p = inputParser;
p.addRequired('glmRezList', @(c) iscell(c) && ~isempty(c));
p.addRequired('trIdList',  @(c) iscell(c) && numel(c)==numel(glmRezList));

p.addParameter('headers', {}, @(c) isempty(c) || iscell(c) || isstring(c));
p.addParameter('OrthMode', "GS",   @(s)ischar(s)||isstring(s));
p.addParameter('SignFix',  "projmean", @(s)ischar(s)||isstring(s));
p.addParameter('SignFixWinSec', [0 4], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('SignFixTrialPolicy', "GoNoGoByGroup", @(s)ischar(s)||isstring(s));

p.addParameter('Eps', 1e-10, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('MultiDimGroups', {'GoToneOn','NoGoToneOn','ToneOffGo','ToneOffNoGo'}, @(c)iscell(c)||isstring(c));
p.addParameter('MultiDimK', 3, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('PriorityNames', {'GoToneOn','NoGoToneOn','ToneOffGo','ToneOffNoGo','Lick'}, @(c)iscell(c)||isstring(c));
p.addParameter('Verbose', true, @(x)islogical(x)&&isscalar(x));
p.parse(glmRezList, trIdList, varargin{:});
opt = p.Results;

opt.MultiDimGroups = cellstr(string(opt.MultiDimGroups));
opt.PriorityNames  = cellstr(string(opt.PriorityNames));

glmRezList = glmRezList(:);
trIdList   = trIdList(:);
S = numel(glmRezList);

% group list from first session (assumed consistent across sessions)
grp0 = glmRezList{1}.group;
G = local_num_groups_(grp0);

groupNames = strings(1,G);
groupCols  = cell(1,G);
for g = 1:G
    gg = local_get_group_(grp0, g);
    groupNames(g) = string(gg.name);
    groupCols{g}  = unique(round(gg.cols(:)'),'stable');
end

% pooled axes (rows)
Araw_rows = [];
axisNames = {};
axisMeta  = struct('groupIdx',{},'groupName',{},'pcIdx',{},'sv',{},'expl',{},'nSess',{});

for g = 1:G
    gName = groupNames(g);
    cols0 = groupCols{g};
    cols0 = cols0(isfinite(cols0));

    nKeep = 1;
    if any(strcmpi(char(gName), opt.MultiDimGroups)), nKeep = opt.MultiDimK; end

    Ypool = [];
    for s = 1:S
        gr = glmRezList{s};
        P = size(gr.X_design,2);

        cols = cols0(cols0>=1 & cols0<=P);
        if isempty(cols), continue; end

        Xz = (gr.X_design - gr.muX) ./ gr.sdX;
        Xz(~isfinite(Xz)) = 0;

        Yhat = Xz(:,cols) * gr.beta(cols,:); % [M x K]

        % equalize session contributions so long sessions don’t dominate
        Yhat = Yhat / max(norm(Yhat,'fro'), opt.Eps);

        Ypool = [Ypool; Yhat]; %#ok<AGROW>
    end

    if isempty(Ypool) || norm(Ypool,'fro') < opt.Eps
        continue;
    end

    [~, Ssvd, V] = svd(Ypool, 'econ'); % V: [K x r]
    r = size(V,2);
    nKeepEff = min(nKeep, r);

    svals = diag(Ssvd);
    denom = max(sum(svals.^2), opt.Eps);

    for pc = 1:nKeepEff
        v = V(:,pc)'; % 1 x K
        v = v / max(norm(v), opt.Eps);

        % sign fix
        if strcmpi(string(opt.SignFix),"projmean")
            sc = nan(S,1);
            for ss = 1:S
                sc(ss) = local_projmean_score_glm_(glmRezList{ss}, trIdList{ss}, v, gName, opt.SignFixWinSec, opt.SignFixTrialPolicy, opt.Eps);
            end
            if mean(sc,'omitnan') < 0, v = -v; end
        elseif strcmpi(string(opt.SignFix),"maxabs")
            v = sign_fix_axis(v, "maxabs");
        end

        Araw_rows(end+1,:) = v; %#ok<AGROW>
        axisNames{end+1}   = char(gName + "_" + string(pc)); %#ok<AGROW>

        axisMeta(end+1).groupIdx = g; %#ok<AGROW>
        axisMeta(end).groupName  = char(gName);
        axisMeta(end).pcIdx      = pc;
        axisMeta(end).sv         = svals(pc);
        axisMeta(end).expl       = (svals(pc)^2)/denom;
        axisMeta(end).nSess      = S;
    end
end

if isempty(Araw_rows)
    axesOut = [];
    return;
end

% order pooled axes by PriorityNames then pcIdx
rankPerGroup = 1:G;
for i = 1:numel(opt.PriorityNames)
    gi = find(strcmpi(cellstr(groupNames), opt.PriorityNames{i}), 1, 'first');
    if ~isempty(gi), rankPerGroup(gi) = -1000 + i; end
end

nAxis = size(Araw_rows,1);
key1 = nan(nAxis,1); key2 = nan(nAxis,1);
for a = 1:nAxis
    key1(a) = rankPerGroup(axisMeta(a).groupIdx);
    key2(a) = axisMeta(a).pcIdx;
end
[~, ord] = sortrows([key1 key2],[1 2]);

Araw_ord     = Araw_rows(ord,:);
names_ord    = axisNames(ord);
axisMeta_ord = axisMeta(ord);

% optional GS across groups
if strcmpi(string(opt.OrthMode),"none")
    A = Araw_ord;
    keep = true(1,size(Araw_ord,1));
    namesFinal = names_ord;
    metaFinal  = axisMeta_ord;
else
    [A, keep] = gs_orth_rows(Araw_ord, opt.Eps);
    namesFinal = names_ord(keep);
    metaFinal  = axisMeta_ord(keep);
end

axesOut = struct();
axesOut.Araw_ord      = Araw_ord;
axesOut.names_ord     = names_ord;
axesOut.axisMeta_ord  = axisMeta_ord;

axesOut.A             = A;
axesOut.keepAfterGS   = keep;
axesOut.names         = namesFinal;
axesOut.axisMeta_final= metaFinal;

axesOut.headers       = opt.headers;
axesOut.opt           = opt;

if opt.Verbose
    fprintf('[buildPooledAxes_concatYhatSVD] axes=%d (kept after GS=%d) from %d sessions\n', ...
        size(Araw_ord,1), size(A,1), S);
end
end

function score = local_projmean_score_glm_(glmRez, trId, axisRow, groupName, winSec, policy, epsVal)
M = size(glmRez.Yz,1);
nW = numel(glmRez.decBins.time);
if rem(M,nW)~=0, score = NaN; return; end
N = M / nW;

timeVec = glmRez.decBins.time(:);
tMask = (timeVec >= winSec(1)) & (timeVec <= winSec(2));
if ~any(tMask), score = NaN; return; end

trialMask = local_trialMask_fromGroup_(trId, string(groupName), N, policy);

z = glmRez.Yz * axisRow(:);
Z = reshape(z, [N, nW]);
score = mean(Z(trialMask, tMask), 'all', 'omitnan');

if ~isfinite(score) || abs(score) < epsVal
    % caller uses mean(...,'omitnan'), so NaNs just reduce influence
end
end

function [uniqHdr, glmFlatFirst, trFlatFirst] = buildHeaderLookup_(headerC, glmRezC, trIdC)
hdrFlat = headerC(:);
glmFlat = glmRezC(:);
trFlat  = trIdC(:);

isHdr = ~cellfun(@isempty, hdrFlat);
hdrFlat = hdrFlat(isHdr);
glmFlat = glmFlat(isHdr);
trFlat  = trFlat(isHdr);

hdrFlatS = string(hdrFlat);
[uniqHdr, ~, ic] = unique(hdrFlatS, 'stable');

if numel(uniqHdr) < numel(hdrFlatS)
    counts = accumarray(ic, 1);
    dup = uniqHdr(counts > 1);
    warning('collectPerMouseAndGlobalGlmPrjAxes:DuplicateHeaders', ...
        'Duplicate headers in headerC; matching uses first occurrence. Example(s): %s', ...
        strjoin(cellstr(dup(1:min(5,end))), ', '));
end

glmFlatFirst = cell(numel(uniqHdr),1);
trFlatFirst  = cell(numel(uniqHdr),1);
for i = 1:numel(uniqHdr)
    ii = find(hdrFlatS == uniqHdr(i), 1, 'first');
    glmFlatFirst{i} = glmFlat{ii};
    trFlatFirst{i}  = trFlat{ii};
end
end

function [glmRez, trId, found] = fetchByHeaderFromLookup_(hdr, uniqHdr, glmFlatFirst, trFlatFirst)
hdr = string(hdr);
j = find(uniqHdr == hdr, 1, 'first');
if isempty(j), glmRez = []; trId = []; found = false; return; end
glmRez = glmFlatFirst{j};
trId   = trFlatFirst{j};
found  = ~isempty(glmRez) && isstruct(glmRez) && ~isempty(trId) && isstruct(trId);
end

function hdrMat = normalizeHeaderMatrix_(hdrIn)
if isempty(hdrIn), hdrMat = cell(0,0); return; end
if isstring(hdrIn), hdrIn = cellstr(hdrIn); end
if iscell(hdrIn) && isvector(hdrIn), hdrIn = hdrIn(:); end

hdrMat = cell(size(hdrIn));
for i = 1:numel(hdrIn)
    x = hdrIn{i};
    if isempty(x)
        hdrMat{i} = [];
    elseif isstring(x) || ischar(x)
        s = char(string(x));
        if strlength(string(s))==0, hdrMat{i} = []; else, hdrMat{i} = s; end
    else
        hdrMat{i} = [];
    end
end
end

function mouseId = inferMouseIdFromRow_(hdrRow)
mouseId = '';
for k = 1:numel(hdrRow)
    h = hdrRow{k};
    if isempty(h), continue; end
    tok = regexp(string(h), '(m\d{3,5})', 'tokens', 'once');
    if ~isempty(tok)
        mouseId = char(string(tok{1}));
        return;
    end
end
end

function [glmEndList, trEndList, hdrEndList] = pickLastSessionPerMouse_(hdrMat, uniqHdr, glmFlatFirst, trFlatFirst)
[nM, nS] = size(hdrMat);
glmEndList = {};
trEndList  = {};
hdrEndList = {};

for iM = 1:nM
    hLast = "";
    for iS = nS:-1:1
        if ~isempty(hdrMat{iM,iS})
            hLast = string(hdrMat{iM,iS});
            break;
        end
    end
    if strlength(hLast)==0, continue; end

    [gr, tr, ok] = fetchByHeaderFromLookup_(hLast, uniqHdr, glmFlatFirst, trFlatFirst);
    if ~ok, continue; end

    glmEndList{end+1,1} = gr; %#ok<AGROW>
    trEndList{end+1,1}  = tr; %#ok<AGROW>
    hdrEndList{end+1,1} = char(hLast); %#ok<AGROW>
end
end

% -------------------- small utilities --------------------
function G = local_num_groups_(grp)
if iscell(grp)
    G = numel(grp);
elseif isstruct(grp)
    G = numel(grp);
else
    error('glmRez.group must be cell or struct array, got: %s', class(grp));
end
end

function gg = local_get_group_(grp, g)
if iscell(grp)
    gg = grp{g};
else
    gg = grp(g);
end
assert(isstruct(gg) && isfield(gg,'name') && isfield(gg,'cols'), ...
    'glmRez.group(%d) must have fields .name and .cols', g);
end

function v = sign_fix_axis(v, mode)
mode = lower(string(mode));
switch mode
    case "maxabs"
        [~,idx] = max(abs(v));
        if v(idx) < 0, v = -v; end
    otherwise
        % "none" or unknown -> no-op
end
end

function [Q, keepIdx] = gs_orth_rows(A, epsVal)
Q = [];
keepIdx = false(1,size(A,1));
for i = 1:size(A,1)
    v = A(i,:);
    if norm(v) < epsVal, continue; end
    for j = 1:size(Q,1)
        v = v - (v*Q(j,:)') * Q(j,:);
    end
    nv = norm(v);
    if nv < epsVal, continue; end
    Q(end+1,:) = v / nv; %#ok<AGROW>
    keepIdx(i) = true;
end
end

function trialMask = local_trialMask_fromGroup_(trId, groupName, N, policy)
policy = string(policy);
groupName = string(groupName);

trialMask = true(N,1);

want = "all";
if strcmpi(policy, "GoNoGoByGroup")
    if contains(lower(groupName), "nogo")
        want = "nogo";
    else
        want = "go";
    end
end

if want == "all" || isempty(trId) || ~isstruct(trId)
    return;
end

if want == "go"
    trialMask = local_makeMask_(trId, 'goI', N);
else
    trialMask = local_makeMask_(trId, 'nogoI', N);
end

if ~any(trialMask)
    warning('local_trialMask_fromGroup_:EmptyMask', 'Empty trial mask for %s; using all trials', groupName);
    trialMask = true(N,1);
end
end

function m = local_makeMask_(trId, field, N)
m = true(N,1);
if ~isfield(trId, field) || isempty(trId.(field))
    warning('local_makeMask_:MissingField', 'trId.%s missing/empty; using all trials', field);
    return;
end
x = trId.(field);

if islogical(x)
    x = x(:);
    if numel(x) ~= N
        warning('local_makeMask_:BadLength', 'trId.%s length %d != N=%d; using all trials', field, numel(x), N);
        return;
    end
    m = x;
else
    idx = x(:);
    idx = idx(isfinite(idx));
    idx = unique(round(idx));
    idx = idx(idx>=1 & idx<=N);
    m = false(N,1);
    m(idx) = true;
end
end

function [headerC, glmRezC, trIdC] = trimTrailingAllEmptyCols_(headerC, glmRezC, trIdC)
% Trim trailing columns that are entirely empty ([]) in whichever arrays have them.
% This prevents size-mismatch errors when one input has an extra empty column at the end.

assert(iscell(headerC) && iscell(glmRezC) && iscell(trIdC), 'Inputs must be cell arrays.');
[J1,S1] = size(headerC);
[J2,S2] = size(glmRezC);
[J3,S3] = size(trIdC);

% If row counts disagree, don't guess.
if ~(J1==J2 && J1==J3)
    error('Row counts differ: headerC=%dx%d, glmRezC=%dx%d, trIdC=%dx%d', J1,S1,J2,S2,J3,S3);
end

Smax = max([S1 S2 S3]);

% Pad smaller ones with empty so we can test columns uniformly
if S1 < Smax, headerC(:,end+1:Smax) = {[]}; end
if S2 < Smax, glmRezC(:,end+1:Smax) = {[]}; end
if S3 < Smax, trIdC(:,end+1:Smax)   = {[]}; end

% Find last column that has ANY non-empty content in ANY of the three inputs
keepLast = 0;
for s = 1:Smax
    anyNonEmpty = any(~cellfun(@isempty, headerC(:,s))) || ...
                  any(~cellfun(@isempty, glmRezC(:,s)))  || ...
                  any(~cellfun(@isempty, trIdC(:,s)));
    if anyNonEmpty
        keepLast = s;
    end
end

if keepLast == 0
    % all empty everywhere; keep a 0-column shape consistent with original rows
    headerC = cell(J1,0);
    glmRezC = cell(J1,0);
    trIdC   = cell(J1,0);
else
    headerC = headerC(:,1:keepLast);
    glmRezC = glmRezC(:,1:keepLast);
    trIdC   = trIdC(:,1:keepLast);
end
end