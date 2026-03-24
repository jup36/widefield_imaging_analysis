function anchor = collectPerMouseAndGlobalGlmPrjAxes_maxabs(headerC, glmRezC, varargin)
%collectPerMouseAndGlobalGlmPrjAxes  Build expert per-mouse + expert global GLM-TDR axes,
%AND (NEW) build axes for EVERY session (per-session), optionally pooled-all-sessions axes.
%
% DROP-IN: replaces your current collectPerMouseAndGlobalGlmPrjAxes.m
%
% REQUIRED INPUTS
%   headerC : cell (J x S) session headers
%   glmRezC : cell (J x S) glmRez structs (same size as headerC)
%
% NAME-VALUE
%   'expertHeaders_perMouse' : cell matrix (nMouse x nSess) headers (empties allowed)
%   'expertHeaders_global'   : cell matrix (nMouse x nSess) headers (empties allowed)
%
%   'ComputeAllSessions'      : true (default)
%   'ComputeAllSessionsPooled': true (default)
%
%   % Passed to targetedDimRed_from_glmRez / buildAnchorAxes_from_glmRezList:
%   'OrthMode'        : "GS" (default) | "none"      (for POOLED axes)
%   'SignFix'         : "maxabs" (default) | "none"
%   'Eps'             : 1e-10
%   'MultiDimGroups'  : {'GoToneOn','NoGoToneOn','ToneOffGo','ToneOffNoGo'}
%   'MultiDimK'       : 3
%   'PriorityNames'   : {'GoToneOn','NoGoToneOn','ToneOffGo','ToneOffNoGo','Lick'}
%   'Verbose'         : true
%
%   % NEW:
%   'PerSessionGS'    : true (default) -> store per-session GS axes as .A (in addition to .Araw_ord)
%     (Backward compat: 'SavePerSessionGS' also accepted by buildAnchorAxes_from_glmRezList)
%
% OUTPUT (anchor struct)
%   .perMouse.axesByMouse{iMouse}   : pooled expert axes for that mouse
%   .global.axes                    : pooled expert axes across mice
%   .allSessions.perSessionMat{j,s} : per-session axes struct (or [])
%   .allSessions.perSessionList(k)  : flat list of per-session structs
%   .allSessions.pooled             : pooled axes across ALL sessions (optional)

% -------------------- parse --------------------
p = inputParser;
p.addRequired('headerC', @(c) iscell(c));
p.addRequired('glmRezC', @(c) iscell(c) && isequal(size(c), size(headerC)));

p.addParameter('expertHeaders_perMouse', {}, @(c) iscell(c) || isstring(c));
p.addParameter('expertHeaders_global',   {}, @(c) iscell(c) || isstring(c));

p.addParameter('ComputeAllSessions', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('ComputeAllSessionsPooled', true, @(x)islogical(x)&&isscalar(x));

% pass-through options
p.addParameter('OrthMode',"GS",@(s)ischar(s)||isstring(s));      % pooled only
p.addParameter('SignFix',"maxabs",@(s)ischar(s)||isstring(s));
p.addParameter('Eps',1e-10,@(x)isscalar(x)&&x>0);
p.addParameter('MultiDimGroups', {'GoToneOn','NoGoToneOn','ToneOffGo','ToneOffNoGo'}, @(c)iscell(c)||isstring(c));
p.addParameter('MultiDimK', 3, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('PriorityNames', {'GoToneOn','NoGoToneOn','ToneOffGo','ToneOffNoGo','Lick'}, @(c)iscell(c)||isstring(c));
p.addParameter('Verbose', true, @(x)islogical(x)&&isscalar(x));

% NEW: per-session GS storage
p.addParameter('PerSessionGS', true, @(x)islogical(x)&&isscalar(x));

p.parse(headerC, glmRezC, varargin{:});
opt = p.Results;

% -------------------- normalize expert header matrices --------------------
hdrPM_mat = normalizeHeaderMatrix_(opt.expertHeaders_perMouse);
hdrG_mat  = normalizeHeaderMatrix_(opt.expertHeaders_global);

% -------------------- build lookup map from headerC -> glmRezC --------------------
[uniqHdr, glmFlatFirst] = buildHeaderLookup_(headerC, glmRezC);

% -------------------- per-mouse pooling (experts only) --------------------
[nMousePM, nSessPM] = size(hdrPM_mat);
pm_foundMat     = false(nMousePM, nSessPM);
pm_glmMat       = cell(nMousePM, nSessPM);
pm_mouseIds     = cell(nMousePM, 1);
pm_axesByMouse  = cell(nMousePM, 1);
pm_headersUsed  = cell(nMousePM, 1);

for iM = 1:nMousePM
    hdrRow = hdrPM_mat(iM,:);
    pm_mouseIds{iM} = inferMouseIdFromRow_(hdrRow);

    glmPool = {};
    usedHdr = strings(0,1);

    for iS = 1:nSessPM
        h = hdrRow{iS};
        if isempty(h), continue; end

        [g, ok] = fetchByHeaderFromLookup_(h, uniqHdr, glmFlatFirst);
        pm_foundMat(iM,iS) = ok;
        pm_glmMat{iM,iS}   = g;

        if ok
            glmPool{end+1,1} = g; %#ok<AGROW>
            usedHdr(end+1,1) = string(h); %#ok<AGROW>
        end
    end

    pm_headersUsed{iM} = cellstr(usedHdr);

    if isempty(glmPool)
        pm_axesByMouse{iM} = [];
        continue;
    end

    pm_axesByMouse{iM} = buildAnchorAxes_from_glmRezList(glmPool, ...
        'headers', cellstr(usedHdr), ...
        'OrthMode', opt.OrthMode, ...
        'SignFix',  opt.SignFix, ...
        'Eps',      opt.Eps, ...
        'MultiDimGroups', cellstr(string(opt.MultiDimGroups)), ...
        'MultiDimK', opt.MultiDimK, ...
        'PriorityNames', cellstr(string(opt.PriorityNames)), ...
        'Verbose',  opt.Verbose, ...
        'PerSessionGS', opt.PerSessionGS);
end

% -------------------- global pooling (experts only) --------------------
[nMouseG, nSessG] = size(hdrG_mat);
g_foundMat = false(nMouseG, nSessG);
g_glmMat   = cell(nMouseG, nSessG);

glmPoolGlobal = {};
hdrUsedGlobal = strings(0,1);

for iM = 1:nMouseG
    for iS = 1:nSessG
        h = hdrG_mat{iM,iS};
        if isempty(h), continue; end

        [g, ok] = fetchByHeaderFromLookup_(h, uniqHdr, glmFlatFirst);
        g_foundMat(iM,iS) = ok;
        g_glmMat{iM,iS}   = g;

        if ok
            glmPoolGlobal{end+1,1} = g; %#ok<AGROW>
            hdrUsedGlobal(end+1,1) = string(h); %#ok<AGROW>
        end
    end
end

if isempty(glmPoolGlobal)
    warning('collectPerMouseAndGlobalGlmPrjAxes:NoGlobalExpertsFound', ...
        'None of expertHeaders_global matched headerC/glmRezC. anchor.global.axes will be empty.');
    globalAxes = [];
else
    globalAxes = buildAnchorAxes_from_glmRezList(glmPoolGlobal, ...
        'headers', cellstr(hdrUsedGlobal), ...
        'OrthMode', opt.OrthMode, ...
        'SignFix',  opt.SignFix, ...
        'Eps',      opt.Eps, ...
        'MultiDimGroups', cellstr(string(opt.MultiDimGroups)), ...
        'MultiDimK', opt.MultiDimK, ...
        'PriorityNames', cellstr(string(opt.PriorityNames)), ...
        'Verbose',  opt.Verbose, ...
        'PerSessionGS', opt.PerSessionGS);
end

% ============================================================
% NEW: axes for EVERY session (per-session), not just experts
% ============================================================
allSess = struct();
allSess.perSessionMat  = cell(size(headerC));
allSess.perSessionList = struct('j',{},'s',{},'header',{},'tdr',{}, ...
    'Araw_ord',{},'names_ord',{},'axisMeta_ord',{}, ...
    'A',{},'names',{},'keepAfterGS',{});
allSess.pooled         = [];

if opt.ComputeAllSessions
    [perMat, perList, allHdrFlat, allGlmFlat] = buildAllSessionAxes_(headerC, glmRezC, opt);

    allSess.perSessionMat  = perMat;
    allSess.perSessionList = perList;

    if opt.ComputeAllSessionsPooled && ~isempty(allGlmFlat)
        allSess.pooled = buildAnchorAxes_from_glmRezList(allGlmFlat, ...
            'headers', allHdrFlat, ...
            'OrthMode', opt.OrthMode, ...
            'SignFix',  opt.SignFix, ...
            'Eps',      opt.Eps, ...
            'MultiDimGroups', cellstr(string(opt.MultiDimGroups)), ...
            'MultiDimK', opt.MultiDimK, ...
            'PriorityNames', cellstr(string(opt.PriorityNames)), ...
            'Verbose',  opt.Verbose, ...
            'PerSessionGS', opt.PerSessionGS);
    end
end

% -------------------- pack output --------------------
anchor = struct();

anchor.perMouse = struct();
anchor.perMouse.headersMat    = hdrPM_mat;
anchor.perMouse.mouseIds      = pm_mouseIds;
anchor.perMouse.foundMat      = pm_foundMat;
anchor.perMouse.glmRezMat     = pm_glmMat;
anchor.perMouse.axesByMouse   = pm_axesByMouse;
anchor.perMouse.headersUsed   = pm_headersUsed;

anchor.global = struct();
anchor.global.headersMat      = hdrG_mat;
anchor.global.foundMat        = g_foundMat;
anchor.global.glmRezMat       = g_glmMat;
anchor.global.axes            = globalAxes;
anchor.global.headersUsed     = cellstr(hdrUsedGlobal);

anchor.allSessions = allSess;
anchor.opt = opt;

if opt.Verbose
    fprintf('[anchorMatch] perMouse pooled=%d/%d mice | global pooled=%d sessions | allSessions per-session=%d\n', ...
        sum(~cellfun(@isempty, pm_axesByMouse)), nMousePM, numel(glmPoolGlobal), numel(anchor.allSessions.perSessionList));
end

end

%% ======================================================================
% Helpers
% ======================================================================

function [perMat, perList, hdrFlatOut, glmFlatOut] = buildAllSessionAxes_(headerC, glmRezC, opt)
% Build targetedDimRed_from_glmRez for EVERY non-empty glmRezC{j,s}
% and also store per-session GS orthonormalized axes as .A.

[J,S] = size(headerC);
perMat = cell(J,S);

hdrFlatOut = {};
glmFlatOut = {};

k = 0;
perList = struct('j',{},'s',{},'header',{},'tdr',{}, ...
    'Araw_ord',{},'names_ord',{},'axisMeta_ord',{}, ...
    'A',{},'names',{},'keepAfterGS',{});

for j = 1:J
    for s = 1:S
        gr = glmRezC{j,s};
        h  = headerC{j,s};

        if isempty(gr) || ~isstruct(gr) || isempty(h)
            perMat{j,s} = [];
            continue;
        end

        try
            tdr = targetedDimRed_from_glmRez(gr, ...
                'OrthMode', "none", ...    % get ordered RAW axes
                'SignFix',  opt.SignFix, ...
                'Eps',      opt.Eps, ...
                'ProjectWhichY', "Yz", ...
                'MultiDimGroups', cellstr(string(opt.MultiDimGroups)), ...
                'MultiDimK', opt.MultiDimK, ...
                'PriorityNames', cellstr(string(opt.PriorityNames)));

            sess = struct();
            sess.header       = h;
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
            perList(k).tdr = tdr;
            perList(k).Araw_ord = sess.Araw_ord;
            perList(k).names_ord = sess.names_ord;
            perList(k).axisMeta_ord = sess.axisMeta_ord;
            perList(k).A = sess.A;
            perList(k).names = sess.names;
            perList(k).keepAfterGS = sess.keepAfterGS;

            hdrFlatOut{end+1,1} = h;  %#ok<AGROW>
            glmFlatOut{end+1,1} = gr; %#ok<AGROW>

        catch ME
            perMat{j,s} = [];
            if opt.Verbose
                warning('buildAllSessionAxes:Fail', 'Skipping %d,%d (%s): %s', j, s, string(h), ME.message);
            end
        end
    end
end
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

function [uniqHdr, glmFlatFirst] = buildHeaderLookup_(headerC, glmRezC)
hdrFlat = headerC(:);
glmFlat = glmRezC(:);

isHdr = ~cellfun(@isempty, hdrFlat);
hdrFlat = hdrFlat(isHdr);
glmFlat = glmFlat(isHdr);

hdrFlatS = string(hdrFlat);
[uniqHdr, ~, ic] = unique(hdrFlatS, 'stable');

if numel(uniqHdr) < numel(hdrFlatS)
    counts = accumarray(ic, 1);
    dup = uniqHdr(counts > 1);
    warning('collectPerMouseAndGlobalGlmPrjAxes:DuplicateHeaders', ...
        'Duplicate headers detected in headerC (matching uses first occurrence). Example(s): %s', ...
        strjoin(cellstr(dup(1:min(5,end))), ', '));
end

glmFlatFirst = cell(numel(uniqHdr),1);
for i = 1:numel(uniqHdr)
    ii = find(hdrFlatS == uniqHdr(i), 1, 'first');
    glmFlatFirst{i} = glmFlat{ii};
end
end

function [glmRez, found] = fetchByHeaderFromLookup_(hdr, uniqHdr, glmFlatFirst)
hdr = string(hdr);
j = find(uniqHdr == hdr, 1, 'first');
if isempty(j), glmRez = []; found = false; return; end
glmRez = glmFlatFirst{j};
found = ~isempty(glmRez) && isstruct(glmRez);
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

function axesOut = buildAnchorAxes_from_glmRezList(glmRezList, varargin)
%BUILDANCHORAXES_FROM_GLMREZLIST  Build pooled ("anchor") GLM-TDR axes from many glmRez,
%and ALSO store per-session axes including optional per-session GS.
%
% axesOut = buildAnchorAxes_from_glmRezList(glmRezList, 'Name', value, ...)
%
% Strategy:
%   1) For each session: run targetedDimRed_from_glmRez(gr, OrthMode="none") to obtain Araw_ord, names_ord.
%   2) Optionally GS-orthonormalize per-session Araw_ord -> session(i).A (PerSessionGS).
%   3) Pool across sessions per axis-name:
%        - gather that axis vector from each session (Araw_ord row)
%        - normalize, sign-align to first available
%        - average, renormalize, sign-fix
%   4) Order pooled axes by PriorityNames (group rank) then pcIdx.
%   5) Optionally GS on pooled Araw_ord -> axesOut.A (OrthMode).
%
% OUTPUT (axesOut)
%   .Araw_rows      : pooled, unordered [nAxis x K]
%   .names_raw      : pooled, unordered names
%   .Araw_ord       : pooled, ordered pre-GS
%   .names_ord      : pooled, ordered names pre-GS
%   .axisMeta_ord   : pooled, ordered meta
%   .A              : pooled final axes after optional GS
%   .names          : pooled final names after optional GS keep
%   .axisMeta_final : pooled final meta after optional GS keep
%   .keepAfterGS    : pooled keep mask after GS
%
%   .session(i)     : per-session package (only successful sessions)
%     .header
%     .tdr               (raw tdr output)
%     .Araw_ord
%     .names_ord
%     .axisMeta_ord
%     .A                 (per-session GS result if PerSessionGS=true else ==Araw_ord)
%     .names             (names after per-session keep mask)
%     .keepAfterGS
%
%   .tdrBySession   : cell(S,1) (includes failures as [])
%   .sessionKeep    : logical(S,1)
%   .headers        : headers for successful sessions (maps to .session)
%
% Drop-in: self-contained helper functions included at end.

% -------------------- parse --------------------
p = inputParser;
p.addRequired('glmRezList', @(c) iscell(c) && ~isempty(c));

p.addParameter('headers', {}, @(c) isempty(c) || iscell(c) || isstring(c));

% pooled-axis options
p.addParameter('OrthMode', "GS",   @(s)ischar(s)||isstring(s));   % pooled GS only
p.addParameter('SignFix',  "maxabs",@(s)ischar(s)||isstring(s));
p.addParameter('Eps',      1e-10,  @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('MultiDimGroups', {'GoToneOn','NoGoToneOn','ToneOffGo','ToneOffNoGo'}, @(c)iscell(c)||isstring(c));
p.addParameter('MultiDimK', 3, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('PriorityNames', {'GoToneOn','NoGoToneOn','ToneOffGo','ToneOffNoGo','Lick'}, @(c)iscell(c)||isstring(c));
p.addParameter('Verbose', true, @(x)islogical(x)&&isscalar(x));

% per-session GS option
p.addParameter('PerSessionGS', true, @(x)islogical(x)&&isscalar(x));

p.parse(glmRezList, varargin{:});
opt = p.Results;

opt.MultiDimGroups = cellstr(string(opt.MultiDimGroups));
opt.PriorityNames  = cellstr(string(opt.PriorityNames));

glmRezList = glmRezList(:);
S = numel(glmRezList);

% -------------------- headers normalize --------------------
hdrList = opt.headers;
if isstring(hdrList), hdrList = cellstr(hdrList); end
if isempty(hdrList)
    hdrList = cell(S,1);
else
    hdrList = hdrList(:);
    if numel(hdrList) ~= S
        error('buildAnchorAxes_from_glmRezList:BadHeaders', ...
            'headers must match glmRezList length (%d).', S);
    end
end

% -------------------- run per-session TDR (NO per-session GS inside) --------------------
tdrBySession = cell(S,1);
sessionKeep  = false(S,1);

for s = 1:S
    gr = glmRezList{s};
    if isempty(gr) || ~isstruct(gr)
        tdrBySession{s} = [];
        continue;
    end

    try
        tdrBySession{s} = targetedDimRed_from_glmRez(gr, ...
            'OrthMode',       "none", ...   % IMPORTANT: per-session raw ordered axes
            'SignFix',        opt.SignFix, ...
            'Eps',            opt.Eps, ...
            'ProjectWhichY',  "Yz", ...     % irrelevant for axis extraction
            'MultiDimGroups', opt.MultiDimGroups, ...
            'MultiDimK',      opt.MultiDimK, ...
            'PriorityNames',  opt.PriorityNames);

        sessionKeep(s) = true;

    catch ME
        if opt.Verbose
            warning('buildAnchorAxes_from_glmRezList:SessionFail', ...
                'Skipping session %d/%d (targetedDimRed_from_glmRez failed): %s', s, S, ME.message);
        end
        sessionKeep(s) = false;
        tdrBySession{s} = [];
    end
end

tdrGood = tdrBySession(sessionKeep);
hdrGood = hdrList(sessionKeep);

if isempty(tdrGood)
    error('buildAnchorAxes_from_glmRezList:NoValidSessions', ...
        'No glmRez entries produced valid targetedDimRed outputs.');
end

% -------------------- per-session packaging (includes per-session A) --------------------
sess = repmat(struct( ...
    'header',       [], ...
    'tdr',          [], ...
    'Araw_ord',     [], ...
    'names_ord',    [], ...
    'axisMeta_ord', [], ...
    'A',            [], ...
    'names',        [], ...
    'keepAfterGS',  []), numel(tdrGood), 1);

for i = 1:numel(tdrGood)
    tdr = tdrGood{i};

    sess(i).header       = hdrGood{i};
    sess(i).tdr          = tdr;

    sess(i).Araw_ord     = tdr.Araw_ord;
    sess(i).names_ord    = tdr.names_ord;
    sess(i).axisMeta_ord = tdr.axisMeta_ord;

    if opt.PerSessionGS
        [A_gs, keep_gs] = gs_orth_rows(tdr.Araw_ord, opt.Eps);
        sess(i).A           = A_gs;
        sess(i).keepAfterGS = keep_gs;
        sess(i).names       = tdr.names_ord(keep_gs);
    else
        sess(i).A           = tdr.Araw_ord;
        sess(i).keepAfterGS = true(1,size(tdr.Araw_ord,1));
        sess(i).names       = tdr.names_ord;
    end
end

% -------------------- union of axis names across sessions --------------------
namesAll = {};
for i = 1:numel(tdrGood)
    namesAll = [namesAll, tdrGood{i}.names_ord]; %#ok<AGROW>
end
namesAll = unique(namesAll, 'stable');

K = size(tdrGood{1}.Araw_ord, 2);

% -------------------- pool each axis name --------------------
A_pool = nan(numel(namesAll), K);
meta   = struct('groupName',{},'pcIdx',{},'nSessAvail',{},'nSessUsed',{});

for a = 1:numel(namesAll)
    nm = string(namesAll{a});

    Vc = {};
    for i = 1:numel(tdrGood)
        tdr = tdrGood{i};
        idx = find(strcmp(string(tdr.names_ord), nm), 1, 'first');
        if isempty(idx), continue; end

        v = tdr.Araw_ord(idx,:);
        if any(~isfinite(v)) || norm(v) < opt.Eps, continue; end

        v = v(:)' / max(norm(v), opt.Eps);
        Vc{end+1,1} = v; %#ok<AGROW>
    end

    nAvail = numel(Vc);
    if nAvail == 0, continue; end

    % sign-align to reference
    vref = Vc{1};
    V = zeros(nAvail, K);
    for ii = 1:nAvail
        v = Vc{ii};
        if (v * vref') < 0
            v = -v;
        end
        V(ii,:) = v;
    end

    vbar = mean(V, 1);
    if norm(vbar) < opt.Eps, continue; end

    vbar = vbar / max(norm(vbar), opt.Eps);
    vbar = sign_fix_axis(vbar, opt.SignFix);

    A_pool(a,:) = vbar;

    [gName, pcIdx] = parse_axis_name(nm);
    meta(end+1).groupName = gName; %#ok<AGROW>
    meta(end).pcIdx       = pcIdx;
    meta(end).nSessAvail  = nAvail;
    meta(end).nSessUsed   = nAvail;
end

% keep valid pooled axes
keepAxis = all(isfinite(A_pool), 2) & (vecnorm(A_pool,2,2) > opt.Eps);
A_pool     = A_pool(keepAxis,:);
names_kept = namesAll(keepAxis);
meta       = meta(keepAxis);

if isempty(A_pool)
    error('buildAnchorAxes_from_glmRezList:EmptyPooledAxes', ...
        'No axes survived pooling. Check naming consistency / groups / degeneracy.');
end

% -------------------- order pooled axes: PriorityNames + pcIdx --------------------
key1 = nan(numel(names_kept),1); % group rank
key2 = nan(numel(names_kept),1); % pcIdx
key3 = (1:numel(names_kept))';   % stable tiebreak

for i = 1:numel(names_kept)
    [gName, pcIdx] = parse_axis_name(string(names_kept{i}));
    key2(i) = pcIdx;

    % robust: if not in PriorityNames -> huge rank
    hit = find(strcmpi(opt.PriorityNames, gName), 1, 'first');
    if isempty(hit)
        key1(i) = 1e6;
    else
        key1(i) = hit;
    end
end

[~, ord] = sortrows([key1 key2 key3], [1 2 3]);

Araw_ord  = A_pool(ord,:);
names_ord = names_kept(ord);
meta_ord  = meta(ord);

% -------------------- optional GS on pooled axes --------------------
if strcmpi(string(opt.OrthMode), "none")
    A = Araw_ord;
    keepAfterGS = true(1, size(Araw_ord,1));
    namesFinal  = names_ord;
    metaFinal   = meta_ord;
else
    [A, keepAfterGS] = gs_orth_rows(Araw_ord, opt.Eps);
    namesFinal = names_ord(keepAfterGS);
    metaFinal  = meta_ord(keepAfterGS);
end

% -------------------- pack output --------------------
axesOut = struct();
axesOut.Araw_rows     = A_pool;
axesOut.names_raw     = names_kept;

axesOut.Araw_ord      = Araw_ord;
axesOut.names_ord     = names_ord;
axesOut.axisMeta_ord  = meta_ord;

axesOut.A              = A;
axesOut.names          = namesFinal;
axesOut.axisMeta_final = metaFinal;
axesOut.keepAfterGS    = keepAfterGS;

axesOut.session        = sess;
axesOut.headers        = hdrGood;

axesOut.tdrBySession   = tdrBySession;
axesOut.sessionKeep    = sessionKeep;
axesOut.opt            = opt;

if opt.Verbose
    fprintf('[buildAnchorAxes] pooled axes=%d (kept after pooled GS=%d) from %d/%d sessions | per-session stored=%d | per-session GS=%d\n', ...
        size(Araw_ord,1), size(A,1), sum(sessionKeep), S, numel(sess), opt.PerSessionGS);
end

end



function tdr = targetedDimRed_from_glmRez(glmRez, varargin)
%TARGETEDDIMRED_FROM_GLMREZ  GLM-based targeted DR (2B) + optional GS + projection.
%
% tdr = targetedDimRed_from_glmRez(glmRez, 'Name', value, ...)
%
% CORE IDEA (2B / "effect-as-axis")
%   For each predictor group g (e.g., GoToneOn, Lick), define axes in motif-space
%   from the group-predicted component:
%       Yhat_g = Xz(:, cols_g) * beta(cols_g, :)
%       axis_g_pc = V(:,pc)' from SVD(Yhat_g)   (right-singular vectors; motif-space)
%   Then order axes by group priority + pcIdx, optionally apply Gram–Schmidt (GS),
%   and project motif activity into this low-D space.
%
% IMPORTANT STACKING (confirmed in your pipeline)
%   Rows are TIME-MAJOR:
%       [time bin 1: trials 1..N], [time bin 2: trials 1..N], ..., [time bin nW]
%
% REQUIRED glmRez fields
%   beta        : [P x K]
%   X_design    : [M x P]
%   muX         : [1 x P]
%   sdX         : [1 x P]
%   Yz          : [M x K]
%   group       : EITHER (1xG cell) of structs OR (1xG struct array), each with .name, .cols
%   decBins.time: vector used to infer nW
%
% OPTIONAL glmRez fields
%   Ybig        : [M x K]
%
% NAME-VALUE OPTIONS
%   'OrthMode'            : "GS" (default) or "none"
%   'SignFix'             : "maxabs" (default) or "none"
%   'Eps'                 : numeric (default 1e-10)
%   'ProjectWhichY'       : "Yz" (default) or "Ybig"
%   'MultiDimGroups'      : cellstr group names that should keep multiple PCs
%   'MultiDimK'           : # dims to keep for those groups (default 3)
%   'PriorityNames'       : ordering priority for group blocks (default tone blocks + Lick)
%
% OUTPUT tdr struct (key fields)
%   .Araw_rows     : [nAxis x K] raw axes (after sign-fix), unordered
%   .Araw_ord      : [nAxis x K] raw axes after ordering (pre-GS)
%   .A             : [G' x K] final axes after optional GS
%   .names_ord     : 1 x nAxis ordered axis names (pre-GS)
%   .names         : 1 x G' final axis names (post-GS)
%   .axisMeta      : struct array per axis (groupIdx, groupName, pcIdx, sv, expl)
%   .Xz            : [M x P] reconstructed standardized X
%   .Zrows         : [M x G'] projected rows
%   .Z             : [N x nW x G'] projected trajectories (trial x time x axis)
%   .nW, .N, .K, .P, .M
%
% Drop-in, robust: handles glmRez.group being cell-of-struct OR struct array.

% -------------------- parse --------------------
p = inputParser;
p.addParameter('OrthMode',"GS",@(s)ischar(s)||isstring(s));
p.addParameter('SignFix',"maxabs",@(s)ischar(s)||isstring(s));
p.addParameter('Eps',1e-10,@(x)isscalar(x) && x>0);
p.addParameter('ProjectWhichY',"Yz",@(s)ischar(s)||isstring(s));

p.addParameter('MultiDimGroups', {'GoToneOn','NoGoToneOn','ToneOffGo','ToneOffNoGo'}, @(c)iscell(c)||isstring(c));
p.addParameter('MultiDimK', 3, @(x)isnumeric(x)&&isscalar(x)&&x>=1);

p.addParameter('PriorityNames', {'GoToneOn','NoGoToneOn','ToneOffGo','ToneOffNoGo','Lick'}, @(c)iscell(c)||isstring(c));
p.parse(varargin{:});
opt = p.Results;

% normalize option lists
opt.MultiDimGroups = cellstr(string(opt.MultiDimGroups));
opt.PriorityNames  = cellstr(string(opt.PriorityNames));

% -------------------- validate glmRez --------------------
needFields = {'beta','X_design','muX','sdX','Yz','group','decBins'};
for f = needFields
    assert(isfield(glmRez,f{1}), 'glmRez missing required field: %s', f{1});
end
assert(isfield(glmRez.decBins,'time') && ~isempty(glmRez.decBins.time), ...
    'glmRez.decBins.time is required to infer nW.');

B    = glmRez.beta;        % [P x K]
Xraw = glmRez.X_design;    % [M x P]
muX  = glmRez.muX;         % [1 x P]
sdX  = glmRez.sdX;         % [1 x P]
Yz   = glmRez.Yz;          % [M x K]
grp  = glmRez.group;       % cell-of-struct OR struct array

[M, P]  = size(Xraw);
[Pb, K] = size(B);

assert(P==Pb, 'X_design has %d cols but beta has %d rows.', P, Pb);
assert(size(Yz,1)==M && size(Yz,2)==K, 'Yz must be [M x K] consistent with X_design/beta.');

% -------------------- reconstruct Xz (exact fitting space) --------------------
Xz = (Xraw - muX) ./ sdX;     % implicit expansion
Xz(~isfinite(Xz)) = 0;        % match fitting behavior

% -------------------- determine nW and N (TIME-MAJOR) --------------------
nW = numel(glmRez.decBins.time);
remVal = rem(M, nW);
if remVal ~= 0
    error('M=%d not divisible by nW=%d (remainder=%d). Check stacking/metadata.', M, nW, remVal);
end
N = M / nW;

% -------------------- group names + count --------------------
G = local_num_groups_(grp);
groupNames = cell(1,G);
for g = 1:G
    gg = local_get_group_(grp, g);
    groupNames{g} = char(string(gg.name));
end

% -------------------- build raw axes (2B: SVD of group-predicted component) --------------------
Araw_rows = [];     % [nAxis x K]
axisNames = {};     % 1 x nAxis
axisMeta  = struct('groupIdx',{},'groupName',{},'pcIdx',{},'sv',{},'expl',{});

multiNames = opt.MultiDimGroups;
multiK     = opt.MultiDimK;

for g = 1:G
    gg    = local_get_group_(grp, g);
    gName = string(gg.name);
    cols  = gg.cols(:)';

    % sanity + bound
    cols = cols(isfinite(cols));
    cols = cols(cols>=1 & cols<=P);
    cols = unique(round(cols), 'stable');
    if isempty(cols), continue; end

    Yhat_g = Xz(:,cols) * B(cols,:);   % [M x K]
    if norm(Yhat_g,'fro') < opt.Eps, continue; end

    % dims to keep for this group
    if any(strcmpi(char(gName), multiNames))
        nKeep = multiK;
    else
        nKeep = 1;
    end

    [~,Ssvd,V] = svd(Yhat_g, 'econ');  % V: [K x r]
    r = size(V,2);
    nKeep = min(nKeep, r);

    s = diag(Ssvd);
    denom = max(sum(s.^2), opt.Eps);

    for pc = 1:nKeep
        v = V(:,pc)';                  % 1 x K
        if norm(v) < opt.Eps, continue; end

        v = v / max(norm(v), opt.Eps);
        v = sign_fix_axis_(v, opt.SignFix);

        Araw_rows(end+1,:) = v; %#ok<AGROW>
        axisNames{end+1}   = char(gName + "_" + string(pc)); %#ok<AGROW>

        axisMeta(end+1).groupIdx  = g; %#ok<AGROW>
        axisMeta(end).groupName   = char(gName);
        axisMeta(end).pcIdx       = pc;
        axisMeta(end).sv          = s(pc);
        axisMeta(end).expl        = (s(pc)^2) / denom;
    end
end

if isempty(Araw_rows)
    error('No valid axes were constructed (Araw_rows is empty). Check groups/cols/Xz/beta.');
end

% -------------------- axis ordering: PriorityNames + pcIdx --------------------
% Desired: GoToneOn_1,2,3 -> NoGoToneOn_1,2,3 -> ToneOffGo_1,2,3 -> ToneOffNoGo_1,2,3 -> Lick_1 -> ...
rankPerGroup = 1:G;  % default: natural group order
for i = 1:numel(opt.PriorityNames)
    gIdx = find(strcmpi(groupNames, opt.PriorityNames{i}), 1, 'first');
    if ~isempty(gIdx)
        rankPerGroup(gIdx) = -1000 + i; % force priority groups to front
    end
end

nAxis = size(Araw_rows,1);
key1 = nan(nAxis,1);  % group rank
key2 = nan(nAxis,1);  % pcIdx
for a = 1:nAxis
    key1(a) = rankPerGroup(axisMeta(a).groupIdx);
    key2(a) = axisMeta(a).pcIdx;
end

[~, orderAxes] = sortrows([key1 key2], [1 2]);

Araw_ord     = Araw_rows(orderAxes,:);
names_ord    = axisNames(orderAxes);
axisMeta_ord = axisMeta(orderAxes);

% -------------------- orthogonalize (optional) --------------------
if strcmpi(string(opt.OrthMode),"none")
    A            = Araw_ord;
    keepLocal    = true(1, size(Araw_ord,1));
    namesFinal   = names_ord;
    axisMetaFinal= axisMeta_ord;
else
    [A, keepLocal] = gs_orth_rows_(Araw_ord, opt.Eps);
    namesFinal    = names_ord(keepLocal);
    axisMetaFinal = axisMeta_ord(keepLocal);
end

% -------------------- choose Y for projection --------------------
switch lower(string(opt.ProjectWhichY))
    case "yz"
        Yproj = Yz;
    case "ybig"
        assert(isfield(glmRez,'Ybig') && ~isempty(glmRez.Ybig), ...
            'ProjectWhichY="Ybig" requested but glmRez.Ybig not found/empty.');
        Yproj = glmRez.Ybig;
        assert(all(size(Yproj)==[M K]), 'glmRez.Ybig must be [M x K].');
    otherwise
        error('Unknown ProjectWhichY: %s', string(opt.ProjectWhichY));
end

% -------------------- project into targeted space --------------------
Zrows  = Yproj * A';                  % [M x G']
Gprime = size(A,1);
Z      = reshape(Zrows, [N, nW, Gprime]);  % TIME-MAJOR invert

% -------------------- pack outputs --------------------
tdr = struct();
tdr.M = M; tdr.P = P; tdr.K = K; tdr.nW = nW; tdr.N = N;

tdr.Xz = Xz;

tdr.Araw_rows = Araw_rows;
tdr.axisNames = axisNames;
tdr.axisMeta  = axisMeta;

tdr.Araw_ord     = Araw_ord;
tdr.names_ord    = names_ord;
tdr.axisMeta_ord = axisMeta_ord;

tdr.A              = A;
tdr.names          = namesFinal;
tdr.axisMeta_final = axisMetaFinal;

tdr.Zrows = Zrows;
tdr.Z     = Z;

tdr.orderAxes    = orderAxes;
tdr.keepAfterGS  = keepLocal;

tdr.opt = opt;

end

%% ========================== local helpers ==========================
function G = local_num_groups_(grp)
if iscell(grp)
    G = numel(grp);
elseif isstruct(grp)
    G = numel(grp);
else
    error('glmRez.group must be a cell array or struct array, got: %s', class(grp));
end
end

function gg = local_get_group_(grp, g)
% supports struct array OR cell array of structs
if iscell(grp)
    gg = grp{g};
else
    gg = grp(g);
end
assert(isstruct(gg) && isfield(gg,'name') && isfield(gg,'cols'), ...
    'glmRez.group(%d) must contain a struct with fields .name and .cols.', g);
end


% ========================== local helpers ==========================
function v = sign_fix_axis(v, mode)
mode = lower(string(mode));
switch mode
    case "maxabs"
        [~,idx] = max(abs(v));
        if v(idx) < 0, v = -v; end
    otherwise
        % no-op
end
end

function v = sign_fix_axis_(v, mode)
% Backward-compat shim: older code calls sign_fix_axis_ but newer code uses sign_fix_axis
v = sign_fix_axis(v, mode);
end


function [gName, pcIdx] = parse_axis_name(nm)
nm = char(string(nm));
tok = regexp(nm, '^(.*)_(\d+)$', 'tokens', 'once');
if isempty(tok)
    gName = nm;
    pcIdx = 1;
else
    gName = tok{1};
    pcIdx = str2double(tok{2});
    if ~isfinite(pcIdx) || pcIdx < 1, pcIdx = 1; end
end
end

function [Q, keepIdx] = gs_orth_rows(A, epsVal)
%GS_ORTH_ROWS  Gram–Schmidt on row vectors; returns orthonormal rows.
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
