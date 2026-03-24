function anchor = buildPerMouseLOOAnchors_fromExperts(headerC, glmRezC, trIdC, expertHeaders_perMouse, varargin)
%buildPerMouseLOOAnchors_fromExperts  Per-mouse pooled GLM-TDR anchors with LOO/LOSO.
%
% LEAN VERSION:
%   By default, returns only the essential fields needed for subsequent
%   projection / downstream analysis:
%       - A
%       - names
%   plus a few bookkeeping fields at the mouse level.
%
% NEW NAME-VALUE
%   'LeanOutput'      : true (default)
%       If true, each axes struct is slimmed to:
%           .A
%           .names
%
%   'KeepOpt'         : true (default)
%       Keep anchor.opt in output.
%
%   'KeepFoundMat'    : true (default)
%       Keep foundMat in output.
%
%   'KeepHeadersUsed' : true (default)
%       Keep headersUsed / headersUsedLOO in output.
%
% NOTES
%   - This does NOT change the math or the constructed anchors.
%   - It only trims the returned struct to avoid huge memory / save overhead.

% -------------------- parse --------------------
p = inputParser;
p.addRequired('headerC', @(c) iscell(c));
p.addRequired('glmRezC', @(c) iscell(c) && isequal(size(c), size(headerC)));
p.addRequired('trIdC',   @(c) iscell(c) && isequal(size(c), size(headerC)));
p.addRequired('expertHeaders_perMouse', @(c) iscell(c) || isstring(c));

p.addParameter('OrthMode',"GS",@(s)ischar(s)||isstring(s));
p.addParameter('Eps',1e-10,@(x)isnumeric(x)&&isscalar(x)&&x>0);

p.addParameter('Verbose', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('VerboseEvery', 1, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('VerboseGroups', false, @(x)islogical(x)&&isscalar(x));

p.addParameter('SignFix',"projmean",@(s)ischar(s)||isstring(s));
p.addParameter('SignFixWinSec',[0 4],@(x)isnumeric(x)&&numel(x)==2);
p.addParameter('SignFixTrialPolicy',"GoNoGoByGroup",@(s)ischar(s)||isstring(s));

p.addParameter('MultiDimGroups', {'GoToneOn','NoGoToneOn','ToneOffGo','ToneOffNoGo'}, @(c)iscell(c)||isstring(c));
p.addParameter('MultiDimK', 3, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('PriorityNames', {'GoToneOn','NoGoToneOn','ToneOffGo','ToneOffNoGo','Lick'}, @(c)iscell(c)||isstring(c));

p.addParameter('AlignMode',"procrustes_within_group",@(s)ischar(s)||isstring(s));
p.addParameter('AlignRef',"latest_session_only",@(s)ischar(s)||isstring(s));

p.addParameter('SessionWeighting', "fro", @(s)ischar(s)||isstring(s));

% NEW
p.addParameter('LeanOutput', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('KeepOpt', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('KeepFoundMat', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('KeepHeadersUsed', true, @(x)islogical(x)&&isscalar(x));

p.parse(headerC, glmRezC, trIdC, expertHeaders_perMouse, varargin{:});
opt = p.Results;

opt.MultiDimGroups = cellstr(string(opt.MultiDimGroups));
opt.PriorityNames  = cellstr(string(opt.PriorityNames));
opt.SessionWeighting = lower(string(opt.SessionWeighting));
opt.OrthMode = string(opt.OrthMode);
opt.AlignMode = lower(string(opt.AlignMode));
opt.AlignRef  = lower(string(opt.AlignRef));

% -------------------- normalize expert header matrix --------------------
hdrPM = normalizeHeaderMatrix_(expertHeaders_perMouse);
[nM, nK] = size(hdrPM);

% -------------------- build header lookup from headerC -> glmRezC/trIdC --------------------
[uniqHdr, glmFirst, trFirst] = buildHeaderLookup_(headerC, glmRezC, trIdC);

% -------------------- prealloc outputs --------------------
foundMat         = false(nM, nK);
mouseIds         = cell(nM,1);
axesFull         = cell(nM,1);
axesLOO          = cell(nM,nK);
headersUsed      = cell(nM,1);
headersUsedLOO   = cell(nM,nK);

% ============================================================
% main loop over mice (rows)
% ============================================================
t0 = tic;

for m = 1:nM
    row = hdrPM(m,:);
    mouseIds{m} = inferMouseIdFromRow_(row);

    % Gather found expert sessions for this mouse
    hdrFound = strings(0,1);
    glmList  = {};
    trList   = {};
    kIdxList = [];

    for k = 1:nK
        h = row{k};
        if isempty(h), continue; end
        hs = string(h);

        [gr, tr, ok] = fetchByHeaderFromLookup_(hs, uniqHdr, glmFirst, trFirst);
        foundMat(m,k) = ok;

        if ok
            hdrFound(end+1,1) = hs; %#ok<AGROW>
            glmList{end+1,1}  = gr; %#ok<AGROW>
            trList{end+1,1}   = tr; %#ok<AGROW>
            kIdxList(end+1,1) = k;  %#ok<AGROW>
        end
    end

    headersUsed{m} = cellstr(hdrFound);

    if opt.Verbose && (mod(m, opt.VerboseEvery)==0 || m==1 || m==nM)
        fprintf('[LOOAnchors] mouse %d/%d (%s): found %d/%d expert sessions\n', ...
            m, nM, string(mouseIds{m}), numel(glmList), nK);
    end

    % --- build axesFull (pool all found experts) ---
    if numel(glmList) >= 1
        axFull = buildPooledAxes_concatYhatSVD_(glmList, trList, cellstr(hdrFound), opt);

        refHdr_full = pickLatestHeader_(hdrFound);
        refAxes_full = buildReferenceAxes_fromHeader_(refHdr_full, uniqHdr, glmFirst, trFirst, opt);

        axFull = alignAxesToRef_(axFull, refAxes_full, opt, char(refHdr_full));

        if opt.LeanOutput
            axFull = slimAxesStruct_(axFull);
        end
        axesFull{m} = axFull;
    else
        axesFull{m} = [];
    end

    % --- build axesLOO{k}: exclude the k-th expert (if that k was found) ---
    if numel(glmList) >= 2
        for kk = 1:nK
            if ~foundMat(m,kk)
                axesLOO{m,kk} = [];
                headersUsedLOO{m,kk} = {};
                continue;
            end

            keep = (kIdxList ~= kk);
            if sum(keep) < 1
                axesLOO{m,kk} = [];
                headersUsedLOO{m,kk} = {};
                continue;
            end

            glmSub = glmList(keep);
            trSub  = trList(keep);
            hdrSub = hdrFound(keep);

            headersUsedLOO{m,kk} = cellstr(hdrSub);

            pooled = buildPooledAxes_concatYhatSVD_(glmSub, trSub, cellstr(hdrSub), opt);

            refHdr_loo = pickLatestHeader_(hdrSub);
            refAxes_loo = buildReferenceAxes_fromHeader_(refHdr_loo, uniqHdr, glmFirst, trFirst, opt);

            axLoo = alignAxesToRef_(pooled, refAxes_loo, opt, char(refHdr_loo));

            if opt.LeanOutput
                axLoo = slimAxesStruct_(axLoo);
            end
            axesLOO{m,kk} = axLoo;

            if opt.Verbose && opt.VerboseGroups && (mod(m,opt.VerboseEvery)==0 || m==1 || m==nM)
                fprintf('  - LOO k=%d: pooled %d sessions | ref=%s\n', kk, sum(keep), char(refHdr_loo));
            end
        end
    else
        for kk = 1:nK
            axesLOO{m,kk} = [];
            headersUsedLOO{m,kk} = {};
        end
    end
end

% -------------------- pack output --------------------
anchor = struct();
anchor.perMouse = struct();

anchor.perMouse.headersMat = hdrPM;
anchor.perMouse.mouseIds   = mouseIds;
anchor.perMouse.axesFull   = axesFull;
anchor.perMouse.axesLOO    = axesLOO;

if opt.KeepFoundMat
    anchor.perMouse.foundMat = foundMat;
end

if opt.KeepHeadersUsed
    anchor.perMouse.headersUsed    = headersUsed;
    anchor.perMouse.headersUsedLOO = headersUsedLOO;
end

if opt.KeepOpt
    anchor.opt = opt;
end

if opt.Verbose
    fprintf('[LOOAnchors] done: mice=%d | elapsed=%.1fs\n', nM, toc(t0));
end
end

%% ======================================================================
% Core pooled builder (internal)
% ======================================================================

function axesOut = buildPooledAxes_concatYhatSVD_(glmRezList, trIdList, headers, opt)
% Pool across sessions in glmRezList (cell), per group:
%   - compute Yhat_g = Xz(:,cols)*beta(cols,:)
%   - optional session weighting
%   - vertical concat across sessions -> SVD once -> V(:,pc) gives motif-space axis
%   - sign-fix (projmean across sessions) or maxabs/none
%   - order by PriorityNames then pcIdx
%   - optional GS (OrthMode)

glmRezList = glmRezList(:);
trIdList   = trIdList(:);
S = numel(glmRezList);

grp0 = glmRezList{1}.group;
G = local_num_groups_(grp0);

groupNames = strings(1,G);
groupCols  = cell(1,G);
for g = 1:G
    gg = local_get_group_(grp0, g);
    groupNames(g) = string(gg.name);
    groupCols{g}  = unique(round(gg.cols(:)'),'stable');
end

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
        assert(isfield(gr,'X_design') && isfield(gr,'beta') && isfield(gr,'muX') && isfield(gr,'sdX'), ...
            'glmRez missing required GLM fields');

        P = size(gr.X_design,2);
        cols = cols0(cols0>=1 & cols0<=P);
        if isempty(cols), continue; end

        % memory-friendlier than standardizing full X_design
        Xsub = gr.X_design(:,cols);
        muSub = gr.muX(cols);
        sdSub = gr.sdX(cols);

        Xzsub = (Xsub - muSub) ./ sdSub;
        Xzsub(~isfinite(Xzsub)) = 0;

        Yhat = Xzsub * gr.beta(cols,:); % [M x K]

        if opt.SessionWeighting == "fro"
            Yhat = Yhat / max(norm(Yhat,'fro'), opt.Eps);
        elseif opt.SessionWeighting == "none"
            % no-op
        else
            error('Unknown SessionWeighting: %s', opt.SessionWeighting);
        end

        if norm(Yhat,'fro') >= opt.Eps
            Ypool = [Ypool; Yhat]; %#ok<AGROW>
        end
    end

    if isempty(Ypool) || norm(Ypool,'fro') < opt.Eps
        continue;
    end

    [~, Ssvd, V] = svd(Ypool, 'econ');
    r = size(V,2);
    nKeepEff = min(nKeep, r);

    svals = diag(Ssvd);
    denom = max(sum(svals.^2), opt.Eps);

    for pc = 1:nKeepEff
        v = V(:,pc)';
        v = v / max(norm(v), opt.Eps);

        if strcmpi(string(opt.SignFix),"projmean")
            sc = nan(S,1);
            for ss = 1:S
                sc(ss) = local_projmean_score_glm_(glmRezList{ss}, trIdList{ss}, v, gName, opt.SignFixWinSec, opt.SignFixTrialPolicy, opt.Eps);
            end
            if mean(sc,'omitnan') < 0, v = -v; end
        elseif strcmpi(string(opt.SignFix),"maxabs")
            v = sign_fix_axis_(v, "maxabs");
        elseif strcmpi(string(opt.SignFix),"none")
            % no-op
        else
            error('Unknown SignFix: %s', string(opt.SignFix));
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

if strcmpi(string(opt.OrthMode),"none")
    A = Araw_ord;
    keep = true(1,size(Araw_ord,1));
    namesFinal = names_ord;
    metaFinal  = axisMeta_ord;
elseif strcmpi(string(opt.OrthMode),"gs")
    [A, keep] = gs_orth_rows_(Araw_ord, opt.Eps);
    namesFinal = names_ord(keep);
    metaFinal  = axisMeta_ord(keep);
else
    error('Unknown OrthMode: %s', string(opt.OrthMode));
end

axesOut = struct();
axesOut.Araw_ord       = Araw_ord;
axesOut.names_ord      = names_ord;
axesOut.axisMeta_ord   = axisMeta_ord;

axesOut.A              = A;
axesOut.keepAfterGS    = keep;
axesOut.names          = namesFinal;
axesOut.axisMeta_final = metaFinal;

axesOut.headers        = headers;
axesOut.opt            = opt;
end

%% ======================================================================
% Small helpers (standalone)
% ======================================================================

function ax = slimAxesStruct_(ax)
% Keep only the fields required for downstream projection.
if isempty(ax)
    return;
end
ax = struct( ...
    'A', ax.A, ...
    'names', {ax.names});
end

function hdrMat = normalizeHeaderMatrix_(hdrIn)
if isempty(hdrIn), hdrMat = cell(0,0); return; end
if isstring(hdrIn), hdrIn = cellstr(hdrIn); end
if iscell(hdrIn) && isvector(hdrIn), hdrIn = reshape(hdrIn, [], 1); end

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
    warning('buildPerMouseLOOAnchors:DuplicateHeaders', ...
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

function score = local_projmean_score_glm_(glmRez, trId, axisRow, groupName, winSec, policy, epsVal)
assert(isfield(glmRez,'Yz') && isfield(glmRez,'decBins') && isfield(glmRez.decBins,'time'), ...
    'glmRez missing Yz/decBins.time');
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
end
end

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

function v = sign_fix_axis_(v, mode)
mode = lower(string(mode));
switch mode
    case "maxabs"
        [~,idx] = max(abs(v));
        if v(idx) < 0, v = -v; end
    otherwise
end
end

function [Q, keepIdx] = gs_orth_rows_(A, epsVal)
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

function hLatest = pickLatestHeader_(hdrS)
hdrS = string(hdrS(:));
if isempty(hdrS)
    hLatest = "";
    return;
end

[dtKey, ok] = parseHeaderDatetimeWithSuffix_local_(hdrS);
if ok
    [~, ord] = sort(dtKey, 'ascend');
    hLatest = hdrS(ord(end));
else
    hLatest = hdrS(end);
end
end

function [dtKey, ok] = parseHeaderDatetimeWithSuffix_local_(hdrS)
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
    if numel(tok) >= 2 && ~isempty(tok{2})
        suf = str2double(tok{2});
        if ~isfinite(suf), suf = 0; end
    else
        suf = 0;
    end
    d0 = datetime(mmddyy, 'InputFormat','MMddyy');
    dtKey(i) = d0 + seconds(suf);
end
end

function refAxes = buildReferenceAxes_fromHeader_(refHdr, uniqHdr, glmFirst, trFirst, opt)
refHdr = string(refHdr);
if strlength(refHdr)==0
    refAxes = [];
    return;
end

[gr, tr, ok] = fetchByHeaderFromLookup_(refHdr, uniqHdr, glmFirst, trFirst);
if ~ok
    refAxes = [];
    return;
end

refAxes = buildPooledAxes_concatYhatSVD_({gr}, {tr}, {char(refHdr)}, opt);
end

function axesOut = alignAxesToRef_(axesOut, refAxes, opt, refHdrStr)
if isempty(axesOut) || isempty(refAxes) || opt.AlignMode=="none"
    return;
end

A  = axesOut.Araw_ord;
Ar = refAxes.Araw_ord;

meta  = axesOut.axisMeta_ord;
metar = refAxes.axisMeta_ord;

gNames = unique(string({meta.groupName}));
gNamesR = unique(string({metar.groupName}));
gNames = intersect(gNames, gNamesR, 'stable');

Anew = A;

for gi = 1:numel(gNames)
    gName = gNames(gi);

    idxT = find(strcmpi(string({meta.groupName}), gName));
    idxR = find(strcmpi(string({metar.groupName}), gName));

    if isempty(idxT) || isempty(idxR)
        continue;
    end

    pcT = [meta(idxT).pcIdx];
    pcR = [metar(idxR).pcIdx];

    commonPC = intersect(pcT, pcR, 'stable');
    if isempty(commonPC)
        continue;
    end

    idxT2 = idxT(ismember(pcT, commonPC));
    idxR2 = idxR(ismember(pcR, commonPC));

    [~, oT] = sort([meta(idxT2).pcIdx]);
    [~, oR] = sort([metar(idxR2).pcIdx]);

    idxT2 = idxT2(oT);
    idxR2 = idxR2(oR);

    At = Anew(idxT2,:);
    Br = Ar(idxR2,:);
    k  = size(At,1);

    if k==1
        if (At*Br') < 0
            At = -At;
        end
        Anew(idxT2,:) = At;
        continue;
    end

    switch opt.AlignMode
        case "match_sign_perm"
            Anew(idxT2,:) = match_sign_perm_rows_(At, Br);
        case "procrustes_within_group"
            Anew(idxT2,:) = procrustes_within_group_rows_(At, Br);
        otherwise
            error('Unknown AlignMode: %s', opt.AlignMode);
    end
end

axesOut.Araw_ord = Anew;

if strcmpi(string(opt.OrthMode),"none")
    axesOut.A = axesOut.Araw_ord;
    axesOut.keepAfterGS = true(1,size(axesOut.Araw_ord,1));
    axesOut.names = axesOut.names_ord;
    axesOut.axisMeta_final = axesOut.axisMeta_ord;
else
    [A_gs, keep_gs] = gs_orth_rows_(axesOut.Araw_ord, opt.Eps);
    axesOut.A = A_gs;
    axesOut.keepAfterGS = keep_gs;
    axesOut.names = axesOut.names_ord(keep_gs);
    axesOut.axisMeta_final = axesOut.axisMeta_ord(keep_gs);
end
end

function At_aligned = procrustes_within_group_rows_(At, Br)
M = Br * At';
[U,~,V] = svd(M, 'econ');
R = U * V';
At_aligned = R * At;
end

function At_best = match_sign_perm_rows_(At, Br)
k = size(At,1);
P = perms(1:k);
bestScore = -Inf;
At_best = At;

for i = 1:size(P,1)
    Aperm = At(P(i,:),:);
    sgn = sign(sum(Aperm .* Br, 2));
    sgn(sgn==0) = 1;
    Aperm = Aperm .* sgn;

    sc = sum(sum(Aperm .* Br));
    if sc > bestScore
        bestScore = sc;
        At_best = Aperm;
    end
end
end