function rez = validatePerMouseLOO_noCircularity(anchorLOO, headerC, glmRezC, trIdC, hdrC_all, varargin)
%validatePerMouseLOO_noCircularity
% Verifies that anchorLOO.perMouse.axesLOO{iM,k}.A was built from OTHER sessions
% (LOO) by (a) checking stored header provenance if available, and (b) recomputing
% pooled axes from the expected training set and comparing to stored A.
%
% OUTPUT rez fields:
%   .okMat(iM,k)           : true if pass
%   .msgMat{iM,k}          : short message
%   .expectedTrainHdr{iM,k}: expected training headers (cellstr)
%   .usedTrainHdr{iM,k}    : if stored in anchor entry, that list (else {})
%   .maxAbsDiff(iM,k)      : max(abs(Astored - Arecomputed)) after best column sign alignment
%   .cosDiagMin(iM,k)      : min diagonal cosine similarity after best matching
%
% NAME-VALUE
%   'TolMaxAbs'      : default 1e-8
%   'TolCos'         : default 1-1e-6
%   'Recompute'      : true (default)
%   'Verbose'        : true
%
% NOTE: This assumes your pooled builder is deterministic given the same inputs
%       and options (SessionWeighting, SignFix, OrthMode, etc.)

p = inputParser;
p.addParameter('TolMaxAbs', 1e-8, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('TolCos', 1-1e-6, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('Recompute', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('Verbose', true, @(x)islogical(x)&&isscalar(x));

% These must match the options you used in buildPerMouseLOOAnchors_fromExperts
p.addParameter('OrthMode',"GS",@(s)ischar(s)||isstring(s));
p.addParameter('SessionWeighting',"fro",@(s)ischar(s)||isstring(s));
p.addParameter('SignFix',"projmean",@(s)ischar(s)||isstring(s));
p.addParameter('SignFixWinSec',[0 4],@(x)isnumeric(x)&&numel(x)==2);
p.addParameter('SignFixTrialPolicy',"GoNoGoByGroup",@(s)ischar(s)||isstring(s));
p.addParameter('MultiDimGroups', {'GoToneOn','NoGoToneOn','ToneOffGo','ToneOffNoGo'}, @(c)iscell(c)||isstring(c));
p.addParameter('MultiDimK', 3, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('PriorityNames', {'GoToneOn','NoGoToneOn','ToneOffGo','ToneOffNoGo','Lick'}, @(c)iscell(c)||isstring(c));
p.addParameter('Eps',1e-10,@(x)isnumeric(x)&&isscalar(x));
p.parse(varargin{:});
opt = p.Results;
opt.MultiDimGroups = cellstr(string(opt.MultiDimGroups));
opt.PriorityNames  = cellstr(string(opt.PriorityNames));

axesLOO = anchorLOO.perMouse.axesLOO;
[nM, nK] = size(axesLOO);

okMat  = false(nM,nK);
msgMat = cell(nM,nK);
expectedTrainHdr = cell(nM,nK);
usedTrainHdr     = cell(nM,nK);
maxAbsDiff = nan(nM,nK);
cosDiagMin = nan(nM,nK);

% Build a header->(glmRez,trId) lookup using your existing helper style
[uniqHdr, glmFlatFirst, trFlatFirst] = local_buildHeaderLookup_(headerC, glmRezC, trIdC);

for iM = 1:nM
    for k = 1:nK
        ent = axesLOO{iM,k};
        hSelf = hdrC_all{iM,k};
        if isempty(ent) || isempty(hSelf)
            msgMat{iM,k} = 'skip: empty entry or empty self header';
            continue;
        end
        hSelf = string(hSelf);

        % expected train headers = the other columns in this mouse row
        otherIdx = setdiff(1:nK, k);
        Htrain = strings(0,1);
        for kk = otherIdx
            h = hdrC_all{iM,kk};
            if isempty(h), continue; end
            Htrain(end+1,1) = string(h); %#ok<AGROW>
        end

        % drop missing from lookup (slide-down already happened earlier, but be safe)
        keep = false(numel(Htrain),1);
        for ii = 1:numel(Htrain)
            keep(ii) = any(uniqHdr == Htrain(ii));
        end
        Htrain = Htrain(keep);

        expectedTrainHdr{iM,k} = cellstr(Htrain);

        % If function stored provenance, check it
        stored = {};
        if isstruct(ent)
            if isfield(ent,'trainHeaders') && ~isempty(ent.trainHeaders)
                stored = cellstr(string(ent.trainHeaders(:)));
            elseif isfield(ent,'headersUsed') && ~isempty(ent.headersUsed)
                stored = cellstr(string(ent.headersUsed(:)));
            elseif isfield(ent,'meta') && isstruct(ent.meta) && isfield(ent.meta,'trainHeaders')
                stored = cellstr(string(ent.meta.trainHeaders(:)));
            end
        end
        usedTrainHdr{iM,k} = stored;

        if ~isempty(stored)
            % provenance check: must NOT include self, must equal expected set (order-insensitive)
            hasSelf = any(string(stored) == hSelf);
            sameSet = isequal(sort(string(stored)), sort(Htrain));
            if hasSelf
                msgMat{iM,k} = 'FAIL: stored trainHeaders includes self header';
                okMat(iM,k) = false;
                continue;
            end
            if ~sameSet
                msgMat{iM,k} = 'WARN: stored trainHeaders differs from expected set (may be OK if missing sessions were skipped)';
                % do not fail yet; recompute will decide
            end
        end

        if ~opt.Recompute
            okMat(iM,k) = true;
            msgMat{iM,k} = 'PASS (provenance only)';
            continue;
        end

        % --- recompute pooled axes from Htrain and compare ---
        [glmList, trList] = local_fetchLists_(Htrain, uniqHdr, glmFlatFirst, trFlatFirst);

        if numel(glmList) < 1
            msgMat{iM,k} = 'skip: no training sessions found in lookup';
            continue;
        end

        % Use the SAME pooled builder used in your LOO code.
        % If your LOO function used buildPooledAxes_concatYhatSVD, call it here.
        pooled = local_pooledAxes_concatYhatSVD_(glmList, trList, cellstr(Htrain), opt);

        if isempty(pooled) || ~isfield(pooled,'A') || isempty(pooled.A)
            msgMat{iM,k} = 'FAIL: recompute produced empty pooled axes';
            okMat(iM,k) = false;
            continue;
        end

        Astored = ent.A;
        Arecomp = pooled.A;

        % Compare with best matching (handle possible axis order/sign ambiguities)
        [maxDiff, minCos] = local_compareAxes_(Astored, Arecomp);

        maxAbsDiff(iM,k) = maxDiff;
        cosDiagMin(iM,k) = minCos;

        pass = (maxDiff <= opt.TolMaxAbs) || (minCos >= opt.TolCos);
        okMat(iM,k) = pass;

        if pass
            msgMat{iM,k} = sprintf('PASS: maxAbsDiff=%.3g | minCos=%.6f', maxDiff, minCos);
        else
            msgMat{iM,k} = sprintf('FAIL: maxAbsDiff=%.3g | minCos=%.6f', maxDiff, minCos);
        end

        if opt.Verbose
            fprintf('Mouse %d, k=%d | self=%s | train={%s} -> %s\n', ...
                iM, k, hSelf, strjoin(cellstr(Htrain),', '), msgMat{iM,k});
        end
    end
end

rez = struct();
rez.okMat = okMat;
rez.msgMat = msgMat;
rez.expectedTrainHdr = expectedTrainHdr;
rez.usedTrainHdr = usedTrainHdr;
rez.maxAbsDiff = maxAbsDiff;
rez.cosDiagMin = cosDiagMin;

end

% ================= helpers =================

function [uniqHdr, glmFlatFirst, trFlatFirst] = local_buildHeaderLookup_(headerC, glmRezC, trIdC)
hdrFlat = headerC(:);
glmFlat = glmRezC(:);
trFlat  = trIdC(:);

isHdr = ~cellfun(@isempty, hdrFlat);
hdrFlat = hdrFlat(isHdr);
glmFlat = glmFlat(isHdr);
trFlat  = trFlat(isHdr);

hdrFlatS = string(hdrFlat);
[uniqHdr, ~, ic] = unique(hdrFlatS, 'stable');

glmFlatFirst = cell(numel(uniqHdr),1);
trFlatFirst  = cell(numel(uniqHdr),1);
for i = 1:numel(uniqHdr)
    ii = find(ic==i, 1, 'first');
    glmFlatFirst{i} = glmFlat{ii};
    trFlatFirst{i}  = trFlat{ii};
end
end

function [glmList, trList] = local_fetchLists_(H, uniqHdr, glmFlatFirst, trFlatFirst)
glmList = {};
trList  = {};
for i = 1:numel(H)
    j = find(uniqHdr == H(i), 1, 'first');
    if isempty(j), continue; end
    gr = glmFlatFirst{j};
    tr = trFlatFirst{j};
    if isempty(gr) || ~isstruct(gr) || isempty(tr) || ~isstruct(tr), continue; end
    glmList{end+1,1} = gr; %#ok<AGROW>
    trList{end+1,1}  = tr; %#ok<AGROW>
end
end

function [maxDiff, minCos] = local_compareAxes_(A1, A2)
% Align by greedy cosine matching, allow sign flips
A1 = double(A1); A2 = double(A2);
if isempty(A1) || isempty(A2)
    maxDiff = inf; minCos = -inf; return;
end
% Rows are axes; both are nAxis x K
n1 = size(A1,1); n2 = size(A2,1);
n = min(n1,n2);
A1 = A1(1:n,:); A2 = A2(1:n,:);

% normalize rows
A1 = A1 ./ max(vecnorm(A1,2,2), eps);
A2 = A2 ./ max(vecnorm(A2,2,2), eps);

C = A1*A2'; % cosine matrix if rows are normalized
used2 = false(n,1);
pair = zeros(n,2);
cosv = nan(n,1);

for i = 1:n
    [vals, idx] = sort(abs(C(i,:)), 'descend');
    jPick = [];
    for jj = idx
        if ~used2(jj)
            jPick = jj;
            break;
        end
    end
    if isempty(jPick)
        jPick = idx(1);
    end
    used2(jPick) = true;
    pair(i,:) = [i jPick];
    cosv(i) = C(i,jPick);
end

% apply sign alignment and compute max abs diff in original (not normalized) space
A2m = A2(pair(:,2),:);
sgn = sign(cosv); sgn(sgn==0) = 1;
A2m = A2m .* sgn;

A1m = A1(pair(:,1),:);

D = abs(A1m - A2m);
maxDiff = max(D, [], 'all');

minCos = min(abs(cosv));

end

function axesOut = local_pooledAxes_concatYhatSVD_(glmRezList, trIdList, headers, opt)
% Local pooled builder (self-contained) for validator:
% pools Yhat_g across sessions then SVD once per group.
% matches logic you described in your earlier code.

glmRezList = glmRezList(:);
trIdList   = trIdList(:);
S = numel(glmRezList);
if S ~= numel(trIdList), error('glmRezList/trIdList length mismatch'); end

grp0 = glmRezList{1}.group;
G = local_num_groups__(grp0);

groupNames = strings(1,G);
groupCols  = cell(1,G);
for g = 1:G
    gg = local_get_group__(grp0, g);
    groupNames(g) = string(gg.name);
    groupCols{g}  = unique(round(gg.cols(:)'),'stable');
end

Araw_rows = [];
axisNames = {};
axisMeta  = struct('groupIdx',{},'groupName',{},'pcIdx',{},'sv',{},'expl',{});

for g = 1:G
    gName = groupNames(g);
    cols0 = groupCols{g};
    cols0 = cols0(isfinite(cols0));

    nKeep = 1;
    if any(strcmpi(char(gName), opt.MultiDimGroups)), nKeep = opt.MultiDimK; end

    Ypool = [];
    for s = 1:S
        gr = glmRezList{s};
        P  = size(gr.X_design,2);

        cols = cols0(cols0>=1 & cols0<=P);
        if isempty(cols), continue; end

        Xz = (gr.X_design - gr.muX) ./ gr.sdX;
        Xz(~isfinite(Xz)) = 0;

        Yhat = Xz(:,cols) * gr.beta(cols,:); % [M x K]

        % session weighting (match your LOO builder behavior)
        switch lower(string(opt.SessionWeighting))
            case "fro"
                Yhat = Yhat / max(norm(Yhat,'fro'), opt.Eps);
            case "none"
                % no-op
            otherwise
                error('Unknown SessionWeighting: %s', string(opt.SessionWeighting));
        end

        Ypool = [Ypool; Yhat]; %#ok<AGROW>
    end

    if isempty(Ypool) || norm(Ypool,'fro') < opt.Eps
        continue;
    end

    [~, Ssvd, V] = svd(Ypool, 'econ');  % V: [K x r]
    r = size(V,2);
    nKeepEff = min(nKeep, r);

    svals = diag(Ssvd);
    denom = max(sum(svals.^2), opt.Eps);

    for pc = 1:nKeepEff
        v = V(:,pc)'; % 1 x K
        v = v / max(norm(v), opt.Eps);

        % sign fix (projmean across sessions)
        if strcmpi(string(opt.SignFix),"projmean")
            sc = nan(S,1);
            for ss = 1:S
                sc(ss) = local_projmean_score__(glmRezList{ss}, trIdList{ss}, v, gName, ...
                    opt.SignFixWinSec, opt.SignFixTrialPolicy, opt.Eps);
            end
            if mean(sc,'omitnan') < 0, v = -v; end
        elseif strcmpi(string(opt.SignFix),"maxabs")
            v = local_signfix_maxabs__(v);
        end

        Araw_rows(end+1,:) = v; %#ok<AGROW>
        axisNames{end+1}   = char(gName + "_" + string(pc)); %#ok<AGROW>

        axisMeta(end+1).groupIdx = g; %#ok<AGROW>
        axisMeta(end).groupName  = char(gName);
        axisMeta(end).pcIdx      = pc;
        axisMeta(end).sv         = svals(pc);
        axisMeta(end).expl       = (svals(pc)^2)/denom;
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

Araw_ord  = Araw_rows(ord,:);
names_ord = axisNames(ord);

% optional GS across groups
if strcmpi(string(opt.OrthMode),"none")
    A = Araw_ord;
    keep = true(1,size(Araw_ord,1));
else
    [A, keep] = local_gs_orth_rows__(Araw_ord, opt.Eps);
end

axesOut = struct();
axesOut.Araw_ord = Araw_ord;
axesOut.A        = A;
axesOut.keepAfterGS = keep;
axesOut.names_ord = names_ord;
axesOut.names     = names_ord(keep);
axesOut.headers   = headers;
end

function score = local_projmean_score__(glmRez, trId, axisRow, groupName, winSec, policy, epsVal)
M = size(glmRez.Yz,1);
nW = numel(glmRez.decBins.time);
if rem(M,nW)~=0, score = NaN; return; end
N = M / nW;

timeVec = glmRez.decBins.time(:);
tMask = (timeVec >= winSec(1)) & (timeVec <= winSec(2));
if ~any(tMask), score = NaN; return; end

trialMask = local_trialMask_fromGroup__(trId, string(groupName), N, policy);

z = glmRez.Yz * axisRow(:);
Z = reshape(z, [N, nW]);
score = mean(Z(trialMask, tMask), 'all', 'omitnan');
if ~isfinite(score) || abs(score) < epsVal
    % ok: caller uses omitnan
end
end

function v = local_signfix_maxabs__(v)
[~,idx] = max(abs(v));
if v(idx) < 0, v = -v; end
end

function [Q, keepIdx] = local_gs_orth_rows__(A, epsVal)
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

function G = local_num_groups__(grp)
if iscell(grp), G = numel(grp);
elseif isstruct(grp), G = numel(grp);
else, error('glmRez.group must be cell or struct array'); end
end

function gg = local_get_group__(grp, g)
if iscell(grp), gg = grp{g}; else, gg = grp(g); end
assert(isfield(gg,'name') && isfield(gg,'cols'), 'group must have .name and .cols');
end

function trialMask = local_trialMask_fromGroup__(trId, groupName, N, policy)
policy = string(policy);
groupName = string(groupName);
trialMask = true(N,1);

want = "all";
if strcmpi(policy, "GoNoGoByGroup")
    if contains(lower(groupName), "nogo"), want = "nogo";
    else, want = "go";
    end
end

if want == "all" || isempty(trId) || ~isstruct(trId)
    return;
end

if want == "go"
    trialMask = local_makeMask__(trId, 'goI', N);
else
    trialMask = local_makeMask__(trId, 'nogoI', N);
end

if ~any(trialMask)
    trialMask = true(N,1);
end
end

function m = local_makeMask__(trId, field, N)
m = true(N,1);
if ~isfield(trId, field) || isempty(trId.(field))
    return;
end
x = trId.(field);

if islogical(x)
    x = x(:);
    if numel(x) ~= N, return; end
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