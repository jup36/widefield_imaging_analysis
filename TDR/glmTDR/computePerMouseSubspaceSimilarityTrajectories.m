function rezSim = computePerMouseSubspaceSimilarityTrajectories(rezGroupSubspace, expertHeaders_perMouse, varargin)
%COMPUTEPERMOUSESUBSPACESIMILARITYTRAJECTORIES
% Build per-mouse LOO-safe expert reference subspaces and compute:
%   (1) session-wise similarity to expert reference subspace
%   (2) within-session pairwise subspace similarity
%
% NEW IN THIS VERSION
%   Also builds combined tone-epoch subspaces by UNION-SPAN:
%       GoTone   = orth([GoToneOn,   ToneOffGo])
%       NoGoTone = orth([NoGoToneOn, ToneOffNoGo])
%
%   These combined subspaces are treated exactly like the others for:
%       - similarity relative to expert reference
%       - within-session cross-trial-type similarity
%       - all-pairs similarity
%
% IMPORTANT
%   - Only sessions on/after the mouse-specific Day4 date are included if
%     'day4MarkC' is provided.
%   - "Similarity to expert" uses a leave-one-out reference when the target
%     session itself is part of the expert set.
%   - Since the order/sign of within-group singular vectors is arbitrary,
%     all comparisons are done at the SUBSPACE level, not vector-by-vector.
%
% INPUTS
%   rezGroupSubspace       : struct from extractPerSessionGroupSubspaces
%   expertHeaders_perMouse : cell/string matrix of expert session headers
%                            (rows = mice, columns = expert sessions)
%
% NAME-VALUE
%   'GroupNames'           : default = {'GoToneOn','NoGoToneOn','ToneOffGo','ToneOffNoGo'}
%   'IncludeCombinedTone'  : true (default)
%   'CombinedToneDefs'     : default =
%                              {
%                                'GoTone',   {'GoToneOn','ToneOffGo'}
%                                'NoGoTone', {'NoGoToneOn','ToneOffNoGo'}
%                              }
%   'day4MarkC'            : [] (default) or Nx2 cell/string:
%                              {mouseId, "MMDDYY"; ...}
%   'MinRefSessions'       : 2 (default)
%                            If LOO leaves too few expert sessions, supplement
%                            with nearest eligible sessions from the same mouse.
%   'CrossTrialPairList'   : default = {
%                               {'GoToneOn','NoGoToneOn'}
%                               {'ToneOffGo','ToneOffNoGo'}
%                               {'GoTone','NoGoTone'}   % only used if IncludeCombinedTone=true
%                             }
%   'DoAllPairs'           : true (default)
%                            Also compute all pairwise within-session similarities
%                            among active group names.
%   'Verbose'              : true (default)
%   'VerboseEvery'         : 1 (default)
%
% OUTPUT
%   rezSim.perMouse{m}
%       .mouseId
%       .sessions.headers
%       .sessions.date
%       .sessions.nSess
%       .sessions.isEligible
%       .sessions.day4Date
%       .sessions.isExpertEligible
%       .groupNamesBase
%       .groupNamesActive
%       .combinedDefs
%       .combinedSubspace{s}.(combinedName)
%           .U
%           .sourceNames
%           .dim
%       .refFull.(groupName)
%           .U
%           .nSourceSess
%           .sourceHeaders
%       .refLOO{s}.(groupName)
%           .U
%           .nSourceSess
%           .sourceHeaders
%           .didSupplement
%       .relativeToExpert.(groupName)
%           .overlap          : [nSess x 1]
%           .principalAngles  : {nSess x 1}
%       .withinSession.crossTrialType
%           .pairNames        : cellstr
%           .overlap          : [nSess x nPairs]
%           .principalAngles  : cell(nSess, nPairs)
%       .withinSession.allPairs
%           .pairNames        : cellstr
%           .overlap          : [nSess x nPairs]
%           .principalAngles  : cell(nSess, nPairs)

%% -------------------- parse --------------------
ip = inputParser;
ip.FunctionName = mfilename;

addRequired(ip, 'rezGroupSubspace', @(x) isstruct(x) && isfield(x,'perMouse'));
addRequired(ip, 'expertHeaders_perMouse', @(x) iscell(x) || isstring(x));

addParameter(ip, 'GroupNames', {'GoToneOn','NoGoToneOn','ToneOffGo','ToneOffNoGo'}, ...
    @(x) iscell(x) || isstring(x));

addParameter(ip, 'IncludeCombinedTone', true, @(x) islogical(x) && isscalar(x));
addParameter(ip, 'CombinedToneDefs', ...
    {'GoTone', {'GoToneOn','ToneOffGo'}; ...
     'NoGoTone', {'NoGoToneOn','ToneOffNoGo'}}, ...
    @(x) iscell(x) && size(x,2)==2);

addParameter(ip, 'day4MarkC', [], @(x) isempty(x) || iscell(x) || isstring(x));
addParameter(ip, 'MinRefSessions', 2, @(x) isnumeric(x) && isscalar(x) && x >= 1);

addParameter(ip, 'CrossTrialPairList', ...
    {{'GoToneOn','NoGoToneOn'}; {'ToneOffGo','ToneOffNoGo'}; {'GoTone','NoGoTone'}}, ...
    @(x) iscell(x));

addParameter(ip, 'DoAllPairs', true, @(x) islogical(x) && isscalar(x));
addParameter(ip, 'Verbose', true, @(x) islogical(x) && isscalar(x));
addParameter(ip, 'VerboseEvery', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);

parse(ip, rezGroupSubspace, expertHeaders_perMouse, varargin{:});
P = ip.Results;

groupNamesBase = cellstr(string(P.GroupNames(:)'));
combinedDefs = normalizeCombinedDefs_sim_(P.CombinedToneDefs);
crossTrialPairList = normalizePairList_sim_(P.CrossTrialPairList);

% Active groups = base + combined, if requested
if P.IncludeCombinedTone
    combinedNames = combinedDefs(:,1)';
    groupNamesActive = [groupNamesBase, combinedNames];
else
    combinedNames = {};
    groupNamesActive = groupNamesBase;
end

% Remove any crossTrial pairs that reference inactive groups
crossTrialPairList = filterPairListToActiveGroups_sim_(crossTrialPairList, groupNamesActive);

%% -------------------- normalize expert header matrix --------------------
expertHdrMat = normalizeHeaderMatrix_sim_(expertHeaders_perMouse);
expertMouseIds = cell(size(expertHdrMat,1),1);
for i = 1:size(expertHdrMat,1)
    expertMouseIds{i} = inferMouseIdFromRow_sim_(expertHdrMat(i,:));
end

%% -------------------- normalize day4 map --------------------
day4Map = containers.Map('KeyType','char','ValueType','char');
if ~isempty(P.day4MarkC)
    d4 = P.day4MarkC;
    if isstring(d4), d4 = cellstr(d4); end
    assert(iscell(d4) && size(d4,2)==2, ...
        'day4MarkC must be Nx2 cell/string array: {mouseId, "MMDDYY"; ...}.');
    for i = 1:size(d4,1)
        k = char(string(d4{i,1}));
        v = char(string(d4{i,2}));
        if ~isempty(k) && ~isempty(v)
            day4Map(k) = v;
        end
    end
end

%% -------------------- output shell --------------------
nMouse = numel(rezGroupSubspace.perMouse);
rezSim = struct();
rezSim.opt = P;
rezSim.perMouse = cell(nMouse,1);

%% ============================================================
% main loop over mice
% ============================================================
for m = 1:nMouse
    pm = rezGroupSubspace.perMouse{m};
    mouseId = char(string(pm.mouseId));

    if P.Verbose && (mod(m,P.VerboseEvery)==0 || m==1 || m==nMouse)
        fprintf('[computePerMouseSubspaceSimilarityTrajectories] mouse %d/%d (%s)\n', ...
            m, nMouse, mouseId);
    end

    headers = string(pm.sessions.headers(:));
    nSess   = numel(headers);

    if isfield(pm.sessions,'date') && numel(pm.sessions.date)==nSess
        sessDate = pm.sessions.date(:);
    else
        [sessDate, okDt] = parseHeaderDatesOnly_sim_(headers);
        if ~okDt
            sessDate = NaT(nSess,1);
        end
    end

    % ---------- day4 eligibility ----------
    if isKey(day4Map, mouseId)
        day4Date = datetime(day4Map(mouseId), 'InputFormat', 'MMddyy');
        isEligible = sessDate >= day4Date;
    else
        day4Date = NaT;
        isEligible = true(nSess,1);
    end

    % ---------- expert headers for this mouse ----------
    expertRow = {};
    iExpertRow = find(strcmpi(expertMouseIds, mouseId), 1, 'first');
    if ~isempty(iExpertRow)
        expertRow = expertHdrMat(iExpertRow,:);
    end
    expertHdr = string(expertRow);
    expertHdr = expertHdr(strlength(expertHdr)>0);

    % Keep only expert headers on/after day4 if day4 exists
    if ~isnat(day4Date) && ~isempty(expertHdr)
        [expertDt, okEDt] = parseHeaderDatesOnly_sim_(expertHdr);
        if okEDt
            expertHdr = expertHdr(expertDt >= day4Date);
        end
    end

    isExpertEligible = ismember(headers, expertHdr) & isEligible;

    outM = struct();
    outM.mouseId = mouseId;
    outM.groupNamesBase = groupNamesBase;
    outM.groupNamesActive = groupNamesActive;
    outM.combinedDefs = combinedDefs;

    outM.sessions = struct();
    outM.sessions.headers = cellstr(headers);
    outM.sessions.date = sessDate;
    outM.sessions.nSess = nSess;
    outM.sessions.isEligible = isEligible;
    outM.sessions.day4Date = day4Date;
    outM.sessions.isExpertEligible = isExpertEligible;

    %% ---------- build combined session subspaces ----------
    outM.combinedSubspace = cell(nSess,1);
    if P.IncludeCombinedTone
        for s = 1:nSess
            outM.combinedSubspace{s} = struct();
            for c = 1:size(combinedDefs,1)
                cName = combinedDefs{c,1};
                srcNames = combinedDefs{c,2};

                Ulist = cell(numel(srcNames),1);
                for j = 1:numel(srcNames)
                    Ulist{j} = getUfromSessionSubspace_sim_(pm.subspace, s, srcNames{j});
                end

                Uc = combineSubspaces_union_sim_(Ulist);

                outM.combinedSubspace{s}.(cName) = struct( ...
                    'U', Uc, ...
                    'sourceNames', {srcNames}, ...
                    'dim', size(Uc,2));
            end
        end
    else
        for s = 1:nSess
            outM.combinedSubspace{s} = struct();
        end
    end

    %% ---------- build full expert reference per active group ----------
    outM.refFull = struct();

    for g = 1:numel(groupNamesActive)
        gName = groupNamesActive{g};
        Ucell = collectGroupUcell_active_sim_(pm.subspace, outM.combinedSubspace, gName, isExpertEligible, combinedNames);
        refU = averageSubspaces_projection_sim_(Ucell);

        srcHdr = headers(isExpertEligible);
        outM.refFull.(gName) = struct( ...
            'U', refU, ...
            'nSourceSess', numel(Ucell), ...
            'sourceHeaders', cellstr(srcHdr(:)) );
    end

    %% ---------- build session-wise LOO refs + expert similarity ----------
    outM.refLOO = cell(nSess,1);
    outM.relativeToExpert = struct();

    for g = 1:numel(groupNamesActive)
        gName = groupNamesActive{g};
        overlapVec = nan(nSess,1);
        thetaCell  = cell(nSess,1);

        for s = 1:nSess
            outM.refLOO{s}.(gName) = struct( ...
                'U', [], ...
                'nSourceSess', 0, ...
                'sourceHeaders', {{}}, ...
                'didSupplement', false);
        end

        for s = 1:nSess
            if ~isEligible(s)
                continue;
            end

            targetU = getUactiveFromSession_sim_(pm.subspace, outM.combinedSubspace, s, gName, combinedNames);
            if isempty(targetU)
                continue;
            end

            trainMask = isExpertEligible;

            % LOO exclusion if target session itself is in expert set
            if isExpertEligible(s)
                trainMask(s) = false;
            end

            didSupplement = false;

            % If too few ref sessions remain, supplement with nearest eligible sessions
            if sum(trainMask) < P.MinRefSessions
                supplementMask = nearestEligibleSupplementMask_sim_( ...
                    s, isEligible, trainMask, P.MinRefSessions);
                if any(supplementMask)
                    trainMask = trainMask | supplementMask;
                    didSupplement = true;
                end
            end

            Ucell = collectGroupUcell_active_sim_(pm.subspace, outM.combinedSubspace, gName, trainMask, combinedNames);
            refU = averageSubspaces_projection_sim_(Ucell);

            srcHdr = headers(trainMask);

            outM.refLOO{s}.(gName) = struct( ...
                'U', refU, ...
                'nSourceSess', numel(Ucell), ...
                'sourceHeaders', cellstr(srcHdr(:))', ...
                'didSupplement', didSupplement);

            overlapVec(s) = subspaceOverlap_sim_(targetU, refU);
            thetaCell{s}  = principalAngles_sim_(targetU, refU);
        end

        outM.relativeToExpert.(gName) = struct( ...
            'overlap', overlapVec, ...
            'principalAngles', {thetaCell});
    end

    %% ---------- within-session cross-trial-type similarity ----------
    outM.withinSession = struct();

    ctPairNames = pairNamesFromList_sim_(crossTrialPairList);
    nCtPairs = size(crossTrialPairList,1);
    ctOverlap = nan(nSess, nCtPairs);
    ctTheta   = cell(nSess, nCtPairs);

    for s = 1:nSess
        if ~isEligible(s)
            continue;
        end

        for p = 1:nCtPairs
            g1 = crossTrialPairList{p,1};
            g2 = crossTrialPairList{p,2};

            U1 = getUactiveFromSession_sim_(pm.subspace, outM.combinedSubspace, s, g1, combinedNames);
            U2 = getUactiveFromSession_sim_(pm.subspace, outM.combinedSubspace, s, g2, combinedNames);

            ctOverlap(s,p) = subspaceOverlap_sim_(U1, U2);
            ctTheta{s,p} = principalAngles_sim_(U1, U2);
        end
    end

    outM.withinSession.crossTrialType = struct( ...
        'pairNames', {ctPairNames}, ...
        'overlap', ctOverlap, ...
        'principalAngles', {ctTheta});

    %% ---------- all-pairs within-session similarity ----------
    if P.DoAllPairs
        allPairList = makeAllPairs_sim_(groupNamesActive);
        allPairNames = pairNamesFromList_sim_(allPairList);
        nAllPairs = size(allPairList,1);

        allOverlap = nan(nSess, nAllPairs);
        allTheta   = cell(nSess, nAllPairs);

        for s = 1:nSess
            if ~isEligible(s)
                continue;
            end

            for p = 1:nAllPairs
                g1 = allPairList{p,1};
                g2 = allPairList{p,2};

                U1 = getUactiveFromSession_sim_(pm.subspace, outM.combinedSubspace, s, g1, combinedNames);
                U2 = getUactiveFromSession_sim_(pm.subspace, outM.combinedSubspace, s, g2, combinedNames);

                allOverlap(s,p) = subspaceOverlap_sim_(U1, U2);
                allTheta{s,p} = principalAngles_sim_(U1, U2);
            end
        end

        outM.withinSession.allPairs = struct( ...
            'pairNames', {allPairNames}, ...
            'overlap', allOverlap, ...
            'principalAngles', {allTheta});
    end

    rezSim.perMouse{m} = outM;
end

end

%% ======================================================================
% Helpers
% ======================================================================

function defs = normalizeCombinedDefs_sim_(defsIn)
assert(iscell(defsIn) && size(defsIn,2)==2, ...
    'CombinedToneDefs must be Nx2 cell array: {combinedName, {source1 source2 ...}; ...}');
defs = cell(size(defsIn,1),2);
for i = 1:size(defsIn,1)
    defs{i,1} = char(string(defsIn{i,1}));
    src = defsIn{i,2};
    assert(iscell(src) || isstring(src), 'Second column of CombinedToneDefs must be cell/string array.');
    defs{i,2} = cellstr(string(src(:)'));
end
end

function pairList = filterPairListToActiveGroups_sim_(pairListIn, activeNames)
pairList = cell(0,2);
activeNames = cellstr(string(activeNames(:)'));
for i = 1:size(pairListIn,1)
    g1 = char(string(pairListIn{i,1}));
    g2 = char(string(pairListIn{i,2}));
    if ismember(g1, activeNames) && ismember(g2, activeNames)
        pairList(end+1,1:2) = {g1, g2}; %#ok<AGROW>
    end
end
end

function pairList = normalizePairList_sim_(pairIn)
if isempty(pairIn)
    pairList = cell(0,2);
    return;
end

assert(iscell(pairIn), 'CrossTrialPairList must be a cell array.');

if size(pairIn,2) == 1
    pairList = cell(size(pairIn,1),2);
    for i = 1:size(pairIn,1)
        row = pairIn{i};
        assert(iscell(row) && numel(row)==2, ...
            'Each row of CrossTrialPairList must contain two group names.');
        pairList{i,1} = char(string(row{1}));
        pairList{i,2} = char(string(row{2}));
    end
elseif size(pairIn,2) == 2
    pairList = cell(size(pairIn,1),2);
    for i = 1:size(pairIn,1)
        pairList{i,1} = char(string(pairIn{i,1}));
        pairList{i,2} = char(string(pairIn{i,2}));
    end
else
    error('CrossTrialPairList must be Nx1 of 1x2 cells, or Nx2 cell array.');
end
end

function names = pairNamesFromList_sim_(pairList)
n = size(pairList,1);
names = cell(n,1);
for i = 1:n
    names{i} = sprintf('%s__vs__%s', char(string(pairList{i,1})), char(string(pairList{i,2})));
end
end

function pairList = makeAllPairs_sim_(groupNames)
n = numel(groupNames);
pairList = cell(n*(n-1)/2, 2);
c = 0;
for i = 1:n-1
    for j = i+1:n
        c = c + 1;
        pairList{c,1} = groupNames{i};
        pairList{c,2} = groupNames{j};
    end
end
end

function U = getUfromSessionSubspace_sim_(subspaceCell, s, groupName)
U = [];
if s > numel(subspaceCell) || isempty(subspaceCell{s})
    return;
end
if ~isfield(subspaceCell{s}, groupName)
    return;
end
sub = subspaceCell{s}.(groupName);
if isempty(sub) || ~isfield(sub,'U') || isempty(sub.U)
    return;
end
U = sub.U;
end

function U = getUfromCombinedSubspace_sim_(combinedCell, s, groupName)
U = [];
if s > numel(combinedCell) || isempty(combinedCell{s})
    return;
end
if ~isfield(combinedCell{s}, groupName)
    return;
end
sub = combinedCell{s}.(groupName);
if isempty(sub) || ~isfield(sub,'U') || isempty(sub.U)
    return;
end
U = sub.U;
end

function U = getUactiveFromSession_sim_(baseSubspaceCell, combinedCell, s, groupName, combinedNames)
if ismember(groupName, combinedNames)
    U = getUfromCombinedSubspace_sim_(combinedCell, s, groupName);
else
    U = getUfromSessionSubspace_sim_(baseSubspaceCell, s, groupName);
end
end

function Ucell = collectGroupUcell_active_sim_(baseSubspaceCell, combinedCell, groupName, mask, combinedNames)
Ucell = {};
idx = find(mask(:)');
for i = 1:numel(idx)
    s = idx(i);
    U = getUactiveFromSession_sim_(baseSubspaceCell, combinedCell, s, groupName, combinedNames);
    if ~isempty(U)
        Ucell{end+1,1} = U; %#ok<AGROW>
    end
end
end

function Uc = combineSubspaces_union_sim_(Ulist)
% Combine multiple subspaces by union-span: orth([U1 U2 ...])
Uc = [];
for i = 1:numel(Ulist)
    U = Ulist{i};
    if isempty(U)
        continue;
    end
    if isempty(Uc)
        Uc = U;
    else
        Uc = [Uc U]; %#ok<AGROW>
    end
end

if isempty(Uc)
    return;
end

Uc = orth(Uc);
end

function Ubar = averageSubspaces_projection_sim_(Ucell)
% Average subspaces by averaging projection matrices and taking top eigvecs.
% Target dimensionality is determined adaptively from the median dimension.
if isempty(Ucell)
    Ubar = [];
    return;
end

K = size(Ucell{1},1);
dims = zeros(numel(Ucell),1);
Psum = zeros(K,K);
nUsed = 0;

for i = 1:numel(Ucell)
    U = Ucell{i};
    if isempty(U)
        continue;
    end
    dims(i) = size(U,2);
    Psum = Psum + (U * U');
    nUsed = nUsed + 1;
end

if nUsed == 0 || norm(Psum,'fro') < eps
    Ubar = [];
    return;
end

dims = dims(dims>0);
kTarget = max(1, round(median(dims)));

Pavg = Psum / nUsed;
Pavg = (Pavg + Pavg') / 2;

[V,D] = eig(Pavg);
[d, ord] = sort(diag(D), 'descend');
V = V(:, ord);

kEff = min([kTarget, size(V,2), sum(d > eps)]);
if kEff < 1
    Ubar = [];
    return;
end

Ubar = V(:,1:kEff);
end

function supplementMask = nearestEligibleSupplementMask_sim_(s, isEligible, trainMask, minRefSessions)
supplementMask = false(size(trainMask));

if sum(trainMask) >= minRefSessions
    return;
end

nSess = numel(trainMask);
dist = abs((1:nSess)' - s);

candidateMask = isEligible & ~trainMask;
candidateMask(s) = false;

cand = find(candidateMask);
if isempty(cand)
    return;
end

[~, ord] = sort(dist(cand), 'ascend');
cand = cand(ord);

need = minRefSessions - sum(trainMask);
need = max(need, 0);
if need > 0
    cand = cand(1:min(need, numel(cand)));
    supplementMask(cand) = true;
end
end

function sim = subspaceOverlap_sim_(U1, U2)
if isempty(U1) || isempty(U2)
    sim = NaN;
    return;
end
d = min(size(U1,2), size(U2,2));
M = U1' * U2;
sim = norm(M, 'fro')^2 / d;
end

function theta = principalAngles_sim_(U1, U2)
if isempty(U1) || isempty(U2)
    theta = NaN;
    return;
end
s = svd(U1' * U2);
s = min(max(s, -1), 1);
theta = acos(s);
end

function hdrMat = normalizeHeaderMatrix_sim_(hdrIn)
if isempty(hdrIn)
    hdrMat = cell(0,0);
    return;
end
if isstring(hdrIn)
    hdrIn = cellstr(hdrIn);
end
if iscell(hdrIn) && isvector(hdrIn)
    hdrIn = reshape(hdrIn, [], 1);
end

hdrMat = cell(size(hdrIn));
for i = 1:numel(hdrIn)
    x = hdrIn{i};
    if isempty(x)
        hdrMat{i} = [];
    elseif isstring(x) || ischar(x)
        s = char(string(x));
        if strlength(string(s))==0
            hdrMat{i} = [];
        else
            hdrMat{i} = s;
        end
    else
        hdrMat{i} = [];
    end
end
end

function mouseId = inferMouseIdFromRow_sim_(hdrRow)
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

function [dtOnly, ok] = parseHeaderDatesOnly_sim_(hdrS)
hdrS = string(hdrS(:));
n = numel(hdrS);
dtOnly = NaT(n,1);
ok = true;

for i = 1:n
    h = char(hdrS(i));
    tok = regexp(h, '_(\d{6})(?:-\d+)?$', 'tokens', 'once');
    if isempty(tok)
        ok = false;
        return;
    end
    try
        dtOnly(i) = datetime(tok{1}, 'InputFormat','MMddyy');
    catch
        ok = false;
        return;
    end
end
end