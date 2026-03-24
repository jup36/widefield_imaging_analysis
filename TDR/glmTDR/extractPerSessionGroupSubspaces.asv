function rez = extractPerSessionGroupSubspaces(headerC, glmRezC, trIdC, varargin)
%EXTRACTPERSESSIONGROUPSUBSPACES
% Extract session-specific motif-space subspaces for selected GLM predictor groups.
%
% rez = extractPerSessionGroupSubspaces(headerC, glmRezC, trIdC, 'Name', value, ...)
%
% REQUIRED
%   headerC : cell (J x S) of session headers
%   glmRezC : cell (J x S) of glmRez structs
%   trIdC   : cell (J x S) of trial-type structs
%
% NAME-VALUE
%   'GroupNames'         : default = {'GoToneOn','NoGoToneOn','ToneOffGo','ToneOffNoGo'}
%   'KPerGroup'          : default = 3
%   'TrialPolicy'        : "GoNoGoByGroup" (default) | "all"
%   'SessionWeighting'   : "none" (default) | "fro"
%   'Verbose'            : true
%   'VerboseEvery'       : 1
%   'selectMice'         : [] (default) or cell/string array of mouse IDs
%                          e.g. {'m1044','m1045','m1048','m1049','m1092','m1094','m1613','m1859','m1873'}
%                          If provided, ONLY these mice are processed, and
%                          rez.perMouse is returned in exactly this order.
%
% OUTPUT
%   rez.perMouse{m}
%       .mouseId
%       .sessions.headers
%       .sessions.date
%       .sessions.nSess
%       .groupNames
%       .kPerGroup
%       .subspace{s}.(groupName)
%           .U           : [Kmotif x kEff] right singular vectors
%           .sval        : singular values kept
%           .explained   : variance explained within group Yhat
%           .nTrial      : number of selected trials
%           .nTime       : number of time bins
%           .cols        : predictor cols used
%           .groupName
%
% NOTE
%   The columns of U define a subspace. Their order/sign should NOT be
%   interpreted across sessions without alignment. Use subspace metrics later.
%

%% Parse
ip = inputParser;
ip.FunctionName = mfilename;

addRequired(ip, 'headerC', @(x) iscell(x));
addRequired(ip, 'glmRezC', @(x) iscell(x) && isequal(size(x), size(headerC)));
addRequired(ip, 'trIdC',   @(x) iscell(x) && isequal(size(x), size(headerC)));

addParameter(ip, 'GroupNames', {'GoToneOn','NoGoToneOn','ToneOffGo','ToneOffNoGo'}, ...
    @(x) iscell(x) || isstring(x));
addParameter(ip, 'KPerGroup', 3, @(x) isnumeric(x) && isscalar(x) && x>=1);
addParameter(ip, 'TrialPolicy', "GoNoGoByGroup", @(x) ischar(x) || isstring(x));
addParameter(ip, 'SessionWeighting', "none", @(x) ischar(x) || isstring(x));
addParameter(ip, 'Verbose', true, @(x) islogical(x) && isscalar(x));
addParameter(ip, 'VerboseEvery', 1, @(x) isnumeric(x) && isscalar(x) && x>=1);
addParameter(ip, 'selectMice', [], @(x) isempty(x) || iscell(x) || isstring(x));

parse(ip, headerC, glmRezC, trIdC, varargin{:});
P = ip.Results;

groupNames = cellstr(string(P.GroupNames));
kPerGroup  = P.KPerGroup;
trialPolicy = string(P.TrialPolicy);
sessionWeighting = lower(string(P.SessionWeighting));

if isempty(P.selectMice)
    selectMice = {};
else
    selectMice = cellstr(string(P.selectMice(:)));
end

%% Trim trailing empties
[headerC, glmRezC, trIdC] = trimTrailingAllEmptyCols_sub_(headerC, glmRezC, trIdC);

%% Build lookup
[uniqHdr, glmFirst, trFirst] = buildHeaderLookup_sub_(headerC, glmRezC, trIdC);

%% Determine mice and order
if isempty(selectMice)
    mouseIds = inferAllMouseIds_sub_(uniqHdr);
else
    mouseIds = selectMice(:);
end
nMouse = numel(mouseIds);

rez = struct();
rez.opt = P;
rez.perMouse = cell(nMouse,1);

for m = 1:nMouse
    mouseId = mouseIds{m};

    if P.Verbose && (mod(m, P.VerboseEvery)==0 || m==1 || m==nMouse)
        fprintf('[extractPerSessionGroupSubspaces] mouse %d/%d (%s)\n', m, nMouse, mouseId);
    end

    [sessHdr, sessGlm, sessTr] = collectMouseSessions_sub_(mouseId, uniqHdr, glmFirst, trFirst);
    nSess = numel(sessHdr);

    outM = struct();
    outM.mouseId = mouseId;
    outM.groupNames = groupNames;
    outM.kPerGroup = kPerGroup;

    outM.sessions = struct();
    outM.sessions.headers = cellstr(sessHdr(:));
    outM.sessions.nSess = nSess;

    [sessDt, okDt] = parseHeaderDatesOnly_sub_(sessHdr);
    if okDt
        outM.sessions.date = sessDt;
    else
        outM.sessions.date = NaT(nSess,1);
    end

    outM.subspace = cell(nSess,1);

    for s = 1:nSess
        gr = sessGlm{s};
        tr = sessTr{s};

        outS = struct();

        for g = 1:numel(groupNames)
            gName = groupNames{g};

            try
                sub = extractOneGroupSubspace_sub_(gr, tr, gName, kPerGroup, trialPolicy, sessionWeighting);
            catch ME
                warning('Failed subspace extraction | mouse=%s session=%s group=%s | %s', ...
                    mouseId, char(string(sessHdr(s))), gName, ME.message);
                sub = [];
            end

            outS.(gName) = sub;
        end

        outM.subspace{s} = outS;
    end

    rez.perMouse{m} = outM;
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% HELPER FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sub = extractOneGroupSubspace_sub_(glmRez, trId, groupName, kPerGroup, trialPolicy, sessionWeighting)
% Returns a struct with U = [Kmotif x kEff]

assert(isfield(glmRez,'X_design') && isfield(glmRez,'beta') && isfield(glmRez,'muX') && isfield(glmRez,'sdX') ...
    && isfield(glmRez,'group') && isfield(glmRez,'decBins') && isfield(glmRez.decBins,'time'), ...
    'glmRez missing required fields.');

% Find group columns
[groupCols, found] = getGroupCols_sub_(glmRez.group, groupName);
if ~found || isempty(groupCols)
    sub = [];
    return;
end

% Standardized design
Xz = (glmRez.X_design - glmRez.muX) ./ glmRez.sdX;
Xz(~isfinite(Xz)) = 0;

% Group-predicted activity in motif space
Yhat = Xz(:, groupCols) * glmRez.beta(groupCols, :);   % [M x Kmotif]

% Reshape to [trial x time x Kmotif]
timeVec = glmRez.decBins.time(:);
nTime = numel(timeVec);
M = size(Yhat,1);
assert(rem(M, nTime)==0, 'Yhat rows not divisible by nTime.');
nTrial = M / nTime;
Kmotif = size(Yhat,2);

Y3 = reshape(Yhat, [nTrial, nTime, Kmotif]);

% Apply trial selection
trialMask = trialMaskFromGroup_sub_(trId, string(groupName), nTrial, trialPolicy);
if isempty(trialMask) || ~any(trialMask)
    trialMask = true(nTrial,1);
end

Y3 = Y3(trialMask,:,:);
nTrialSel = size(Y3,1);

% Flatten back to rows x motifs for SVD
Y2 = reshape(Y3, [nTrialSel*nTime, Kmotif]);

% Optional session weighting
switch sessionWeighting
    case "fro"
        Y2 = Y2 / max(norm(Y2, 'fro'), eps);
    case "none"
        % no-op
    otherwise
        error('Unknown SessionWeighting: %s', sessionWeighting);
end

% SVD
if isempty(Y2) || norm(Y2,'fro') < eps
    sub = [];
    return;
end

[~, S, V] = svd(Y2, 'econ');
sval = diag(S);

kEff = min(kPerGroup, size(V,2));
U = V(:,1:kEff);    % [Kmotif x kEff]

expl = (sval(1:kEff).^2) / max(sum(sval.^2), eps);

sub = struct();
sub.U = U;
sub.sval = sval(1:kEff);
sub.explained = expl;
sub.nTrial = nTrialSel;
sub.nTime = nTime;
sub.cols = groupCols;
sub.groupName = char(groupName);
end


function sim = subspaceOverlap(U1, U2)
%SUBSPACEOVERLAP Rotation/sign/permutation-invariant overlap between two subspaces.
%
% U1, U2: [K x d1], [K x d2], assumed orthonormal columns
% sim = average squared cosine of principal angles
%
% Returns scalar in [0,1].

if isempty(U1) || isempty(U2)
    sim = NaN;
    return;
end

d = min(size(U1,2), size(U2,2));
M = U1' * U2;
sim = norm(M, 'fro')^2 / d;
end

function theta = principalAngles_sub_(U1, U2)
if isempty(U1) || isempty(U2)
    theta = NaN;
    return;
end

s = svd(U1' * U2);
s = min(max(s, -1), 1);
theta = acos(s);
end

function [groupCols, found] = getGroupCols_sub_(grp, groupName)

found = false;
groupCols = [];

if iscell(grp)
    G = numel(grp);
    getter = @(i) grp{i};
elseif isstruct(grp)
    G = numel(grp);
    getter = @(i) grp(i);
else
    error('glmRez.group must be cell or struct.');
end

for g = 1:G
    gg = getter(g);
    if strcmpi(string(gg.name), string(groupName))
        groupCols = unique(round(gg.cols(:)'),'stable');
        found = true;
        return;
    end
end
end


function trialMask = trialMaskFromGroup_sub_(trId, groupName, nTrial, policy)

policy = lower(string(policy));
groupName = lower(string(groupName));

trialMask = true(nTrial,1);

if isempty(trId) || ~isstruct(trId)
    return;
end

switch policy
    case "all"
        return

    case "gonogobygroup"
        if contains(groupName, "nogo")
            trialMask = makeMask_sub_(trId, 'nogoI', nTrial);
        else
            trialMask = makeMask_sub_(trId, 'goI', nTrial);
        end

    otherwise
        error('Unknown TrialPolicy: %s', policy);
end

if ~any(trialMask)
    trialMask = true(nTrial,1);
end
end


function m = makeMask_sub_(trId, field, N)

m = true(N,1);

if ~isfield(trId, field) || isempty(trId.(field))
    return;
end

x = trId.(field);

if islogical(x)
    x = x(:);
    if numel(x)==N
        m = x;
    end
else
    idx = unique(round(x(:)));
    idx = idx(idx>=1 & idx<=N);
    m = false(N,1);
    m(idx) = true;
end
end


function [sessHdr, sessGlm, sessTr] = collectMouseSessions_sub_(mouseId, uniqHdr, glmFirst, trFirst)

mouseId = string(mouseId);
mask = contains(uniqHdr, mouseId);

sessHdr = uniqHdr(mask);
sessGlm = glmFirst(mask);
sessTr  = trFirst(mask);

[dtKey, ok] = parseHeaderDatetimeWithSuffix_sub_(sessHdr);
if ok
    [~, ord] = sort(dtKey, 'ascend');
    sessHdr = sessHdr(ord);
    sessGlm = sessGlm(ord);
    sessTr  = sessTr(ord);
end
end


function [uniqHdr, glmFlatFirst, trFlatFirst] = buildHeaderLookup_sub_(headerC, glmRezC, trIdC)

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
    warning('Duplicate headers found; using first occurrence. Example: %s', string(dup(1)));
end

glmFlatFirst = cell(numel(uniqHdr),1);
trFlatFirst  = cell(numel(uniqHdr),1);

for i = 1:numel(uniqHdr)
    ii = find(hdrFlatS == uniqHdr(i), 1, 'first');
    glmFlatFirst{i} = glmFlat{ii};
    trFlatFirst{i}  = trFlat{ii};
end
end


function mouseIds = inferAllMouseIds_sub_(uniqHdr)

mouseS = strings(0,1);

for i = 1:numel(uniqHdr)
    tok = regexp(string(uniqHdr(i)), '(m\d{3,5})', 'tokens', 'once');
    if ~isempty(tok)
        mouseS(end+1,1) = string(tok{1}); %#ok<AGROW>
    end
end

mouseS = unique(mouseS, 'stable');
mouseIds = cellstr(mouseS);
end


function [dtKey, ok] = parseHeaderDatetimeWithSuffix_sub_(hdrS)

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

    d0 = datetime(mmddyy, 'InputFormat','MMddyy');
    dtKey(i) = d0 + seconds(suf);
end
end


function [dtOnly, ok] = parseHeaderDatesOnly_sub_(hdrS)

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


function [headerC, glmRezC, trIdC] = trimTrailingAllEmptyCols_sub_(headerC, glmRezC, trIdC)

[J,~] = size(headerC);
Smax = max([size(headerC,2), size(glmRezC,2), size(trIdC,2)]);

if size(headerC,2) < Smax, headerC(:,end+1:Smax) = {[]}; end
if size(glmRezC,2) < Smax, glmRezC(:,end+1:Smax) = {[]}; end
if size(trIdC,2) < Smax, trIdC(:,end+1:Smax) = {[]}; end

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
    headerC = cell(J,0);
    glmRezC = cell(J,0);
    trIdC   = cell(J,0);
else
    headerC = headerC(:,1:keepLast);
    glmRezC = glmRezC(:,1:keepLast);
    trIdC   = trIdC(:,1:keepLast);
end
end