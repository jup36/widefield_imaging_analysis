function axesOut = buildAnchorAxes_from_glmRezList(glmRezList, varargin)
%BUILDANCHORAXES_FROM_GLMREZLIST  Build a single shared ("anchor") GLM-TDR axis set from many glmRez.
%
% axesOut = buildAnchorAxes_from_glmRezList(glmRezList, 'Name', value, ...)
%
% Aggregation rule (per axis name):
%   - build per-session axes (NO per-session GS)
%   - sign-align across sessions
%   - mean + renormalize -> pooled axis
%   - finally order + optional GS

% -------------------- parse --------------------
p = inputParser;
p.addRequired('glmRezList', @(c) iscell(c) && ~isempty(c));

p.addParameter('OrthMode',"GS",@(s)ischar(s)||isstring(s));
p.addParameter('SignFix',"maxabs",@(s)ischar(s)||isstring(s));
p.addParameter('Eps',1e-10,@(x)isscalar(x)&&x>0);
p.addParameter('MultiDimGroups', {'GoToneOn','NoGoToneOn','ToneOffGo','ToneOffNoGo'}, @(c)iscell(c)||isstring(c));
p.addParameter('MultiDimK', 3, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('PriorityNames', {'GoToneOn','NoGoToneOn','ToneOffGo','ToneOffNoGo','Lick'}, @(c)iscell(c)||isstring(c));
p.addParameter('Verbose', true, @(x)islogical(x)&&isscalar(x));

p.parse(glmRezList, varargin{:});
opt = p.Results;

glmRezList = glmRezList(:);
S = numel(glmRezList);

% -------------------- run per-session TDR --------------------
tdrBySession = cell(S,1);
sessionKeep  = false(S,1);

for s = 1:S
    gr = glmRezList{s};
    if isempty(gr) || ~isstruct(gr)
        continue;
    end

    try
        % IMPORTANT: no GS per session; GS happens after pooling
        tdrBySession{s} = targetedDimRed_from_glmRez(gr, ...
            'OrthMode',       "none", ...
            'SignFix',        opt.SignFix, ...
            'Eps',            opt.Eps, ...
            'ProjectWhichY',  "Yz", ...  % irrelevant for axis extraction
            'MultiDimGroups', cellstr(string(opt.MultiDimGroups)), ...
            'MultiDimK',      opt.MultiDimK, ...
            'PriorityNames',  cellstr(string(opt.PriorityNames)));
        sessionKeep(s) = true;

    catch ME
        if opt.Verbose
            warning('buildAnchorAxes_from_glmRezList:SessionFail', ...
                'Skipping session %d/%d (targetedDimRed_from_glmRez failed): %s', s, S, ME.message);
        end
        sessionKeep(s) = false;
    end
end

tdrGood = tdrBySession(sessionKeep);
if isempty(tdrGood)
    error('buildAnchorAxes_from_glmRezList:NoValidSessions', ...
        'No glmRez entries produced valid targetedDimRed outputs.');
end

% -------------------- union of axis names --------------------
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
    for i = 1:nAvail
        v = Vc{i};
        if (v * vref') < 0
            v = -v;
        end
        V(i,:) = v;
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
A_pool    = A_pool(keepAxis,:);
names_kept = namesAll(keepAxis);
meta       = meta(keepAxis);

if isempty(A_pool)
    error('buildAnchorAxes_from_glmRezList:EmptyPooledAxes', ...
        'No axes survived pooling. Check naming consistency / groups / degeneracy.');
end

% -------------------- order axes: PriorityNames + pcIdx --------------------
priorityNames = cellstr(string(opt.PriorityNames));
rankMap = containers.Map('KeyType','char','ValueType','double');
for i = 1:numel(priorityNames)
    rankMap(priorityNames{i}) = i;
end

key1 = nan(numel(names_kept),1); % group rank
key2 = nan(numel(names_kept),1); % pcIdx
key3 = (1:numel(names_kept))';   % stable tiebreak

for i = 1:numel(names_kept)
    [gName, pcIdx] = parse_axis_name(string(names_kept{i}));
    if isKey(rankMap, gName)
        key1(i) = rankMap(gName);
    else
        key1(i) = 1e6;
    end
    key2(i) = pcIdx;
end

[~, ord] = sortrows([key1 key2 key3], [1 2 3]);

Araw_ord  = A_pool(ord,:);
names_ord = names_kept(ord);
meta_ord  = meta(ord);

% -------------------- optional GS on pooled axes --------------------
if strcmpi(string(opt.OrthMode), "none")
    A = Araw_ord;
    keepAfterGS = true(1, size(Araw_ord,1));
    namesFinal = names_ord;
    metaFinal  = meta_ord;
else
    [A, keepAfterGS] = gs_orth_rows(Araw_ord, opt.Eps);
    namesFinal = names_ord(keepAfterGS);
    metaFinal  = meta_ord(keepAfterGS);
end

% -------------------- pack output --------------------
axesOut = struct();
axesOut.Araw_rows    = A_pool;       % pooled, unordered
axesOut.names_raw    = names_kept;
axesOut.Araw_ord     = Araw_ord;     % pooled + ordered
axesOut.names_ord    = names_ord;
axesOut.axisMeta_ord = meta_ord;

axesOut.A              = A;          % final (post-GS)
axesOut.names          = namesFinal;
axesOut.axisMeta_final = metaFinal;
axesOut.keepAfterGS    = keepAfterGS;

axesOut.tdrBySession = tdrBySession;
axesOut.sessionKeep  = sessionKeep;

axesOut.opt = opt;

if opt.Verbose
    fprintf('[buildAnchorAxes] pooled axes=%d (kept after GS=%d) from %d/%d sessions\n', ...
        size(Araw_ord,1), size(A,1), sum(sessionKeep), S);
end

end


function v = sign_fix_axis(v, mode)
mode = lower(string(mode));
switch mode
    case "maxabs"
        [~,idx] = max(abs(v));
        if v(idx) < 0, v = -v; end
    case "none"
        % no-op
    otherwise
        % no-op
end
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
