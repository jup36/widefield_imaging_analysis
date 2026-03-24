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

%% ========================== helpers ==========================
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