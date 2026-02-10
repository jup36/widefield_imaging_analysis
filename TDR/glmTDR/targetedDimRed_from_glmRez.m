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
%   This is consistent with constructions like:
%       Ybig = cell2mat(Hs.Yw')   where Hs.Yw is 1xnW cell, each cell is N x K
%
% REQUIRED glmRez fields
%   beta        : [P x K]   coefficients (standardized X and standardized Y)
%   X_design    : [M x P]   raw design matrix (pre-standardization)
%   muX         : [1 x P]   mean used to build Xz during fitting (post-drop)
%   sdX         : [1 x P]   std  used to build Xz during fitting (post-drop)
%   Yz          : [M x K]   z-scored motif targets used during fitting
%   group       : 1xG struct array, each struct with: .name, .cols
%   decBins.time: vector used to infer nW
%
% OPTIONAL glmRez fields
%   Ybig        : [M x K]   raw motif activity (if you want ProjectWhichY="Ybig")
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
%   .Araw_rows     : [nAxis x K] raw axes (after sign-fix), ordered BEFORE GS
%   .Araw_ord      : [nAxis x K] raw axes after ordering
%   .A             : [G' x K] final axes after optional GS
%   .names_ord     : 1 x nAxis ordered axis names (pre-GS)
%   .names         : 1 x G' final axis names (post-GS)
%   .axisMeta      : struct array per axis (groupIdx, groupName, pcIdx, sv, expl)
%   .Xz            : [M x P] reconstructed standardized X
%   .Zrows         : [M x G'] projected rows
%   .Z             : [N x nW x G'] projected trajectories (trial x time x axis)
%   .nW, .N, .K, .P, .M
%
% Junchol Park lab-style: deterministic + robust + interpretability-preserving.

% -------------------- parse --------------------
p = inputParser;
p.addParameter('OrthMode',"GS",@(s)ischar(s)||isstring(s));
p.addParameter('SignFix',"maxabs",@(s)ischar(s)||isstring(s));
p.addParameter('Eps',1e-10,@(x)isscalar(x) && x>0);
p.addParameter('ProjectWhichY',"Yz",@(s)ischar(s)||isstring(s));

p.addParameter('MultiDimGroups', {'GoToneOn','NoGoToneOn','ToneOffGo','ToneOffNoGo'}, @(c)iscell(c));
p.addParameter('MultiDimK', 3, @(x)isnumeric(x)&&isscalar(x)&&x>=1);

p.addParameter('PriorityNames', {'GoToneOn','NoGoToneOn','ToneOffGo','ToneOffNoGo','Lick'}, @(c)iscell(c));
p.parse(varargin{:});
opt = p.Results;

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
grp  = glmRez.group;       % 1xG struct array

[M, P] = size(Xraw);
[Pb, K] = size(B);

assert(P==Pb, 'X_design has %d cols but beta has %d rows.', P, Pb);
assert(size(Yz,1)==M && size(Yz,2)==K, 'Yz must be [M x K] consistent with X_design/beta.');

% -------------------- reconstruct Xz (exact fitting space) --------------------
% (Even if X_design was already standardized upstream, re-standardizing here
% is safe because we use the SAME muX/sdX that were saved alongside beta.)
Xz = (Xraw - muX) ./ sdX;    % implicit expansion
Xz(~isfinite(Xz)) = 0;       % MUST match fitting behavior (NaNs -> 0 in z-space)

% -------------------- determine nW and N (TIME-MAJOR) --------------------
nW = numel(glmRez.decBins.time);
remVal = rem(M, nW);
if remVal ~= 0
    error('M=%d not divisible by nW=%d (remainder=%d). Check stacking/metadata.', M, nW, remVal);
end
N = M / nW;

% -------------------- build raw axes (2B: SVD of group-predicted component) --------------------
G = numel(grp);
groupNames = cell(1,G);
for g = 1:G
    groupNames{g} = char(string(grp{g}.name));
end


Araw_rows = [];     % [nAxis x K]
axisNames = {};     % 1 x nAxis
axisMeta  = struct('groupIdx',{},'groupName',{},'pcIdx',{},'sv',{},'expl',{});

multiNames = opt.MultiDimGroups;
multiK     = opt.MultiDimK;

for g = 1:G
    gName = string(grp{g}.name);
    cols  = grp{g}.cols(:)';

    cols  = cols(cols>=1 & cols<=P);
    if isempty(cols), continue; end

    Yhat_g = Xz(:,cols) * B(cols,:); % [M x K]
    if norm(Yhat_g,'fro') < opt.Eps, continue; end

    % how many dims to keep for this group?
    if any(strcmpi(gName, multiNames))
        nKeep = multiK;
    else
        nKeep = 1;
    end

    [~,S,V] = svd(Yhat_g, 'econ');   % V: [K x r]
    r = size(V,2);
    nKeep = min(nKeep, r);

    s = diag(S);
    denom = max(sum(s.^2), opt.Eps);

    for pc = 1:nKeep
        v = V(:,pc)';                % 1 x K (motif-space axis)
        if norm(v) < opt.Eps, continue; end

        v = v / max(norm(v), opt.Eps);
        v = sign_fix_axis(v, opt.SignFix);

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

% -------------------- axis ordering: group priority + pcIdx --------------------
% Desired: GoToneOn_1,2,3 -> NoGoToneOn_1,2,3 -> ToneOffGo_1,2,3 -> ToneOffNoGo_1,2,3 -> Lick_1 -> ...
priorityNames = opt.PriorityNames;

rankPerGroup = 1:G;  % default: natural group order
for i = 1:numel(priorityNames)
    gIdx = find(strcmp(groupNames, priorityNames{i}), 1, 'first');
    if ~isempty(gIdx)
        rankPerGroup(gIdx) = -1000 + i; % force priority groups to the front
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

Araw_ord      = Araw_rows(orderAxes,:);
names_ord     = axisNames(orderAxes);
axisMeta_ord  = axisMeta(orderAxes);

% -------------------- orthogonalize (optional) --------------------
if strcmpi(string(opt.OrthMode),"none")
    A = Araw_ord;
    keepLocal = true(1, size(Araw_ord,1));
    namesFinal = names_ord;
    axisMetaFinal = axisMeta_ord;
else
    [A, keepLocal] = gs_orth_rows(Araw_ord, opt.Eps);
    namesFinal = names_ord(keepLocal);
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
Zrows = Yproj * A';                % [M x G']
Gprime = size(A,1);

% TIME-MAJOR reshape (blocks of N rows per time bin)
Z = reshape(Zrows, [N, nW, Gprime]);  % N x T x G'

% -------------------- pack outputs --------------------
tdr = struct();
tdr.M = M; tdr.P = P; tdr.K = K; tdr.nW = nW; tdr.N = N;

tdr.Xz = Xz;

tdr.Araw_rows = Araw_rows;     % as-built (unordered)
tdr.axisNames = axisNames;     % as-built (unordered)
tdr.axisMeta  = axisMeta;      % as-built (unordered)

tdr.Araw_ord  = Araw_ord;      % ordered pre-GS
tdr.names_ord = names_ord;     % ordered pre-GS
tdr.axisMeta_ord = axisMeta_ord;

tdr.A     = A;                 % final axes (post-GS)
tdr.names = namesFinal;        % final axis names (post-GS)
tdr.axisMeta_final = axisMetaFinal;

tdr.Zrows = Zrows;
tdr.Z     = Z;

tdr.orderAxes = orderAxes;
tdr.keepAfterGS = keepLocal;

tdr.opt = opt;

end

% ========================== helpers ==========================
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
