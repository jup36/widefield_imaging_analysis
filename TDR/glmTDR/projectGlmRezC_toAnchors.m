function prj = projectGlmRezC_toAnchors(glmA, headerC, glmRezC, varargin)
%PROJECTGLMREZC_TOANCHORS  Project each session's motif activity onto:
%  (1) global expert anchor axes
%  (2) per-mouse expert anchor axes
%  (3) per-session axes (each session projected to its OWN orthonormalized A)
%
% prj = projectGlmRezC_toAnchors(glmA, headerC, glmRezC, 'Name', value, ...)
%
% INPUTS
%   glmA    : anchor struct from collectPerMouseAndGlobalGlmPrjAxes
%   headerC : (J x S) cell array of session headers
%   glmRezC : (J x S) cell array of glmRez structs (same size)
%
% NAME-VALUE
%   'ProjectWhichY' : "Yz" (default) | "Ybig"
%   'Eps'           : 1e-10 (default)
%   'Verbose'       : true (default)
%
% OUTPUT (prj struct)
%   .global
%       .A, .names         : global axes + names
%       .ZC                : (J x S) cell; [N x nW x nAxis]
%       .okMat             : logical (J x S)
%   .perMouse
%       .ZC                : (J x S) cell; [N x nW x nAxis_mouse]
%       .namesC            : (J x S) cell; names used
%       .okMat             : logical (J x S)
%       .mouseIdMat        : (J x S) cellstr inferred mouseId
%   .perSession
%       .ZC                : (J x S) cell; [N x nW x nAxis_session]
%       .namesC            : (J x S) cell; per-session names (matched to A rows)
%       .ArawC             : (J x S) cell; per-session Araw_ord (optional storage)
%       .okMat             : logical (J x S)
%
% STACKING ASSUMPTION
%   glmRez.Yz (and Ybig) are TIME-MAJOR in rows:
%     [time bin 1: trials 1..N], [time bin 2: trials 1..N], ..., [time bin nW]
%   Invert to trial x time x axis using:
%     Z = reshape(Zrows, [N, nW, nAxis]);
%
% JC lab-style: robust + explicit.

% -------------------- parse --------------------
p = inputParser;
p.addRequired('glmA',   @(s)isstruct(s) && ~isempty(s));
p.addRequired('headerC',@(c)iscell(c));
p.addRequired('glmRezC',@(c)iscell(c) && isequal(size(c), size(headerC)));

p.addParameter('ProjectWhichY',"Yz",@(s)ischar(s)||isstring(s));
p.addParameter('Eps',1e-10,@(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('Verbose',true,@(x)islogical(x)&&isscalar(x));

p.parse(glmA, headerC, glmRezC, varargin{:});
opt = p.Results;

[J,S] = size(headerC);

% -------------------- validate global axes schema --------------------
hasGlobal = isfield(glmA,'global') && isstruct(glmA.global) && ...
            isfield(glmA.global,'axes') && isstruct(glmA.global.axes) && ...
            isfield(glmA.global.axes,'A') && isfield(glmA.global.axes,'names') && ...
            ~isempty(glmA.global.axes.A);

if hasGlobal
    A_global = glmA.global.axes.A;
    names_global = cellstr(string(glmA.global.axes.names));
    nAxisG = size(A_global,1);
else
    A_global = [];
    names_global = {};
    nAxisG = 0;
end

% -------------------- allocate outputs --------------------
prj = struct();

prj.global = struct();
prj.global.A     = A_global;
prj.global.names = names_global;
prj.global.ZC    = cell(J,S);
prj.global.okMat = false(J,S);

prj.perMouse = struct();
prj.perMouse.ZC         = cell(J,S);
prj.perMouse.namesC     = cell(J,S);
prj.perMouse.okMat      = false(J,S);
prj.perMouse.mouseIdMat = cell(J,S);

prj.perSession = struct();
prj.perSession.ZC     = cell(J,S);
prj.perSession.namesC = cell(J,S);
prj.perSession.ArawC  = cell(J,S);   % optional, useful for debugging
prj.perSession.okMat  = false(J,S);

prj.meta = struct();
prj.meta.headerC = headerC;
prj.meta.ProjectWhichY = char(string(opt.ProjectWhichY));

% -------------------- main loop --------------------
nOKg = 0;
nOKm = 0;
nOKs = 0;

for j = 1:J
    for s = 1:S

        hdr = headerC{j,s};
        if isempty(hdr) || (isstring(hdr) && strlength(hdr)==0)
            continue;
        end

        gr = glmRezC{j,s};
        if isempty(gr) || ~isstruct(gr)
            continue;
        end

        % infer mouseId from header
        mouseId = infer_mouse_id_from_header(hdr);
        prj.perMouse.mouseIdMat{j,s} = mouseId;

        % choose Y for projection
        try
            [Yproj, nW, N, K] = get_project_target_Y(gr, opt.ProjectWhichY);
        catch ME
            if opt.Verbose
                warning('projectGlmRezC_toAnchors:YprojFail', ...
                    '[%s] Failed to get %s for projection: %s', char(string(hdr)), char(string(opt.ProjectWhichY)), ME.message);
            end
            continue;
        end

        % ---------------- GLOBAL projection ----------------
        if hasGlobal
            if size(A_global,2) ~= K
                if opt.Verbose
                    warning('projectGlmRezC_toAnchors:GlobalKMismatch', ...
                        '[%s] Global axes K=%d but session Y has K=%d. Skipping global projection.', ...
                        char(string(hdr)), size(A_global,2), K);
                end
            else
                Zg = project_time_major(Yproj, A_global, N, nW, opt.Eps);
                if ~isempty(Zg)
                    prj.global.ZC{j,s}   = Zg;    % [N x nW x nAxisG]
                    prj.global.okMat(j,s)= true;
                    nOKg = nOKg + 1;
                end
            end
        end

        % ---------------- PER-MOUSE projection ----------------
        [Apm, namesPm, okPm] = fetch_perMouse_anchor_axes(glmA, mouseId);
        if okPm
            if size(Apm,2) ~= K
                if opt.Verbose
                    warning('projectGlmRezC_toAnchors:PerMouseKMismatch', ...
                        '[%s] Mouse=%s axes K=%d but session Y has K=%d. Skipping per-mouse projection.', ...
                        char(string(hdr)), mouseId, size(Apm,2), K);
                end
            else
                Zm = project_time_major(Yproj, Apm, N, nW, opt.Eps);
                if ~isempty(Zm)
                    prj.perMouse.ZC{j,s}     = Zm;
                    prj.perMouse.namesC{j,s} = namesPm;
                    prj.perMouse.okMat(j,s)  = true;
                    nOKm = nOKm + 1;
                end
            end
        end

        % ---------------- PER-SESSION projection (NEW) ----------------
        [Asess, namesSess, ArawSess, okSess] = fetch_perSession_axes(glmA, j, s);
        if okSess
            if size(Asess,2) ~= K
                if opt.Verbose
                    warning('projectGlmRezC_toAnchors:PerSessionKMismatch', ...
                        '[%s] Per-session axes K=%d but session Y has K=%d. Skipping per-session projection.', ...
                        char(string(hdr)), size(Asess,2), K);
                end
            else
                Zs = project_time_major(Yproj, Asess, N, nW, opt.Eps);
                if ~isempty(Zs)
                    prj.perSession.ZC{j,s}     = Zs;
                    prj.perSession.namesC{j,s} = namesSess;
                    prj.perSession.ArawC{j,s}  = ArawSess;
                    prj.perSession.okMat(j,s)  = true;
                    nOKs = nOKs + 1;
                end
            end
        end

    end
end

if opt.Verbose
    fprintf('[projectGlmRezC_toAnchors] Global OK: %d | PerMouse OK: %d | PerSession OK: %d | Global axes=%d\n', ...
        nOKg, nOKm, nOKs, nAxisG);
end

end

%% ============================== HELPERS ==============================

function mouseId = infer_mouse_id_from_header(hdr)
hdrS = char(string(hdr));
tok = regexp(hdrS, '(m\d{3,5})', 'tokens', 'once');
if isempty(tok)
    mouseId = '';
else
    mouseId = char(string(tok{1}));
end
end

function [Yproj, nW, N, K] = get_project_target_Y(glmRez, whichY)
assert(isfield(glmRez,'decBins') && isfield(glmRez.decBins,'time') && ~isempty(glmRez.decBins.time), ...
    'glmRez.decBins.time missing; cannot infer nW.');
nW = numel(glmRez.decBins.time);

whichY = lower(string(whichY));
switch whichY
    case "yz"
        assert(isfield(glmRez,'Yz') && ~isempty(glmRez.Yz), 'glmRez.Yz missing/empty.');
        Yproj = glmRez.Yz;
    case "ybig"
        assert(isfield(glmRez,'Ybig') && ~isempty(glmRez.Ybig), 'glmRez.Ybig missing/empty.');
        Yproj = glmRez.Ybig;
    otherwise
        error('Unknown ProjectWhichY: %s', string(whichY));
end

assert(isnumeric(Yproj) && ndims(Yproj)==2, 'Yproj must be [M x K] numeric.');
[M,K] = size(Yproj);

remVal = rem(M, nW);
assert(remVal==0, 'M=%d not divisible by nW=%d (remainder=%d).', M, nW, remVal);
N = M / nW;
end

function Z = project_time_major(Yproj, A, N, nW, epsVal)
% Project rows (M x K) onto axes (nAxis x K) and reshape to [N x nW x nAxis]
if isempty(Yproj) || isempty(A)
    Z = [];
    return;
end

% normalize A rows defensively
nr = vecnorm(A,2,2);
nr(nr < epsVal) = 1;
A2 = A ./ nr;

Zrows = Yproj * A2';  % [M x nAxis]
nAxis = size(A2,1);
Z = reshape(Zrows, [N, nW, nAxis]);
end

function [A, names, ok] = fetch_perMouse_anchor_axes(glmA, mouseId)
A = [];
names = {};
ok = false;

if ~isfield(glmA,'perMouse') || ~isstruct(glmA.perMouse)
    return;
end
pm = glmA.perMouse;

if ~isfield(pm,'mouseIds') || ~isfield(pm,'axesByMouse') || isempty(mouseId)
    return;
end

mouseIds = cellstr(string(pm.mouseIds(:)));
idx = find(strcmp(mouseIds, char(string(mouseId))), 1, 'first');
if isempty(idx)
    return;
end

ax = pm.axesByMouse{idx};
if isempty(ax) || ~isstruct(ax) || ~isfield(ax,'A') || ~isfield(ax,'names') || isempty(ax.A)
    return;
end

A = ax.A;
names = cellstr(string(ax.names));
ok = true;
end

function [A, names, Araw, ok] = fetch_perSession_axes(glmA, j, s)
% Fetch per-session orthonormalized axes A for session (j,s)
A = [];
names = {};
Araw = [];
ok = false;

if ~isfield(glmA,'allSessions') || ~isstruct(glmA.allSessions) || ...
   ~isfield(glmA.allSessions,'perSessionMat')
    return;
end

perMat = glmA.allSessions.perSessionMat;
if ~iscell(perMat) || j<1 || s<1 || j>size(perMat,1) || s>size(perMat,2)
    return;
end

sess = perMat{j,s};
if isempty(sess) || ~isstruct(sess)
    return;
end

% prefer orthonormalized A if present; otherwise fall back to Araw_ord
if isfield(sess,'A') && ~isempty(sess.A)
    A = sess.A;
    if isfield(sess,'names') && ~isempty(sess.names)
        names = cellstr(string(sess.names));
    elseif isfield(sess,'names_ord') && isfield(sess,'keepAfterGS') && ~isempty(sess.keepAfterGS)
        names = cellstr(string(sess.names_ord(sess.keepAfterGS)));
    elseif isfield(sess,'names_ord')
        names = cellstr(string(sess.names_ord));
    end
else
    if isfield(sess,'Araw_ord') && ~isempty(sess.Araw_ord)
        A = sess.Araw_ord;
        if isfield(sess,'names_ord') && ~isempty(sess.names_ord)
            names = cellstr(string(sess.names_ord));
        end
    else
        return;
    end
end

if isfield(sess,'Araw_ord') && ~isempty(sess.Araw_ord)
    Araw = sess.Araw_ord;
end

ok = true;
end