function prj = projectGlmRezC_toAnchors(glmA, headerC, glmRezC, varargin)
%PROJECTGLMREZC_TOANCHORS  Project each session's motif activity onto per-mouse + global anchor axes.
%
% prj = projectGlmRezC_toAnchors(glmA, headerC, glmRezC, 'Name', value, ...)
%
% INPUTS
%   glmA    : anchor struct produced by collectPerMouseAndGlobalGlmPrjAxes (your current schema)
%             - global axes live at:  glmA.global.axes.A / glmA.global.axes.names
%             - per-mouse axes live at: glmA.perMouse.axesByMouse{m}.A / .names
%   headerC : (J x S) cell array of session headers (e.g. 'm1045_122424')
%   glmRezC : (J x S) cell array of glmRez structs (same size as headerC)
%
% NAME-VALUE
%   'ProjectWhichY' : "Yz" (default) | "Ybig"
%   'Eps'           : numeric (default 1e-10)
%   'Verbose'       : true/false (default true)
%
% OUTPUT (prj struct)
%   .global
%       .A, .names         : global anchor axes + names (copied from glmA.global.axes)
%       .ZC                : (J x S) cell; each = [N x nW x nAxis] projected trajectories
%       .okMat             : (J x S) logical; projection success
%   .perMouse
%       .ZC                : (J x S) cell; each = [N x nW x nAxis_mouse] using that mouse's anchor
%       .namesC            : (J x S) cell; each = axis names used for that session (typically same within mouse)
%       .okMat             : (J x S) logical
%       .mouseIdMat        : (J x S) cellstr mouseId inferred from header
%   .meta
%       .headerC           : copy of headerC
%       .ProjectWhichY     : option used
%
% IMPORTANT STACKING ASSUMPTION (matches your pipeline)
%   glmRez.Yz is TIME-MAJOR in rows:
%     [time bin 1: trials 1..N], [time bin 2: trials 1..N], ..., [time bin nW]
%   Reshape back to trial x time x feature with:
%     Z = reshape(Zrows, [N, nW, nAxis]);
%
% Junchol Park lab-style: deterministic + robust + explicit.

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
assert(isfield(glmA,'global') && isstruct(glmA.global) && ...
       isfield(glmA.global,'axes') && isstruct(glmA.global.axes) && ...
       isfield(glmA.global.axes,'A') && isfield(glmA.global.axes,'names'), ...
       'glmA.global.axes.A / glmA.global.axes.names not found. Did you run global pooling?');

A_global = glmA.global.axes.A;
names_global = glmA.global.axes.names;

assert(isnumeric(A_global) && ndims(A_global)==2, 'glmA.global.axes.A must be 2D numeric.');
assert(iscell(names_global) || isstring(names_global), 'glmA.global.axes.names must be cellstr or string array.');
names_global = cellstr(string(names_global));

nAxisG = size(A_global,1);

% -------------------- allocate outputs --------------------
prj = struct();

prj.global = struct();
prj.global.A     = A_global;
prj.global.names = names_global;
prj.global.ZC    = cell(J,S);
prj.global.okMat = false(J,S);

prj.perMouse = struct();
prj.perMouse.ZC        = cell(J,S);
prj.perMouse.namesC    = cell(J,S);
prj.perMouse.okMat     = false(J,S);
prj.perMouse.mouseIdMat= cell(J,S);

prj.meta = struct();
prj.meta.headerC = headerC;
prj.meta.ProjectWhichY = char(string(opt.ProjectWhichY));

% -------------------- main loop --------------------
nOKg = 0;
nOKm = 0;

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
        if size(A_global,2) ~= K
            if opt.Verbose
                warning('projectGlmRezC_toAnchors:GlobalKMismatch', ...
                    '[%s] Global axes K=%d but session Y has K=%d. Skipping global projection.', char(string(hdr)), size(A_global,2), K);
            end
        else
            Zg = project_time_major(Yproj, A_global, N, nW, opt.Eps);
            if ~isempty(Zg)
                prj.global.ZC{j,s}  = Zg;   % [N x nW x nAxisG]
                prj.global.okMat(j,s)= true;
                nOKg = nOKg + 1;
            end
        end

        % ---------------- PER-MOUSE projection ----------------
        [Apm, namesPm, okPm] = fetch_perMouse_anchor_axes(glmA, mouseId);
        if ~okPm
            continue;
        end

        if size(Apm,2) ~= K
            if opt.Verbose
                warning('projectGlmRezC_toAnchors:PerMouseKMismatch', ...
                    '[%s] Mouse=%s axes K=%d but session Y has K=%d. Skipping per-mouse projection.', ...
                    char(string(hdr)), mouseId, size(Apm,2), K);
            end
            continue;
        end

        Zm = project_time_major(Yproj, Apm, N, nW, opt.Eps);
        if ~isempty(Zm)
            prj.perMouse.ZC{j,s}     = Zm;      % [N x nW x nAxisMouse]
            prj.perMouse.namesC{j,s} = namesPm; % names used for this session
            prj.perMouse.okMat(j,s)  = true;
            nOKm = nOKm + 1;
        end

    end
end

if opt.Verbose
    fprintf('[projectGlmRezC_toAnchors] Global OK: %d | PerMouse OK: %d | Global axes=%d\n', ...
        nOKg, nOKm, nAxisG);
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% HELPERS (END OF FILE) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mouseId = infer_mouse_id_from_header(hdr)
% Infer m#### or m##### from header string. Returns '' if not found.
hdrS = char(string(hdr));
tok = regexp(hdrS, '(m\d{3,5})', 'tokens', 'once');
if isempty(tok)
    mouseId = '';
else
    mouseId = char(string(tok{1}));
end
end

function [Yproj, nW, N, K] = get_project_target_Y(glmRez, whichY)
% Get Yproj (M x K) from glmRez and return (nW, N, K).
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

% normalize A rows defensively (should already be normalized, but be safe)
nr = vecnorm(A,2,2);
nr(nr < epsVal) = 1;
A2 = A ./ nr;

Zrows = Yproj * A2';  % [M x nAxis]
M = size(Yproj,1);

if size(Zrows,1) ~= M
    Z = [];
    return;
end

nAxis = size(A2,1);
Z = reshape(Zrows, [N, nW, nAxis]);  % TIME-MAJOR invert
end

function [A, names, ok] = fetch_perMouse_anchor_axes(glmA, mouseId)
%FETCH_PERMOUSE_ANCHOR_AXES  Return per-mouse anchor axes for a mouseId.
%
% Expected schema (your current glmA):
%   glmA.perMouse.mouseIds    : cellstr (nMouse x 1)
%   glmA.perMouse.axesByMouse : cell (nMouse x 1), each cell is a struct:
%                               .A      [nAxis x K]
%                               .names  {1 x nAxis}

A = [];
names = {};
ok = false;

if ~isfield(glmA,'perMouse') || ~isstruct(glmA.perMouse)
    return;
end
pm = glmA.perMouse;

if ~isfield(pm,'mouseIds') || ~isfield(pm,'axesByMouse')
    return;
end

if isempty(mouseId)
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
