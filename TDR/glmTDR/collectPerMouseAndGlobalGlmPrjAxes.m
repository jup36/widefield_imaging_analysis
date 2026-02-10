function anchor = collectPerMouseAndGlobalGlmPrjAxes(headerC, glmRezC, varargin)
%collectPerMouseAndGlobalGlmPrjAxes  Match expert header lists to glmRezC and build anchor axes.
%
% anchor = collectPerMouseAndGlobalGlmPrjAxes(headerC, glmRezC, 'Name', value, ...)
%
% REQUIRED INPUTS
%   headerC : cell array (J x S) of session headers (e.g., 'm1045_122424')
%   glmRezC : cell array (J x S) of glmRez structs (same size as headerC)
%
% NAME-VALUE
%   'expertHeaders_perMouse' : cell array (nMouse x nSess) of headers (empties allowed)
%   'expertHeaders_global'   : cell array (nMouse x nSess) of headers (empties allowed)
%
%   % Passed through to buildAnchorAxes_from_glmRezList:
%   'OrthMode'        : "GS" (default) | "none"
%   'SignFix'         : "maxabs" (default) | "none"
%   'Eps'             : 1e-10
%   'MultiDimGroups'  : {'GoToneOn','NoGoToneOn','ToneOffGo','ToneOffNoGo'}
%   'MultiDimK'       : 3
%   'PriorityNames'   : {'GoToneOn','NoGoToneOn','ToneOffGo','ToneOffNoGo','Lick'}
%   'Verbose'         : true
%
% OUTPUT (anchor struct)
%   .perMouse.headersMat      : original headers matrix (nMouse x nSess)
%   .perMouse.mouseIds        : cellstr (nMouse x 1)
%   .perMouse.foundMat        : logical (nMouse x nSess)
%   .perMouse.glmRezMat       : cell (nMouse x nSess) matched glmRez or []
%   .perMouse.axesByMouse     : cell (nMouse x 1) each = pooled anchor axes struct (or [])
%   .perMouse.headersUsed     : cell (nMouse x 1) headers actually pooled (found only)
%
%   .global.headersMat        : original headers matrix
%   .global.foundMat          : logical matrix
%   .global.glmRezMat         : cell matrix matched glmRez or []
%   .global.axes              : pooled anchor axes struct across all found global headers
%   .global.headersUsed       : cellstr vector of headers pooled (found only)
%
% Notes
% - Matching is exact string match on header.
% - If duplicates exist in headerC, we warn and use first occurrence.
% - Empty entries ([], '', 0x0 double) inside expertHeaders_* are ignored.

% -------------------- parse --------------------
p = inputParser;
p.addRequired('headerC', @(c) iscell(c));
p.addRequired('glmRezC', @(c) iscell(c) && isequal(size(c), size(headerC)));

p.addParameter('expertHeaders_perMouse', {}, @(c) iscell(c) || isstring(c));
p.addParameter('expertHeaders_global',   {}, @(c) iscell(c) || isstring(c));

% pass-through options for buildAnchorAxes_from_glmRezList
p.addParameter('OrthMode',"GS",@(s)ischar(s)||isstring(s));
p.addParameter('SignFix',"maxabs",@(s)ischar(s)||isstring(s));
p.addParameter('Eps',1e-10,@(x)isscalar(x)&&x>0);
p.addParameter('MultiDimGroups', {'GoToneOn','NoGoToneOn','ToneOffGo','ToneOffNoGo'}, @(c)iscell(c)||isstring(c));
p.addParameter('MultiDimK', 3, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('PriorityNames', {'GoToneOn','NoGoToneOn','ToneOffGo','ToneOffNoGo','Lick'}, @(c)iscell(c)||isstring(c));
p.addParameter('Verbose', true, @(x)islogical(x)&&isscalar(x));

p.parse(headerC, glmRezC, varargin{:});
opt = p.Results;

% -------------------- normalize expert header matrices --------------------
hdrPM_mat = normalizeHeaderMatrix_(opt.expertHeaders_perMouse);
hdrG_mat  = normalizeHeaderMatrix_(opt.expertHeaders_global);

% -------------------- build lookup map from headerC -> glmRezC --------------------
[uniqHdr, idxFirst, glmFlatFirst] = buildHeaderLookup_(headerC, glmRezC);

% -------------------- per-mouse pooling --------------------
[nMousePM, nSessPM] = size(hdrPM_mat);
pm_foundMat  = false(nMousePM, nSessPM);
pm_glmMat    = cell(nMousePM, nSessPM);
pm_mouseIds  = cell(nMousePM, 1);
pm_axesByMouse = cell(nMousePM, 1);
pm_headersUsed = cell(nMousePM, 1);

for iM = 1:nMousePM
    hdrRow = hdrPM_mat(iM,:);

    % infer mouse id from first non-empty header in row
    pm_mouseIds{iM} = inferMouseIdFromRow_(hdrRow);

    glmPool = {};
    usedHdr = strings(0,1);

    for iS = 1:nSessPM
        h = hdrRow{iS};
        if isempty(h), continue; end

        [g, ok] = fetchByHeaderFromLookup_(h, uniqHdr, idxFirst, glmFlatFirst);

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
        'OrthMode', opt.OrthMode, ...
        'SignFix',  opt.SignFix, ...
        'Eps',      opt.Eps, ...
        'MultiDimGroups', cellstr(string(opt.MultiDimGroups)), ...
        'MultiDimK', opt.MultiDimK, ...
        'PriorityNames', cellstr(string(opt.PriorityNames)), ...
        'Verbose',  opt.Verbose);
end

% -------------------- global pooling --------------------
[nMouseG, nSessG] = size(hdrG_mat);
g_foundMat = false(nMouseG, nSessG);
g_glmMat   = cell(nMouseG, nSessG);

glmPoolGlobal = {};
hdrUsedGlobal = strings(0,1);

for iM = 1:nMouseG
    for iS = 1:nSessG
        h = hdrG_mat{iM,iS};
        if isempty(h), continue; end

        [g, ok] = fetchByHeaderFromLookup_(h, uniqHdr, idxFirst, glmFlatFirst);

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
        'OrthMode', opt.OrthMode, ...
        'SignFix',  opt.SignFix, ...
        'Eps',      opt.Eps, ...
        'MultiDimGroups', cellstr(string(opt.MultiDimGroups)), ...
        'MultiDimK', opt.MultiDimK, ...
        'PriorityNames', cellstr(string(opt.PriorityNames)), ...
        'Verbose',  opt.Verbose);
end

% -------------------- pack output --------------------
anchor = struct();

anchor.perMouse = struct();
anchor.perMouse.headersMat   = hdrPM_mat;
anchor.perMouse.mouseIds     = pm_mouseIds;
anchor.perMouse.foundMat     = pm_foundMat;
anchor.perMouse.glmRezMat    = pm_glmMat;
anchor.perMouse.axesByMouse  = pm_axesByMouse;
anchor.perMouse.headersUsed  = pm_headersUsed;

anchor.global = struct();
anchor.global.headersMat     = hdrG_mat;
anchor.global.foundMat       = g_foundMat;
anchor.global.glmRezMat      = g_glmMat;
anchor.global.axes           = globalAxes;
anchor.global.headersUsed    = cellstr(hdrUsedGlobal);

anchor.opt = opt;

if opt.Verbose
    fprintf('[anchorMatch] perMouse pooled=%d/%d mice | global pooled=%d sessions\n', ...
        sum(~cellfun(@isempty, pm_axesByMouse)), nMousePM, numel(glmPoolGlobal));
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% HELPERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hdrMat = normalizeHeaderMatrix_(hdrIn)
% Convert input (cell or string array) to cell matrix of char headers; empty -> [].

if isempty(hdrIn)
    hdrMat = cell(0,0);
    return;
end

if isstring(hdrIn)
    hdrIn = cellstr(hdrIn);
end

% If user passed a vector, make it a column.
if iscell(hdrIn) && isvector(hdrIn)
    hdrIn = hdrIn(:);
end

hdrMat = cell(size(hdrIn));
for i = 1:numel(hdrIn)
    x = hdrIn{i};

    % common "empty cell shows as 0x0 double" case
    if isempty(x)
        hdrMat{i} = [];
        continue;
    end

    if isstring(x) || ischar(x)
        s = char(string(x));
        if strlength(string(s))==0
            hdrMat{i} = [];
        else
            hdrMat{i} = s;
        end
    else
        % anything else (e.g. numeric empty) -> treat as empty
        hdrMat{i} = [];
    end
end
end

function [uniqHdr, idxFirst, glmFlatFirst] = buildHeaderLookup_(headerC, glmRezC)
% Flatten headerC/glmRezC and build unique header list + first index mapping.

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

idxFirst = zeros(numel(uniqHdr),1);
glmFlatFirst = cell(numel(uniqHdr),1);
for i = 1:numel(uniqHdr)
    ii = find(hdrFlatS == uniqHdr(i), 1, 'first');
    idxFirst(i) = ii;
    glmFlatFirst{i} = glmFlat{ii};
end
end

function [glmRez, found] = fetchByHeaderFromLookup_(hdr, uniqHdr, idxFirst, glmFlatFirst)
% Exact match fetch. Returns glmRez struct or [].

hdr = string(hdr);
j = find(uniqHdr == hdr, 1, 'first');
if isempty(j)
    glmRez = [];
    found = false;
    return;
end

% idxFirst is kept for debugging/consistency, but glmFlatFirst is already "first"
glmRez = glmFlatFirst{j};

found = ~isempty(glmRez) && isstruct(glmRez);
end

function mouseId = inferMouseIdFromRow_(hdrRow)
% Get mouse ID from first non-empty header in a row.

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
