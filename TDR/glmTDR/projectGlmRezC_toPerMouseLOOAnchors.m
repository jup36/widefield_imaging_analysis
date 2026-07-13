function prj = projectGlmRezC_toPerMouseLOOAnchors(anchorLOO, headerC, glmRezC, trIdC, expertHeaders_perMouse, varargin)
%projectGlmRezC_toPerMouseLOOAnchors
%
% Project each session's GLM motif activity onto per-mouse LOO-safe anchors.
%
% This is the new-pipeline replacement for the old:
%
%   prj_glmA = projectGlmRezC_toAnchors(glmA, headerC, glmRezC, ...)
%
% but uses:
%
%   anchorLOO.perMouse.axesFull{m}
%   anchorLOO.perMouse.axesLOO{m,kk}
%
% from buildPerMouseLOOAnchors_fromExperts.
%
% -------------------------------------------------------------------------
% REQUIRED INPUTS
% -------------------------------------------------------------------------
%   anchorLOO
%       Output of buildPerMouseLOOAnchors_fromExperts.
%
%   headerC
%       Cell array [mouse/session rows x session columns] of session headers.
%       Example: headerC{j,s} = 'm1045_122424'
%
%   glmRezC
%       Cell array same size as headerC.
%       Each non-empty cell should contain glmRez with fields:
%           .Yz
%           .decBins.time
%
%   trIdC
%       Cell array same size as headerC.
%       Included for bookkeeping. Not required for projection itself.
%
%   expertHeaders_perMouse
%       Mouse x expert-session cell matrix, same object used to build anchorLOO.
%       This is used to map expert header -> LOO column index kk.
%
% -------------------------------------------------------------------------
% NAME-VALUE OPTIONS
% -------------------------------------------------------------------------
%   'WhichAxes'
%       "A" default, or "Araw_ord".
%       Usually use "A", the final GS-orthogonalized axes.
%
%   'ProjectWhichY'
%       "Yz" default.
%       Field in glmRez to project.
%
%   'UseLOOForExpertSessions'
%       true default.
%       If true, expert sessions are projected onto their corresponding
%       leave-one-out axes when available.
%
%   'StoreAsSingle'
%       true default.
%       Store projected trajectories as single to reduce memory.
%
%   'Verbose'
%       true default.
%
%   'VerboseEvery'
%       25 default.
%
% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
%   prj.perMouse.ZC{j,s}
%       Projected trajectories: [nTrial x nTime x nAxis]
%
%   prj.perMouse.namesC{j,s}
%       Axis names for that session.
%
%   prj.perMouse.usedLOO(j,s)
%       True if this session was projected using axesLOO.
%
%   prj.perMouse.mouseIdxC{j,s}
%       Row index into anchorLOO.perMouse.
%
%   prj.perMouse.sourceC{j,s}
%       "axesFull", "axesLOO", or reason for empty.
%
%   prj.perMouse.timeC{j,s}
%       Time vector from glmRez.decBins.time.
%
%   prj.perMouse.headerC
%       Copy of input headerC.
%
%   prj.opt
%       Options used.
%
% -------------------------------------------------------------------------
% NOTES
% -------------------------------------------------------------------------
%   - This function intentionally does NOT compute statistics.
%   - It stores full trial x time x axis projections for plotting.
%   - It follows the same projection convention as
%     computePerMouseLOOProjectionAndStats:
%
%         Zrows = Y * Ause';
%         Z     = reshape(Zrows, [N, nW, nAxis]);
%
% -------------------------------------------------------------------------

%% -------------------- parse inputs --------------------
assert(isstruct(anchorLOO) && isfield(anchorLOO, 'perMouse'), ...
    'anchorLOO must be a struct with field .perMouse.');

assert(iscell(headerC) && iscell(glmRezC), ...
    'headerC and glmRezC must be cell arrays.');

assert(isequal(size(headerC), size(glmRezC)), ...
    'headerC and glmRezC must be the same size.');

if nargin < 4 || isempty(trIdC)
    trIdC = cell(size(headerC));
end

assert(iscell(trIdC) && isequal(size(trIdC), size(headerC)), ...
    'trIdC must be a cell array the same size as headerC.');

assert(iscell(expertHeaders_perMouse) || isstring(expertHeaders_perMouse), ...
    'expertHeaders_perMouse must be a cell or string array.');

% Defaults
opt = struct();
opt.WhichAxes               = "A";
opt.ProjectWhichY           = "Yz";
opt.UseLOOForExpertSessions = true;
opt.StoreAsSingle           = true;
opt.Verbose                 = true;
opt.VerboseEvery            = 25;

% Lightweight name-value parsing
if ~isempty(varargin)
    assert(mod(numel(varargin),2)==0, 'Name-value arguments must come in pairs.');
    for i = 1:2:numel(varargin)
        k = char(string(varargin{i}));
        v = varargin{i+1};
        assert(isfield(opt, k), 'Unknown name-value option: %s', k);
        opt.(k) = v;
    end
end

opt.WhichAxes     = string(opt.WhichAxes);
opt.ProjectWhichY = string(opt.ProjectWhichY);

%% -------------------- normalize expert header matrix --------------------
expertHdrMat = normalizeHeaderMatrix_prj_(expertHeaders_perMouse);
[nExpertMouse, nExpertMax] = size(expertHdrMat);

% Mouse IDs for rows of expertHdrMat
expertMouseIds = cell(nExpertMouse,1);
for m = 1:nExpertMouse
    expertMouseIds{m} = inferMouseIdFromRow_prj_(expertHdrMat(m,:));
end

%% -------------------- get anchor mouse IDs --------------------
anchorMouseIds = getAnchorMouseIds_prj_(anchorLOO, expertHdrMat);

nAnchorMouse = numel(anchorMouseIds);

%% -------------------- preallocate output --------------------
prj = struct();
prj.opt = opt;

prj.perMouse = struct();
prj.perMouse.headerC   = headerC;
prj.perMouse.ZC        = cell(size(headerC));
prj.perMouse.namesC    = cell(size(headerC));
prj.perMouse.timeC     = cell(size(headerC));
prj.perMouse.usedLOO   = false(size(headerC));
prj.perMouse.mouseIdxC = cell(size(headerC));
prj.perMouse.sourceC   = cell(size(headerC));
prj.perMouse.nTrial    = nan(size(headerC));
prj.perMouse.nTime     = nan(size(headerC));
prj.perMouse.nAxis     = nan(size(headerC));

prj.perMouse.anchorMouseIds = anchorMouseIds;

%% -------------------- projection loop --------------------
nTotal = numel(headerC);
nDone = 0;
nProjected = 0;
nSkipped = 0;
nUsedLOO = 0;

for lin = 1:nTotal

    [j, s] = ind2sub(size(headerC), lin);
    hdr = headerC{j,s};

    if isempty(hdr)
        continue;
    end

    hdrS = string(hdr);
    gr = glmRezC{j,s};

    nDone = nDone + 1;

    if isempty(gr) || ~isstruct(gr)
        prj.perMouse.sourceC{j,s} = 'missing_glmRez';
        nSkipped = nSkipped + 1;
        continue;
    end

    if ~isfield(gr, char(opt.ProjectWhichY)) || isempty(gr.(char(opt.ProjectWhichY)))
        prj.perMouse.sourceC{j,s} = sprintf('missing_%s', char(opt.ProjectWhichY));
        nSkipped = nSkipped + 1;
        continue;
    end

    if ~isfield(gr, 'decBins') || ~isfield(gr.decBins, 'time') || isempty(gr.decBins.time)
        prj.perMouse.sourceC{j,s} = 'missing_decBins_time';
        nSkipped = nSkipped + 1;
        continue;
    end

    mouseId = inferMouseIdFromHeader_prj_(hdrS);
    if strlength(mouseId)==0
        prj.perMouse.sourceC{j,s} = 'could_not_infer_mouseId';
        nSkipped = nSkipped + 1;
        continue;
    end

    % Find this mouse in anchorLOO
    mAnchor = find(strcmpi(string(anchorMouseIds), mouseId), 1, 'first');
    if isempty(mAnchor)
        prj.perMouse.sourceC{j,s} = sprintf('mouse_not_in_anchorLOO_%s', char(mouseId));
        nSkipped = nSkipped + 1;
        continue;
    end

    % Get full axes
    axesFull = [];
    if isfield(anchorLOO.perMouse, 'axesFull') && numel(anchorLOO.perMouse.axesFull) >= mAnchor
        axesFull = anchorLOO.perMouse.axesFull{mAnchor};
    end

    [Afull, namesFull] = getAxisMatrix_prj_(axesFull, opt.WhichAxes);

    if isempty(Afull)
        prj.perMouse.sourceC{j,s} = 'empty_axesFull';
        nSkipped = nSkipped + 1;
        continue;
    end

    % Default: full per-mouse axes
    Ause = Afull;
    namesUse = namesFull;
    usedLOO = false;
    sourceStr = 'axesFull';

    % Try LOO if this session is one of the expert sessions
    if opt.UseLOOForExpertSessions

        kk = findExpertColumnForHeader_prj_(hdrS, expertHdrMat, expertMouseIds, mouseId);

        if ~isempty(kk)
            axesLOO = [];
            if isfield(anchorLOO.perMouse, 'axesLOO')
                if size(anchorLOO.perMouse.axesLOO,1) >= mAnchor && ...
                        size(anchorLOO.perMouse.axesLOO,2) >= kk
                    axesLOO = anchorLOO.perMouse.axesLOO{mAnchor, kk};
                end
            end

            [Aloo, namesLoo] = getAxisMatrix_prj_(axesLOO, opt.WhichAxes);

            if ~isempty(Aloo)
                Ause = Aloo;
                namesUse = namesLoo;
                usedLOO = true;
                sourceStr = sprintf('axesLOO_col%d', kk);
            else
                sourceStr = sprintf('axesFull_LOO_missing_col%d', kk);
            end
        end
    end

    % Project
    Y = gr.(char(opt.ProjectWhichY));    % [M x Kmotif]
    timeVec = gr.decBins.time(:);
    nW = numel(timeVec);
    M = size(Y,1);

    if rem(M, nW) ~= 0
        prj.perMouse.sourceC{j,s} = sprintf('bad_rows_M_%d_nW_%d', M, nW);
        nSkipped = nSkipped + 1;
        continue;
    end

    if size(Y,2) ~= size(Ause,2)
        prj.perMouse.sourceC{j,s} = sprintf( ...
            'dimension_mismatch_Ycols_%d_Acols_%d', size(Y,2), size(Ause,2));
        nSkipped = nSkipped + 1;
        continue;
    end

    N = M / nW;
    nAxis = size(Ause,1);

    Zrows = Y * Ause';                  % [M x nAxis]
    Z = reshape(Zrows, [N, nW, nAxis]); % [trial x time x axis]

    if opt.StoreAsSingle
        Z = single(Z);
    end

    % Store
    prj.perMouse.ZC{j,s}        = Z;
    prj.perMouse.namesC{j,s}    = namesUse;
    prj.perMouse.timeC{j,s}     = timeVec;
    prj.perMouse.usedLOO(j,s)   = usedLOO;
    prj.perMouse.mouseIdxC{j,s} = mAnchor;
    prj.perMouse.sourceC{j,s}   = sourceStr;
    prj.perMouse.nTrial(j,s)    = N;
    prj.perMouse.nTime(j,s)     = nW;
    prj.perMouse.nAxis(j,s)     = nAxis;

    nProjected = nProjected + 1;
    if usedLOO
        nUsedLOO = nUsedLOO + 1;
    end

    if opt.Verbose && (nDone==1 || mod(nDone,opt.VerboseEvery)==0)
        fprintf('[projectLOO] %d/%d checked | projected=%d | LOO=%d | skipped=%d | latest=%s\n', ...
            nDone, nTotal, nProjected, nUsedLOO, nSkipped, char(hdrS));
    end
end

%% -------------------- summary --------------------
prj.summary = struct();
prj.summary.nCheckedNonEmpty = nDone;
prj.summary.nProjected       = nProjected;
prj.summary.nUsedLOO         = nUsedLOO;
prj.summary.nSkipped         = nSkipped;

if opt.Verbose
    fprintf('[projectLOO] done | checked=%d | projected=%d | LOO=%d | skipped=%d\n', ...
        nDone, nProjected, nUsedLOO, nSkipped);
end

end

%% ========================================================================
% Local helpers
% ========================================================================

function hdrMat = normalizeHeaderMatrix_prj_(hdrIn)

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
    elseif ischar(x) || isstring(x)
        s = char(string(x));
        if strlength(string(s)) == 0
            hdrMat{i} = [];
        else
            hdrMat{i} = s;
        end
    else
        hdrMat{i} = [];
    end
end

end

function mouseIds = getAnchorMouseIds_prj_(anchorLOO, expertHdrMat)

mouseIds = {};

if isfield(anchorLOO.perMouse, 'mouseIds') && ~isempty(anchorLOO.perMouse.mouseIds)
    mouseIds = anchorLOO.perMouse.mouseIds(:);
    mouseIds = cellfun(@char, cellstr(string(mouseIds)), 'UniformOutput', false);
    return;
end

% Fallback: infer from expert header rows
nM = size(expertHdrMat,1);
mouseIds = cell(nM,1);
for m = 1:nM
    mouseIds{m} = inferMouseIdFromRow_prj_(expertHdrMat(m,:));
end

end

function mouseId = inferMouseIdFromRow_prj_(hdrRow)

mouseId = '';

for k = 1:numel(hdrRow)
    h = hdrRow{k};
    if isempty(h)
        continue;
    end

    tok = regexp(string(h), '(m\d{3,5})', 'tokens', 'once');
    if ~isempty(tok)
        mouseId = char(string(tok{1}));
        return;
    end
end

end

function mouseId = inferMouseIdFromHeader_prj_(hdr)

tok = regexp(string(hdr), '(m\d{3,5})', 'tokens', 'once');

if isempty(tok)
    mouseId = "";
else
    mouseId = string(tok{1});
end

end

function kk = findExpertColumnForHeader_prj_(hdr, expertHdrMat, expertMouseIds, mouseId)

kk = [];

if isempty(expertHdrMat)
    return;
end

mRow = find(strcmpi(string(expertMouseIds), string(mouseId)), 1, 'first');

if isempty(mRow)
    return;
end

row = string(expertHdrMat(mRow,:));
hit = find(row == string(hdr), 1, 'first');

if ~isempty(hit)
    kk = hit;
end

end

function [A, names] = getAxisMatrix_prj_(axesStruct, whichAxes)

A = [];
names = {};

if isempty(axesStruct) || ~isstruct(axesStruct)
    return;
end

whichAxes = string(whichAxes);

switch whichAxes
    case "A"
        if isfield(axesStruct, 'A') && ~isempty(axesStruct.A)
            A = axesStruct.A;
        end

        if isfield(axesStruct, 'names') && ~isempty(axesStruct.names)
            names = axesStruct.names;
        elseif isfield(axesStruct, 'names_ord') && ~isempty(axesStruct.names_ord)
            names = axesStruct.names_ord;
        end

    case "Araw_ord"
        if isfield(axesStruct, 'Araw_ord') && ~isempty(axesStruct.Araw_ord)
            A = axesStruct.Araw_ord;
        end

        if isfield(axesStruct, 'names_ord') && ~isempty(axesStruct.names_ord)
            names = axesStruct.names_ord;
        elseif isfield(axesStruct, 'names') && ~isempty(axesStruct.names)
            names = axesStruct.names;
        end

    otherwise
        error('WhichAxes must be "A" or "Araw_ord".');
end

if ~isempty(A) && isempty(names)
    names = arrayfun(@(i) sprintf('axis_%d', i), 1:size(A,1), 'UniformOutput', false);
end

% Ensure row cellstr
if isstring(names)
    names = cellstr(names);
end
names = names(:)';

end