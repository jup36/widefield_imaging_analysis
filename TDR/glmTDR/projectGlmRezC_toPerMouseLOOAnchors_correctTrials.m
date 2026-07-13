function prj = projectGlmRezC_toPerMouseLOOAnchors_correctTrials(anchorLOO, headerC, glmRezC, trIdC, expertHeaders_perMouse, varargin)
%projectGlmRezC_toPerMouseLOOAnchors_correctTrials
%
% Project each session's GLM motif activity onto per-mouse LOO-safe anchors,
% MASKING NON-CORRECT TRIALS TO NaN ON A PER-AXIS BASIS BY DEFAULT.
%
% This is a drop-in variant of projectGlmRezC_toPerMouseLOOAnchors. The
% projection math (Zrows = Y * Ause'; reshape to [N, nW, nAxis]) is
% IDENTICAL. The only addition is: after computing Z, trials that do not
% match the "correct" outcome appropriate to a given axis's predictor
% group are set to NaN for that axis (all time bins), so that any
% downstream trial-averaging (e.g. mean(...,'omitnan'), Spearman/slope
% stats computed with 'omitnan') automatically excludes them without
% requiring changes to the stats code itself.
%
% WHY PER-AXIS NaN-MASKING INSTEAD OF DROPPING TRIALS
%   Different axes correspond to different predictor groups (e.g.
%   "GoToneOn_1..3" vs "NoGoToneOn_1..3"), and the "correct" trial
%   definition differs by group (Hit for Go-related groups, CR for
%   NoGo-related groups). Since Z is a single [nTrial x nTime x nAxis]
%   array shared across all axes, trials can't be dropped per-axis
%   without breaking the shared trial dimension. NaN-masking keeps the
%   array shape identical to the original function's output while making
%   incorrect trials invisible to nan-aware downstream reductions.
%
% -------------------------------------------------------------------------
% REQUIRED INPUTS  (identical to projectGlmRezC_toPerMouseLOOAnchors)
% -------------------------------------------------------------------------
%   anchorLOO               : output of buildPerMouseLOOAnchors_fromExperts
%                              (or buildPerMouseLOOAnchors_fromExperts_correctTrials)
%   headerC                 : cell [rows x cols] of session headers
%   glmRezC                 : cell, same size as headerC, glmRez structs
%   trIdC                   : cell, same size as headerC, trial-type structs
%                              with fields goI/nogoI/hitI/missI/crI/faI
%                              (REQUIRED here, unlike the original function,
%                              since it is now actually used for masking)
%   expertHeaders_perMouse  : mouse x expert-session header cell matrix
%
% -------------------------------------------------------------------------
% NAME-VALUE OPTIONS
% -------------------------------------------------------------------------
%   'WhichAxes'                : "A" (default) | "Araw_ord"
%   'ProjectWhichY'             : "Yz" (default)
%   'UseLOOForExpertSessions'   : true (default)
%   'StoreAsSingle'             : true (default)
%   'Verbose'                   : true (default)
%   'VerboseEvery'              : 25 (default)
%
%   NEW:
%   'TrialPolicy'               : "CorrectOnlyByGroup" (default)
%                                  Determines which trials are considered
%                                  "valid" for a given axis, based on the
%                                  axis's predictor-group name:
%                                    "CorrectOnlyByGroup" -> Hit trials for
%                                       Go-related groups (e.g. GoToneOn,
%                                       ToneOffGo), CR trials for NoGo-
%                                       related groups (e.g. NoGoToneOn,
%                                       ToneOffNoGo).
%                                    "GoNoGoByGroup" -> Go trials / NoGo
%                                       trials (legacy stimulus-identity
%                                       masking; all outcomes included).
%                                    "all" -> no masking (identical
%                                       behavior to the original
%                                       projectGlmRezC_toPerMouseLOOAnchors).
%                                  Axes whose group name does not match
%                                  Go/NoGo-related patterns (e.g. "Lick")
%                                  are always left unmasked ("all trials"),
%                                  regardless of policy.
%
%   'MaskIncorrectTrials'       : true (default)
%                                  If false, Z is NOT masked (kept
%                                  identical to legacy behavior) but the
%                                  trial-validity mask is still computed
%                                  and returned, so you can inspect it or
%                                  apply it downstream yourself.
%
% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
%   prj.perMouse.ZC{j,s}
%       Projected trajectories: [nTrial x nTime x nAxis].
%       Entries for (trial, :, axis) are NaN wherever that trial is not
%       "valid" for that axis's group under TrialPolicy (unless
%       MaskIncorrectTrials=false).
%
%   prj.perMouse.namesC{j,s}          : axis names (as before)
%   prj.perMouse.usedLOO(j,s)         : as before
%   prj.perMouse.mouseIdxC{j,s}       : as before
%   prj.perMouse.sourceC{j,s}         : as before
%   prj.perMouse.timeC{j,s}           : as before
%   prj.perMouse.headerC              : as before
%
%   NEW:
%   prj.perMouse.trialValidC{j,s}
%       Logical [nTrial x nAxis]. true = trial counted as "correct"/valid
%       for that axis's group under TrialPolicy. Always computed, even if
%       MaskIncorrectTrials=false, so you can apply/inspect masking
%       independently of what was actually applied to ZC.
%
%   prj.perMouse.nValidC{j,s}
%       [1 x nAxis] double. Count of valid trials per axis, i.e.
%       sum(trialValidC{j,s},1). Useful for spotting sessions/axes with
%       too few correct trials to trust the trajectory estimate.
%
%   prj.opt
%       Options used (includes TrialPolicy, MaskIncorrectTrials).
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

% NEW
opt.TrialPolicy             = "CorrectOnlyByGroup";
opt.MaskIncorrectTrials     = true;

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
opt.TrialPolicy   = string(opt.TrialPolicy);

%% -------------------- normalize expert header matrix --------------------
expertHdrMat = normalizeHeaderMatrix_prj_(expertHeaders_perMouse);
[nExpertMouse, ~] = size(expertHdrMat);

% Mouse IDs for rows of expertHdrMat
expertMouseIds = cell(nExpertMouse,1);
for m = 1:nExpertMouse
    expertMouseIds{m} = inferMouseIdFromRow_prj_(expertHdrMat(m,:));
end

%% -------------------- get anchor mouse IDs --------------------
anchorMouseIds = getAnchorMouseIds_prj_(anchorLOO, expertHdrMat);

%% -------------------- preallocate output --------------------
prj = struct();
prj.opt = opt;

prj.perMouse = struct();
prj.perMouse.headerC     = headerC;
prj.perMouse.ZC          = cell(size(headerC));
prj.perMouse.namesC      = cell(size(headerC));
prj.perMouse.timeC       = cell(size(headerC));
prj.perMouse.usedLOO     = false(size(headerC));
prj.perMouse.mouseIdxC   = cell(size(headerC));
prj.perMouse.sourceC     = cell(size(headerC));
prj.perMouse.nTrial      = nan(size(headerC));
prj.perMouse.nTime       = nan(size(headerC));
prj.perMouse.nAxis       = nan(size(headerC));

% NEW
prj.perMouse.trialValidC = cell(size(headerC));
prj.perMouse.nValidC     = cell(size(headerC));

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
    tr = trIdC{j,s};

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

    % ---- NEW: per-axis trial-validity mask + optional NaN-masking ----
    trialValid = true(N, nAxis);
    for a = 1:nAxis
        gName = extractGroupNameFromAxisName_(namesUse{a});
        trialValid(:,a) = buildTrialValidMaskForAxisGroup_(tr, gName, N, opt.TrialPolicy);
    end

    if opt.MaskIncorrectTrials
        for a = 1:nAxis
            invalidTrials = ~trialValid(:,a);
            if any(invalidTrials)
                Z(invalidTrials, :, a) = NaN;
            end
        end
    end
    % --------------------------------------------------------------------

    if opt.StoreAsSingle
        Z = single(Z);
    end

    % Store
    prj.perMouse.ZC{j,s}          = Z;
    prj.perMouse.namesC{j,s}      = namesUse;
    prj.perMouse.timeC{j,s}       = timeVec;
    prj.perMouse.usedLOO(j,s)     = usedLOO;
    prj.perMouse.mouseIdxC{j,s}   = mAnchor;
    prj.perMouse.sourceC{j,s}     = sourceStr;
    prj.perMouse.nTrial(j,s)      = N;
    prj.perMouse.nTime(j,s)       = nW;
    prj.perMouse.nAxis(j,s)       = nAxis;

    prj.perMouse.trialValidC{j,s} = trialValid;
    prj.perMouse.nValidC{j,s}     = sum(trialValid, 1);

    nProjected = nProjected + 1;
    if usedLOO
        nUsedLOO = nUsedLOO + 1;
    end

    if opt.Verbose && (nDone==1 || mod(nDone,opt.VerboseEvery)==0)
        fprintf('[projectLOO:correctTrials] %d/%d checked | projected=%d | LOO=%d | skipped=%d | latest=%s\n', ...
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
    fprintf('[projectLOO:correctTrials] done | checked=%d | projected=%d | LOO=%d | skipped=%d\n', ...
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

function gName = extractGroupNameFromAxisName_(axisName)
% Axis names are constructed elsewhere as "<groupName>_<pcIdx>", e.g.
% "GoToneOn_2". Strip the trailing "_<digits>" to recover the group name.
% If the pattern doesn't match (e.g. a fallback "axis_3" name with no
% known group), return the name unchanged; buildTrialValidMaskForAxisGroup_
% will then treat it as "all trials" since it won't match Go/NoGo patterns.

s = char(string(axisName));
tok = regexp(s, '^(.*)_(\d+)$', 'tokens', 'once');
if ~isempty(tok)
    gName = tok{1};
else
    gName = s;
end

end

function trialValid = buildTrialValidMaskForAxisGroup_(trId, groupName, N, policy)
% Returns logical [N x 1]: true = trial counted as valid for this axis's
% group under the given policy.
%
% Group-name matching mirrors the convention used when building the
% correct-trials-only anchors: any group name containing "nogo"
% (case-insensitive) is treated as NoGo-related; everything else that is
% a recognized Go/NoGo-family group name is treated as Go-related. Axes
% whose group name doesn't look Go/NoGo-related at all (e.g. "Lick",
% or fallback "axis_k" names) are left unmasked ("all trials"), since
% there is no well-defined "correct" trial-type restriction for them.

policy = string(policy);
groupName = string(groupName);

trialValid = true(N,1);

if isempty(trId) || ~isstruct(trId)
    return;
end

isNoGoFamily = contains(lower(groupName), "nogo");
isGoFamily   = contains(lower(groupName), "go") && ~isNoGoFamily;
% covers "GoToneOn", "ToneOffGo" (contains "go"), "NoGoToneOn",
% "ToneOffNoGo" (contains "nogo", caught above first)

if ~isGoFamily && ~isNoGoFamily
    % Not a Go/NoGo-related group (e.g. "Lick") -> no outcome-based
    % restriction; all trials valid regardless of policy.
    return;
end

want = "all";
if strcmpi(policy, "GoNoGoByGroup")
    want = ternary_(isNoGoFamily, "nogo", "go");
elseif strcmpi(policy, "CorrectOnlyByGroup")
    want = ternary_(isNoGoFamily, "cr", "hit");
end

if want == "all"
    return;
end

switch want
    case "go",   trialValid = local_makeMask_prj_(trId, 'goI',   N);
    case "nogo", trialValid = local_makeMask_prj_(trId, 'nogoI', N);
    case "hit",  trialValid = local_makeMask_prj_(trId, 'hitI',  N);
    case "cr",   trialValid = local_makeMask_prj_(trId, 'crI',   N);
end

end

function out = ternary_(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end

function m = local_makeMask_prj_(trId, field, N)
m = true(N,1);
if ~isfield(trId, field) || isempty(trId.(field))
    warning('local_makeMask_prj_:MissingField', ...
        'trId.%s missing/empty; leaving trials unmasked (all valid) for this axis.', field);
    return;
end
x = trId.(field);

if islogical(x)
    x = x(:);
    if numel(x) ~= N
        warning('local_makeMask_prj_:BadLength', ...
            'trId.%s length %d != N=%d; leaving trials unmasked (all valid) for this axis.', ...
            field, numel(x), N);
        return;
    end
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