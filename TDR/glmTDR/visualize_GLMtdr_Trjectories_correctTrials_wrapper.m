%% =========================================================================
%  GLM-TDR: quick load + single-session trial-avg trajectories + across-session axis timecourses
%  --- CORRECT TRIALS ONLY VARIANT ---
%
%  Axes (anchorLOO_correctTrials) and collapsed trend statistics
%  (rezProjStats_correctTrials) are loaded from the correct-trials-only
%  pipeline outputs. The single-session / per-trial-type projection used
%  for visualization (prj_anchorLOO_onCorrectAxes) is deliberately
%  UNMASKED (built with the original projectGlmRezC_toPerMouseLOOAnchors,
%  not the NaN-masking _correctTrials variant) so that within-session
%  trial-type comparisons (e.g. CR vs Hit on the same axis) still show
%  both trial types. See the loading block below for details.
%
%  SYNOPSIS
%  This script demonstrates two common readouts from GLM-based targeted
%  dimensionality reduction (GLM-TDR) results saved in glmTdrRezCollection_020626.mat:
%
%   (1) Single-session visualization:
%       - Extract a session by header (e.g., m1045_122424)
%       - Pull projected trial x time x axis trajectories (global or perMouse)
%       - Plot trial-averaged trajectories for selected trial types (e.g., CR vs Hit)
%
%   (2) Across-session visualization:
%       - For a chosen mouse, plot how a single GLM axis projection (1-D) evolves
%         across sessions (e.g., GoToneOn_3), for a given trial type (e.g., hitI)
%
%  KEY VARIABLES (from the saved collection)
%    headerC   : [mouse x session] headers (e.g., 'm1045_122424')
%    trIdC     : [mouse x session] trial-type logical indices (hitI, crI, etc.)
%    prj_glmA  : projections onto GLM-TDR anchor axes
%                .global.ZC{j,s}   -> [N x T x Dglobal]
%                .perMouse.ZC{j,s} -> [N x T x Dmouse]
%                .global.names     -> 1xDglobal axis names (e.g., 'NoGoToneOn_1')
%                .perMouse.namesC{j,s} -> 1xDmouse axis names for that session/mouse
%
%  ASSUMPTIONS
%    - Helper plotting functions exist on path:
%        plotTrAvgGlmTrjs, perMouseAcrossSessionPrjScoreTrajectories
%    - You want axis indices by name (glmNameToColumns)
%    - TIME vector used for plotting is consistent across sessions for your use-case
% =========================================================================

%% whereabouts
fileDir = compatiblepath("Z:\Rodent Data\dualImaging_parkj\collectData");

%% Load saved GLM-TDR collection
% NOTE: include glmRezC here if you plan to pull timestamps from glmRezC (below).
%load(fullfile(fileDir, "glmTdrRezCollection_020626.mat"), ...
%    "glmRezPathC", "headerC", "trIdC", "glmA", "prj_glmA", "glmRezC");
%
% headerC/trIdC/glmRezC are the raw per-session GLM inputs and are the
% SAME regardless of trial policy -- no change needed here.
load(fullfile(fileDir, "glmTdrRezCollection_redCal_L10K10_070626.mat"), ...
    "glmRezPathC", "headerC", "trIdC", "glmRezC");

% anchorLOO_correctTrials    : axes built from correct trials only (Hit/CR)
% rezProjStats_correctTrials : muGo/muNoGo/diff/energy + trend stats, also
%                              computed from correct trials only
load(fullfile(fileDir, "glmTdrAxesProjRezCollection_correctTrials_redCal_L10K10_070626.mat"), ...
    "anchorLOO_correctTrials", "rezProjStats_correctTrials", "hdrC_all");

% ---------------------------------------------------------------------
% IMPORTANT: prj_anchorLOO_onCorrectAxes below is built with the ORIGINAL
% (unmasked) projectGlmRezC_toPerMouseLOOAnchors function, NOT the
% correct-trials-only variant that NaN-masks incorrect trials per axis.
%
% Why: the plots in this script deliberately compare DIFFERENT trial
% types on the SAME axis within a single session (e.g. CR vs Hit on a
% NoGo-related axis, to visualize selectivity). If we used the NaN-masked
% projection, Hit trials would show up as NaN on a NoGo-related axis
% (since that axis's "valid" trial type is CR-only under
% CorrectOnlyByGroup), which would silently break exactly this kind of
% within-session comparison plot.
%
% The methodological fix (correct-trials-only) lives in the AXES
% themselves (anchorLOO_correctTrials) -- that is what actually changes
% what each axis means. For within-session trial-type comparisons and
% for across-session single-trial-type trajectories (which already do
% their own explicit trial-type filtering via trIdC, e.g.
% 'trialField', "nogoI"), we want the full, unmasked trial-level
% projection onto those improved axes.
%
% Set hdrC_all to the same expert-header matrix used to build
% anchorLOO_correctTrials. If it is already in your workspace from the
% batch pipeline, this line is a no-op; otherwise load/reconstruct it
% before running this block.
% ---------------------------------------------------------------------
if ~exist('hdrC_all', 'var') || isempty(hdrC_all)
    error(['hdrC_all is not defined. Set it to the expert-header matrix ' ...
        'used to build anchorLOO_correctTrials before running this script.']);
end

prj_anchorLOO_onCorrectAxes = projectGlmRezC_toPerMouseLOOAnchors( ...
    anchorLOO_correctTrials, headerC, glmRezC, trIdC, hdrC_all, ...
    'WhichAxes', "A", ...
    'ProjectWhichY', "Yz", ...
    'UseLOOForExpertSessions', true, ...
    'StoreAsSingle', true, ...
    'Verbose', true, ...
    'VerboseEvery', 25);

%% -------------------------------------------------------------------------
% Pick a session by header and extract its projected trajectories
% -------------------------------------------------------------------------
header = "m1045_122424";

% Find [row, col] for this header in headerC
headerI = findHeaderInHeaderC(headerC, header);
assert(~isempty(headerI), 'Header not found in headerC: %s', header);

% Choose projection type:
%   - global  : uses pooled expert anchor across mice
%   - perMouse: uses per-mouse expert anchor
% useProjType = "perMouse";  % "global" | "perMouse"
% 
% switch useProjType
%     case "global"
%         zTrj = prj_glmA.global.ZC{headerI(1), headerI(2)};     % [N x T x Dglobal]
%         axisNameList = prj_glmA.global.names;                  % 1 x Dglobal
%     case "perMouse"
%         zTrj = prj_glmA.perMouse.ZC{headerI(1), headerI(2)};   % [N x T x Dmouse]
%         axisNameList = prj_glmA.perMouse.namesC{headerI(1), headerI(2)}; % 1 x Dmouse
%     otherwise
%         error('Unknown useProjType: %s', useProjType);
% end
axisNameList = rezProjStats_correctTrials.perMouse{1}.axes.names; 

% Trial-type indices for this session (struct with fields: hitI, crI, faI, missI, ...)
trI = trIdC{headerI(1), headerI(2)};
assert(isstruct(trI), 'trIdC{%d,%d} must be a struct.', headerI(1), headerI(2));

% Time vector for bins (seconds); pulled from glmRezC
% NOTE: This requires glmRezC to be loaded above.
%timestamps = -0.85:0.05:4.95; %glmRezC{headerI(1), headerI(2)}.decBins.time;

%% -------------------------------------------------------------------------
% Pull projected trajectories for selected session from new LOO projection object
% -------------------------------------------------------------------------

header = "m1045_122424";

% Find [row, col] for this header
headerI = findHeaderInHeaderC(headerC, header);
assert(~isempty(headerI), 'Header not found in headerC: %s', header);

% New replacement for old zTrj
zTrj = double(prj_anchorLOO_onCorrectAxes.perMouse.ZC{headerI(1), headerI(2)});
assert(~isempty(zTrj), 'No projected trajectory found for %s.', header);

% Axis names for this session
axisNameList = prj_anchorLOO_onCorrectAxes.perMouse.namesC{headerI(1), headerI(2)};
assert(~isempty(axisNameList), 'No axis names found for %s.', header);

% Time vector
timestamps = prj_anchorLOO_onCorrectAxes.perMouse.timeC{headerI(1), headerI(2)};
timestamps = timestamps(:)';

% Trial identity struct
trI = trIdC{headerI(1), headerI(2)};
assert(isstruct(trI), 'trIdC{%d,%d} must be a struct.', headerI(1), headerI(2));

% Map target axis names to indices
glmTargetCols_nogoToneOff = glmNameToColumns(axisNameList, ...
    {'ToneOffNoGo_1','ToneOffNoGo_2','ToneOffNoGo_3'});

assert(~isempty(glmTargetCols_nogoToneOff), ...
    'Requested ToneOffNoGo axes not found in axisNameList.');

% Optional sanity check: did this session use LOO?
fprintf('%s | usedLOO = %d | source = %s\n', ...
    header, ...
    prj_anchorLOO_onCorrectAxes.perMouse.usedLOO(headerI(1), headerI(2)), ...
    prj_anchorLOO_onCorrectAxes.perMouse.sourceC{headerI(1), headerI(2)});

% Day 4 session mark
day4MarkC = {
    "m1044", "121124"; 
    "m1045", "121124"; 
    "m1048", "121724";
    "m1049", "121324"; 
    "m1092", "100324";
    "m1094", "100224";
    "m1237", "100224";
    "m1613", "042825"; 
    "m1859", "041725"; 
    "m1873", "042125"
    }; 

%% -------------------------------------------------------------------------
% Map GLM axis names -> column indices (by name)
% -------------------------------------------------------------------------
% These are examples; keep whichever you need.
glmTargetCols_goToneOn = glmNameToColumns(axisNameList, ...
    {'GoToneOn_1','GoToneOn_2','GoToneOn_3'});

glmTargetCols_goToneOff = glmNameToColumns(axisNameList, ...
    {'ToneOffGo_1','ToneOffGo_2','ToneOffGo_3'});

glmTargetCols_nogoToneOn = glmNameToColumns(axisNameList, ...
    {'NoGoToneOn_1','NoGoToneOn_2','NoGoToneOn_3'});

glmTargetCols_nogoToneOff = glmNameToColumns(axisNameList, ...
    {'ToneOffNoGo_1','ToneOffNoGo_2','ToneOffNoGo_3'});

% Quick sanity
assert(~isempty(glmTargetCols_nogoToneOff), 'Requested axis names not found in axisNameList.');

%% -------------------------------------------------------------------------
% Plot: single-session trial-averaged trajectory along ONE axis
%   Example: ToneOffNoGo_1, compare CR vs Hit
% -------------------------------------------------------------------------
plotAxis = glmTargetCols_nogoToneOn(1);

plotTrAvgGlmTrjs(zTrj, {trI.crI trI.hitI}, plotAxis, timestamps, ...
    'tBounds', [0 5], ...
    'smoothingFactor', 5, ...
    'FadeToWhite', 1, ...
    'FadeN', 120, ...
    'LegendC', {'CR','Hit'});

%% -------------------------------------------------------------------------
% Plot: across-session evolution of ONE axis (1-D) for a given trial type
% -------------------------------------------------------------------------
% NOTE: perMouseAcrossSessionPrjScoreTrajectories expects "timestamps" as an input;
% if different sessions have different time vectors, you may want to pass the
% appropriate one or modify the function to use per-session timestamps.
perMouseAcrossSessionPrjScoreTrajectories(prj_anchorLOO_onCorrectAxes, headerC, timestamps, ...
    'mouseId', "m1045", ...
    'trIdC', trIdC, ...
    'trialField', "nogoI", ...
    'projType', "perMouse", ...
    'targetName', "NoGoToneOn_1", ...
    'day4MarkC', day4MarkC, ...
    'dateLaterThan', [], ...
    'dateEarlierThan', [], ...
    'tBounds', [], ...
    'lineColor', [0 0 1], ...
    'smoothingFactor', 5, ...
    'figSaveDir', compatiblepath('Z:\Rodent Data\dualImaging_parkj\collectFigure\motifTDR\glmTDR_correctTrials'));

perMouseAcrossSessionCollapsedTrajectories( ...
    rezProjStats_correctTrials, timestamps, ...
    'mouseId', 'm1092', ...
    'targetName', 'ToneOffNoGo_3', ...
    'trialField', 'nogo', ...
    'FlipDiffForNoGoAxes', true, ...
    'day4MarkC', day4MarkC, ...
    'tBounds', [-1 5], ...
    'smoothingFactor', 3, ...
    'MakeFigure', true, ...
    'figSaveDir', compatiblepath('Z:\Rodent Data\dualImaging_parkj\collectFigure\motifTDR\glmTDR_correctTrials'));

perMouseAcrossSessionCollapsedTrajectories( ...
    rezProjStats_correctTrials, timestamps, ...
    'mouseId', 'm1094', ...
    'targetName', 'ToneOffGo_3', ...
    'trialField', 'go', ...
    'day4MarkC', day4MarkC, ...
    'lineColor', [0 0 1], ...
    'tBounds', [-1 5], ...
    'smoothingFactor', 3, ...
    'MakeFigure', true);

%%
imageProjectedGoNogoTrials(prj_anchorLOO_onCorrectAxes, headerC, trIdC, "m1045_122424", "perSession", ...
    "axisName", 'GoToneOn_1', 'CLim', [-2 8])



%% %%%%%%%%%%%%%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function idx = findHeaderInHeaderC(headerC, header)
%FINDHEADERINHEADERC  Find the [row, col] location of a session header in headerC.
%
% idx = findHeaderInHeaderC(headerC, header)
%
% INPUTS
%   headerC : cell array (J x S) of session headers (strings/chars or empty)
%   header  : char or string, e.g. 'm1045_122424'
%
% OUTPUT
%   idx     : [row, col] of the FIRST match (row-major order)
%             [] if no match is found
%
% NOTES
%   - Empty cells in headerC are ignored
%   - Exact string match is used
%   - If the header appears multiple times, the first occurrence is returned

% -------------------- sanity --------------------
if nargin < 2 || isempty(headerC) || isempty(header)
    idx = [];
    return;
end

header = string(header);

% -------------------- flatten + filter empties --------------------
hdrFlat = headerC(:);
isNonEmpty = ~cellfun(@isempty, hdrFlat);
if ~any(isNonEmpty)
    idx = [];
    return;
end

hdrFlat  = hdrFlat(isNonEmpty);
hdrFlatS = string(hdrFlat);

% -------------------- exact match --------------------
hit = find(hdrFlatS == header, 1, 'first');
if isempty(hit)
    idx = [];
    return;
end

% -------------------- map back to [row, col] --------------------
linIdxAll = find(isNonEmpty);     % linear indices into headerC
linIdx    = linIdxAll(hit);

[row, col] = ind2sub(size(headerC), linIdx);
idx = [row, col];
end

function glmTargetCols = glmNameToColumns(glmNameListC, targetNameC)
%GLMNAMETOCOLUMNS  Map target GLM axis names to column indices.
%
% glmTargetCols = glmNameToColumns(glmNameListC, targetNameC)
%
% INPUTS
%   glmNameListC : 1xD (or Dx1) cell array / string array of GLM axis names
%                  e.g. {'GoToneOn_1','GoToneOn_2',...}
%   targetNameC  : cell array / string array of names to find
%                  e.g. {'GoToneOn_1','GoToneOn_2','GoToneOn_3'}
%
% OUTPUT
%   glmTargetCols : numeric row vector of column indices into glmNameListC
%                   (exact-match only; order follows targetNameC)
%
% NOTES
%   - Exact string match only (case-sensitive)
%   - Missing target names are silently ignored (returns fewer indices)
%   - Robust to char/string mixtures

% -------------------- sanity --------------------
if isempty(glmNameListC) || isempty(targetNameC)
    glmTargetCols = [];
    return;
end

% Normalize to string row vectors
glmNames = string(glmNameListC(:))';    % 1 x D
targets  = string(targetNameC(:))';     % 1 x T

glmTargetCols = [];

% -------------------- exact matching --------------------
for i = 1:numel(targets)
    idx = find(glmNames == targets(i), 1, 'first');
    if ~isempty(idx)
        glmTargetCols(end+1) = idx; %#ok<AGROW>
    end
end
end