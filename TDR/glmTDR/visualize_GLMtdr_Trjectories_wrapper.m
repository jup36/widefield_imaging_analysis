%% =========================================================================
%  GLM-TDR: quick load + single-session trial-avg trajectories + across-session axis timecourses
%
%  SYNOPSIS
%  This script demonstrates two common “readouts” from GLM-based targeted
%  dimensionality reduction (GLM-TDR) results saved in glmTdrRezCollection_020626.mat:
%
%   (1) Single-session visualization:
%       - Extract a session by header (e.g., m1045_122424)
%       - Pull projected trial×time×axis trajectories (global or perMouse)
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
load(fullfile(fileDir, "glmTdrRezCollection_020626.mat"), ...
    "glmRezPathC", "headerC", "trIdC", "anchorLOO", "rezProjStats");

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
axisNameList = rezProjStats.perMouse{1}.axes.names; 

% Trial-type indices for this session (struct with fields: hitI, crI, faI, missI, ...)
trI = trIdC{headerI(1), headerI(2)};
assert(isstruct(trI), 'trIdC{%d,%d} must be a struct.', headerI(1), headerI(2));

% Time vector for bins (seconds); pulled from glmRezC
% NOTE: This requires glmRezC to be loaded above.
timestamps = -0.85:0.05:4.95; %glmRezC{headerI(1), headerI(2)}.decBins.time;

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
plotAxis = glmTargetCols_nogoToneOff(1);

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
perMouseAcrossSessionPrjScoreTrajectories(prj_glmA, headerC, timestamps, ...
    'mouseId', "m1045", ...
    'trIdC', trIdC, ...
    'trialField', "goI", ...        % goI/nogoI definitely preferred to hitI/crI
    'projType', "perMouse", ...       % "global" | "perMouse"
    'targetName', "GoToneOn_2", ...    % GoToneOn_1, NoGoToneOn_1, ToneOffGo_1, ToneOffNoGo_1
    'day4MarkC', day4MarkC, ...
    'dateLaterThan', [], ...
    'dateEarlierThan', [], ...
    'tBounds', [], ...
    'lineColor', [0 0 1], ...
    'smoothingFactor', 5, ...
    'figSaveDir', compatiblepath('Z:\Rodent Data\dualImaging_parkj\collectFigure\motifTDR\glmTDR'));

perMouseAcrossSessionCollapsedTrajectories( ...
    rezProjStats, headerC, timestamps, ...
    'mouseId', 'm1045', ...
    'targetName', 'NoGoToneOn_1', ...
    'trialField', 'nogo', ...
    'tBounds', [-1 5], ...
    'smoothingFactor', 3, ...
    'MakeFigure', true);

perMouseAcrossSessionCollapsedTrajectories( ...
    rezProjStats, timestamps, ...
    'mouseId', 'm1044', ...
    'targetName', 'NoGoToneOn_3', ...
    'trialField', 'nogo', ...
    'day4MarkC', day4MarkC, ...
    'lineColor', [0 0 1], ...
    'tBounds', [-1 5], ...
    'smoothingFactor', 3, ...
    'MakeFigure', true);

%%
imageProjectedGoNogoTrials(prj_glmA, headerC, trIdC, "m1045_122424", "perSession", ...
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
