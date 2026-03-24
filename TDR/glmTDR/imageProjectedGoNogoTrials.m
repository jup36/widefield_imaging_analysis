function [fig, ax, info] = imageProjectedGoNogoTrials(prj_glmA, headerC, trIdC, header, whichProj, varargin)
%IMAGEPROJECTEDGONOGOTRIALS  Diagnostic imagesc of projected trajectories for Go vs NoGo.
%
% [fig, ax, info] = imageProjectedGoNogoTrials(prj_glmA, headerC, trIdC, header, whichProj, 'Name', value, ...)
%
% REQUIRED INPUTS
%   prj_glmA   : output of projectGlmRezC_toAnchors (must contain .global/.perMouse/.perSession)
%   headerC    : (J x S) cell array of session headers (used to locate i,j)
%   trIdC      : (J x S) cell array of trial-ID structs (e.g. from trialTypeInfo... per session)
%                Must contain fields at least: goI, nogoI. If correctTrialsOnly=true,
%                must contain hitI and crI.
%   header     : char/string, e.g. "m1045_122424"
%   whichProj  : "global" | "perMouse" | "perSession"
%
% NAME-VALUE
%   'axisName'          : REQUIRED. e.g. "GoToneOn_1" (case-insensitive match)
%   'correctTrialsOnly' : false (default). If true: compare Hit vs CR (hitI vs crI)
%   'CLim'              : [] (default). If provided, sets caxis([lo hi]) for BOTH panels
%   'ShowColorbar'      : true (default)
%   'FigureNamePrefix'  : "imageProjectedGoNogoTrials" (default)
%
% OUTPUTS
%   fig  : figure handle
%   ax   : 1x2 axes handles (left=Go/Hit, right=NoGo/CR)
%   info : struct with diagnostic info (indices, axisIdx, names used, etc.)
%
% Notes
% - Assumes ZC{j,s} is [N x nW x nAxis].
% - Axis name matching uses namesC{j,s} when available (perMouse/perSession),
%   otherwise uses prj_glmA.global.names for global.
%
% JC lab-style: explicit + robust.

% -------------------- parse --------------------
p = inputParser;
p.addRequired('prj_glmA', @(x)isstruct(x) && ~isempty(x));
p.addRequired('headerC', @(x)iscell(x));
p.addRequired('trIdC',   @(x)iscell(x) && isequal(size(x), size(headerC)));
p.addRequired('header',  @(x)ischar(x) || isstring(x));
p.addRequired('whichProj', @(x)ischar(x) || isstring(x));

p.addParameter('axisName', [], @(x)ischar(x)||isstring(x));
p.addParameter('correctTrialsOnly', false, @(x)islogical(x)&&isscalar(x));
p.addParameter('CLim', [], @(x) isempty(x) || (isnumeric(x)&&numel(x)==2&&all(isfinite(x))));
p.addParameter('ShowColorbar', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('FigureNamePrefix', "imageProjectedGoNogoTrials", @(x)ischar(x)||isstring(x));

p.parse(prj_glmA, headerC, trIdC, header, whichProj, varargin{:});
opt = p.Results;

if isempty(opt.axisName)
    error('imageProjectedGoNogoTrials:MissingAxisName', ...
        'You must provide Name-Value ''axisName'', e.g., ''GoToneOn_1''.');
end

hdr = char(string(header));
whichProj = lower(string(whichProj));
axisNameQ = string(opt.axisName);

% -------------------- locate session (j,s) from header --------------------
[j, s, ok] = find_header_in_matrix_(headerC, hdr);
if ~ok
    error('imageProjectedGoNogoTrials:HeaderNotFound', ...
        'Header "%s" not found in headerC.', hdr);
end

% -------------------- fetch Z + axis names for this projection --------------------
[Z, names, okZ] = fetch_projection_Z_(prj_glmA, whichProj, j, s);
if ~okZ
    error('imageProjectedGoNogoTrials:ProjectionMissing', ...
        'Projection "%s" missing at (%d,%d) for header "%s".', whichProj, j, s, hdr);
end

% Z is [N x nW x nAxis]
if ndims(Z) ~= 3
    error('imageProjectedGoNogoTrials:BadZ', ...
        'Expected Z to be [N x nW x nAxis], got ndims=%d.', ndims(Z));
end

% -------------------- axis index by name --------------------
axisIdx = match_axis_name_(names, axisNameQ);
if isempty(axisIdx)
    error('imageProjectedGoNogoTrials:AxisNameNotFound', ...
        'axisName "%s" not found for projection "%s" in session "%s".', ...
        char(axisNameQ), char(whichProj), hdr);
end

% -------------------- trial indices --------------------
tid = trIdC{j,s};
if isempty(tid) || ~isstruct(tid)
    error('imageProjectedGoNogoTrials:MissingTrId', ...
        'trIdC{%d,%d} is empty or not a struct for header "%s".', j, s, hdr);
end

if opt.correctTrialsOnly
    req = {'hitI','crI'};
    assert(all(isfield(tid, req)), 'correctTrialsOnly=true requires fields: hitI and crI.');
    leftI  = tid.hitI(:);
    rightI = tid.crI(:);
    leftLabel  = "Hit";
    rightLabel = "CR";
else
    req = {'goI','nogoI'};
    assert(all(isfield(tid, req)), 'correctTrialsOnly=false requires fields: goI and nogoI.');
    leftI  = tid.goI(:);
    rightI = tid.nogoI(:);
    leftLabel  = "Go";
    rightLabel = "NoGo";
end

% -------------------- extract axis slice --------------------
Zax = Z(:,:,axisIdx);      % [N x nW]
Zleft  = Zax(leftI,  :);
Zright = Zax(rightI, :);

% -------------------- plotting --------------------
figName = sprintf('%s | %s | %s | axis=%s', ...
    char(string(opt.FigureNamePrefix)), hdr, char(whichProj), char(axisNameQ));

fig = figure('Name', figName, 'Color', 'w');
tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

ax = gobjects(1,2);

ax(1) = nexttile;
imagesc(Zleft);
axis tight;
xlabel('time bin');
ylabel('trials');
title(sprintf('%s (%d trials)', leftLabel, size(Zleft,1)), 'Interpreter','none');

ax(2) = nexttile;
imagesc(Zright);
axis tight;
xlabel('time bin');
ylabel('trials');
title(sprintf('%s (%d trials)', rightLabel, size(Zright,1)), 'Interpreter','none');

% shared CLim if requested
if ~isempty(opt.CLim)
    caxis(ax(1), opt.CLim);
    caxis(ax(2), opt.CLim);
end

% shared colormap (default)
colormap(fig, parula);

if opt.ShowColorbar
    cb = colorbar(ax(2));
    cb.Location = 'eastoutside';
end

sgtitle(sprintf('%s | axis: %s (idx=%d)', hdr, char(axisNameQ), axisIdx), 'Interpreter','none');

% -------------------- outputs --------------------
info = struct();
info.header = hdr;
info.j = j;
info.s = s;
info.whichProj = char(whichProj);
info.axisName = char(axisNameQ);
info.axisIdx = axisIdx;
info.names = names;
info.correctTrialsOnly = opt.correctTrialsOnly;
info.leftLabel = char(leftLabel);
info.rightLabel = char(rightLabel);
info.nLeft = size(Zleft,1);
info.nRight = size(Zright,1);

end

%% ============================== HELPERS ==============================

function [j, s, ok] = find_header_in_matrix_(headerC, hdr)
ok = false; j = []; s = [];
hdrS = string(hdr);

J = size(headerC,1);
S = size(headerC,2);

for jj = 1:J
    for ss = 1:S
        h = headerC{jj,ss};
        if isempty(h), continue; end
        if string(h) == hdrS
            j = jj; s = ss; ok = true;
            return;
        end
    end
end
end

function [Z, names, ok] = fetch_projection_Z_(prj, whichProj, j, s)
Z = [];
names = {};
ok = false;

switch lower(string(whichProj))
    case "global"
        if ~isfield(prj,'global') || ~isfield(prj.global,'ZC') || isempty(prj.global.ZC{j,s})
            return;
        end
        Z = prj.global.ZC{j,s};
        if isfield(prj.global,'names') && ~isempty(prj.global.names)
            names = cellstr(string(prj.global.names));
        end
        ok = ~isempty(Z);

    case "permouse"
        if ~isfield(prj,'perMouse') || ~isfield(prj.perMouse,'ZC') || isempty(prj.perMouse.ZC{j,s})
            return;
        end
        Z = prj.perMouse.ZC{j,s};
        if isfield(prj.perMouse,'namesC') && ~isempty(prj.perMouse.namesC{j,s})
            names = cellstr(string(prj.perMouse.namesC{j,s}));
        end
        ok = ~isempty(Z);

    case "persession"
        if ~isfield(prj,'perSession') || ~isfield(prj.perSession,'ZC') || isempty(prj.perSession.ZC{j,s})
            return;
        end
        Z = prj.perSession.ZC{j,s};
        if isfield(prj.perSession,'namesC') && ~isempty(prj.perSession.namesC{j,s})
            names = cellstr(string(prj.perSession.namesC{j,s}));
        end
        ok = ~isempty(Z);

    otherwise
        error('Unknown whichProj="%s". Use "global", "perMouse", or "perSession".', char(whichProj));
end
end

function axisIdx = match_axis_name_(names, axisNameQ)
axisIdx = [];

if isempty(names)
    return;
end

namesS = string(names);
q = lower(string(axisNameQ));

% exact case-insensitive match first
idx = find(lower(namesS) == q, 1, 'first');
if ~isempty(idx)
    axisIdx = idx;
    return;
end

% fallback: contains match (useful if user types without underscore variants)
idx = find(contains(lower(namesS), q), 1, 'first');
if ~isempty(idx)
    axisIdx = idx;
end
end