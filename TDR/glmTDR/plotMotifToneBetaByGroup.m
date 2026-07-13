function h = plotMotifToneBetaByGroup(synthBetaByGroup, synthNames, groupNames, motifId, nAnimalsByGroup, varargin)
%PLOTMOTIFTONEBETABYGROUP  Plot group-mean tone-predictor beta profiles, one panel per group.
%
% SYNOPSIS
%   h = plotMotifToneBetaByGroup(synthBetaByGroup, synthNames, groupNames, ...
%           motifId, nAnimalsByGroup, ...)
%
% DESCRIPTION
%   Plots one panel per group (e.g. fast learners, slow learners) showing
%   that group's mean toneOnGo/toneOnNoGo/toneOffGo/toneOffNoGo beta
%   profile for a single motif -- same tone/response block layout, Go
%   (blue) / NoGo (red) coloring, and onset (filled)/offset (open) style
%   as plot_beta_timeBinPaired.m, since this function is just a thin
%   wrapper that calls it once per group on the group-mean synthetic data.
%
%   Because the point of this figure is to directly compare magnitude and
%   direction of tuning between groups, panels share a common y-axis by
%   default ('matchYLim', true): each group is first plotted in an
%   invisible "draft" pass to measure its natural auto-scaled range, then
%   re-plotted for real using the resulting shared range (passed in via
%   plot_beta_timeBinPaired's 'yLim' option), so the divider line and
%   "Tone Onset"/"Tone Offset" labels are correctly positioned for the
%   FINAL shared range rather than each group's own range.
%
% INPUTS
%   synthBetaByGroup, synthNames, groupNames, nAnimalsByGroup : outputs of
%       computeMotifToneBetaByGroup.m.
%   motifId : the real motif number this data came from (used only for
%       titling -- the underlying synthetic beta has a single column).
%
% NAME-VALUE ARGS
%   'plotStyle'     : 'bar' or 'curve', passed through to
%                      plot_beta_timeBinPaired.m. Default: 'bar'.
%   'includeOffset' : passed through. Default: true.
%   'matchYLim'     : logical, harmonize y-axis across all group panels.
%                      Default: true.
%   'figureScaleFactor', 'figureWidthFactor' : passed through. Defaults
%                      match plot_beta_timeBinPaired.m's own defaults.
%   'visible'       : 'on'/'off' for the FINAL panels (the measurement
%                      draft pass is always invisible regardless of this
%                      setting). Default: 'on'.
%   'figSaveDir', 'figSaveKeyword' : save-to-PDF options; each group's
%                      panel is saved with that group's name as 'header'.
%
% OUTPUT
%   h : struct with h.byGroup.(groupName) = the plot_beta_timeBinPaired
%       output struct for that group (fig/ax/etc.), plus h.sharedYLim (if
%       matchYLim was used) and h.opt.
%
% EXAMPLE
%   groupDefs.fast = fast_learners;
%   groupDefs.slow = slow_learners;
%   [synthBetaByGroup, synthNames, ~, groupNames, nAnimalsByGroup] = ...
%       computeMotifToneBetaByGroup(glmRezC, glmLabelC, groupDefs, 20);
%   h = plotMotifToneBetaByGroup(synthBetaByGroup, synthNames, groupNames, 20, nAnimalsByGroup, ...
%           'figSaveDir', figSaveDir, 'figSaveKeyword', 'toneBetaByGroup');
%
% See also: computeMotifToneBetaByGroup, plot_beta_timeBinPaired

p = inputParser;
p.addParameter('plotStyle', 'bar', @(s) any(strcmpi(s, {'bar','curve'})));
p.addParameter('includeOffset', true, @(x) islogical(x) && isscalar(x));
p.addParameter('matchYLim', true, @(x) islogical(x) && isscalar(x));
p.addParameter('figureScaleFactor', 2, @(x) isnumeric(x) && x>0);
p.addParameter('figureWidthFactor', 2.1, @(x) isnumeric(x) && x>0);
p.addParameter('visible', 'on', @(s) any(strcmpi(s, {'on','off'})));
p.addParameter('figSaveDir', {}, @(x) isempty(x) || ischar(x) || isstring(x) || iscell(x));
p.addParameter('figSaveKeyword', '', @(s) ischar(s) || isstring(s));
p.parse(varargin{:});
opt = p.Results;

figSaveDir = opt.figSaveDir;
if iscell(figSaveDir)
    if isempty(figSaveDir), figSaveDir = ""; else, figSaveDir = string(figSaveDir{1}); end
else
    figSaveDir = string(figSaveDir);
end
figSaveKeyword = string(opt.figSaveKeyword);

validGroups = {};
for gi = 1:numel(groupNames)
    gName = groupNames{gi};
    if isfield(synthBetaByGroup, gName) && ~isempty(synthBetaByGroup.(gName))
        validGroups{end+1} = gName; %#ok<AGROW>
    else
        warning('plotMotifToneBetaByGroup:skipEmptyGroup', 'Skipping group "%s" (no data).', gName);
    end
end
assert(~isempty(validGroups), 'No group has valid data to plot.');

h = struct();
h.byGroup = struct();
h.opt = opt;

nAfun = @(gName) getfield_or_zero_(nAnimalsByGroup, gName);

% -------- pass 1: draft (invisible) to measure each group's natural y-range --------
sharedYLim = [];
if opt.matchYLim
    ylAll = nan(numel(validGroups), 2);
    for gi = 1:numel(validGroups)
        gName = validGroups{gi};
        hDraft = plot_beta_timeBinPaired(synthBetaByGroup.(gName), synthNames, 1, ...
            'plotStyle', opt.plotStyle, ...
            'includeOffset', opt.includeOffset, ...
            'visible', 'off');
        ylAll(gi,:) = ylim(hDraft.ax);
        close(hDraft.fig);
    end
    sharedYLim = [min(ylAll(:,1)), max(ylAll(:,2))];
    h.sharedYLim = sharedYLim;
end

% -------- pass 2: final plots, one per group --------
for gi = 1:numel(validGroups)
    gName = validGroups{gi};
    nA = nAfun(gName);

    plotArgs = { ...
        'plotStyle', opt.plotStyle, ...
        'includeOffset', opt.includeOffset, ...
        'figureScaleFactor', opt.figureScaleFactor, ...
        'figureWidthFactor', opt.figureWidthFactor, ...
        'visible', opt.visible, ...
        'header', gName, ...
        'title', sprintf('%s learners (n=%d): motif %d tone \\beta', gName, nA, motifId), ...
        'figSaveDir', figSaveDir, ...
        'figSaveKeyword', strjoin(strings_nonempty_({figSaveKeyword, "byGroup"}), "_"), ...
        'saveMotifId', motifId ...
        };
    if ~isempty(sharedYLim)
        plotArgs = [plotArgs, {'yLim', sharedYLim}]; %#ok<AGROW>
    end

    h.byGroup.(gName) = plot_beta_timeBinPaired(synthBetaByGroup.(gName), synthNames, 1, plotArgs{:});
end

end % function


% ===== helper: safe struct field lookup, default 0 =====
function v = getfield_or_zero_(s, fld)
if isfield(s, fld)
    v = s.(fld);
else
    v = 0;
end
end


% ===== helper: drop empty strings before joining filename parts =====
function out = strings_nonempty_(parts)
parts = string(parts);
out = parts(strlength(parts) > 0);
end
