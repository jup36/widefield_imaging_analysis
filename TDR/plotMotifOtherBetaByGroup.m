function h = plotMotifOtherBetaByGroup(groupOtherBetaTable, groupNames, motifId, varargin)
%PLOTMOTIFOTHERBETABYGROUP  Plot group-mean "other" (non-tone) predictor beta bars, one panel per group.
%
% SYNOPSIS
%   h = plotMotifOtherBetaByGroup(groupOtherBetaTable, groupNames, motifId, ...)
%
% DESCRIPTION
%   The non-tone companion to plotMotifToneBetaByGroup.m. Plots one panel
%   per group (e.g. fast learners, slow learners) showing that group's
%   mean beta for every non-tone predictor base (lick, water, airpuff,
%   pupil, whisker, ...), collapsed exactly the way
%   plot_beta_groupMeanBars.m does with 'SplitOtherByBase', true,
%   'plotTonePredictors', false -- since this function is a thin wrapper
%   that calls it once per group on synthetic single-value-per-base data.
%
%   Both groups' panels use the SAME set of base names, in the SAME
%   x-axis order (the union across groups, sorted alphabetically unless
%   'baseOrder' is given), so bars are positionally comparable across
%   panels. A base missing for a given group (e.g. a continuous predictor
%   never tracked in any of that group's sessions) is left as a gap
%   (NaN) rather than a misleading zero bar.
%
%   Panels share a common y-axis by default ('matchYLim', true), measured
%   via an invisible draft pass per group (same approach as
%   plotMotifToneBetaByGroup.m) so the shared range is baked in BEFORE
%   each group's real panel is plotted/saved.
%
% INPUTS
%   groupOtherBetaTable, groupNames : outputs of computeMotifOtherBetaByGroup.m.
%   motifId : the real motif number this data came from (used only for
%       titling and the saved filename -- the underlying synthetic beta
%       fed to plot_beta_groupMeanBars has a single column).
%
% NAME-VALUE ARGS
%   'baseOrder'   : cellstr giving an explicit x-axis order for the base
%                   names. Default: [] (alphabetical union of whatever
%                   base names are present in groupOtherBetaTable).
%   'matchYLim'   : logical, harmonize y-axis across all group panels.
%                   Default: true.
%   'visible'     : 'on'/'off' for the FINAL panels (the measurement
%                   draft pass is always invisible regardless of this
%                   setting). Default: 'on'.
%   'figSaveDir', 'figSaveKeyword' : save-to-PDF options; each group's
%                   panel is saved with that group's name as 'header'.
%
% OUTPUT
%   h : struct with h.byGroup.(groupName) = the plot_beta_groupMeanBars
%       output struct for that group (fig/ax/etc.), plus h.sharedYLim (if
%       matchYLim was used), h.baseNames (the x-axis order used), and h.opt.
%
% EXAMPLE
%   groupDefs.fast = fast_learners;
%   groupDefs.slow = slow_learners;
%   [groupOtherBetaTable, ~, groupNames] = ...
%       computeMotifOtherBetaByGroup(glmRezC, glmLabelC, groupDefs, 20);
%   h = plotMotifOtherBetaByGroup(groupOtherBetaTable, groupNames, 20, ...
%           'figSaveDir', figSaveDir, 'figSaveKeyword', 'otherBetaByGroup');
%
% See also: computeMotifOtherBetaByGroup, plot_beta_groupMeanBars, plotMotifToneBetaByGroup

p = inputParser;
p.addParameter('baseOrder', {}, @(x) iscell(x) || isstring(x));
p.addParameter('matchYLim', true, @(x) islogical(x) && isscalar(x));
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

validGroups = intersect(groupNames, unique(groupOtherBetaTable.group), 'stable');
assert(~isempty(validGroups), 'No group in groupNames has data in groupOtherBetaTable.');

if isempty(opt.baseOrder)
    baseNames = sort(unique(groupOtherBetaTable.baseName));
else
    baseNames = cellstr(string(opt.baseOrder));
end
nBases = numel(baseNames);

% -------- build synthetic beta/names per group (NaN where a base is missing) --------
synthBetaByGroup = struct();
for gi = 1:numel(validGroups)
    gName = validGroups{gi};
    vec = nan(nBases,1);
    for bI = 1:nBases
        rowsI = strcmp(groupOtherBetaTable.group, gName) & strcmp(groupOtherBetaTable.baseName, baseNames{bI});
        if any(rowsI)
            vec(bI) = groupOtherBetaTable.meanBeta(rowsI);
        end
    end
    if any(isnan(vec))
        missingBases = baseNames(isnan(vec));
        warning('plotMotifOtherBetaByGroup:missingBase', ...
            'Group "%s" has no data for base(s): %s -- left as a gap.', gName, strjoin(missingBases, ', '));
    end
    synthBetaByGroup.(gName) = vec;
end
synthNames = baseNames(:)';   % literal base names as predictor names -- each classifies as "other" with itself as baseKey

h = struct();
h.opt = opt;
h.baseNames = baseNames;

% -------- pass 1: draft (invisible) to measure each group's natural y-range --------
sharedYLim = [];
if opt.matchYLim
    ylAll = nan(numel(validGroups), 2);
    for gi = 1:numel(validGroups)
        gName = validGroups{gi};
        hDraft = plot_beta_groupMeanBars(synthBetaByGroup.(gName), synthNames, 1, ...
            'SplitOtherByBase', true, 'plotTonePredictors', false, 'visible', 'off');
        ylAll(gi,:) = ylim(hDraft.ax);
        close(hDraft.fig);
    end
    sharedYLim = [min(ylAll(:,1)), max(ylAll(:,2))];
    h.sharedYLim = sharedYLim;
end

% -------- pass 2: final plots, one per group --------
h.byGroup = struct();
for gi = 1:numel(validGroups)
    gName = validGroups{gi};

    plotArgs = { ...
        'SplitOtherByBase', true, ...
        'plotTonePredictors', false, ...
        'visible', opt.visible, ...
        'header', gName, ...
        'title', sprintf('%s learners: motif %d mean \\beta (other predictors)', gName, motifId), ...
        'figSaveDir', figSaveDir, ...
        'figSaveKeyword', strjoin(strings_nonempty_({figSaveKeyword, "otherByGroup"}), "_"), ...
        'saveMotifId', motifId ...
        };
    if ~isempty(sharedYLim)
        plotArgs = [plotArgs, {'yLim', sharedYLim}]; %#ok<AGROW>
    end

    h.byGroup.(gName) = plot_beta_groupMeanBars(synthBetaByGroup.(gName), synthNames, 1, plotArgs{:});
end

end % function


% ===== helper: drop empty strings before joining filename parts =====
function out = strings_nonempty_(parts)
parts = string(parts);
out = parts(strlength(parts) > 0);
end
