function h = plot_cvR2_acrossSessions(glmCvR2C_collect, animalIDs, varargin)
%PLOT_CVR2_ACROSSSESSIONS  Line plot of cross-validated global R^2 across sessions, per animal.
%
% SYNOPSIS
%   h = plot_cvR2_acrossSessions(glmCvR2C_collect, animalIDs, ...)
%
% DESCRIPTION
%   Plots one line per animal, tracing that animal's cvR2_global values
%   across its valid sessions (x = sequential valid-session index, not
%   calendar session number -- see note in collect_cvR2_perAnimal.m).
%   Lines are color-coded by animal ID using a distinguishable colormap
%   (MATLAB's built-in 'lines').
%
% INPUTS
%   glmCvR2C_collect : [nAnimals x 1] cell array of per-animal R^2 vectors,
%                      as produced by collect_cvR2_perAnimal.m.
%   animalIDs        : [nAnimals x 1] cellstr of animal ID labels, same
%                      length as glmCvR2C_collect.
%
% NAME-VALUE ARGS
%   'title'      : axes title. Default: auto.
%   'xlabelStr'  : x-axis label. Default: 'Session (sequential)'.
%   'ylabelStr'  : y-axis label. Default: 'Cross-validated global R^2'.
%   'figureScaleFactor', 'figureWidthFactor', 'visible' : figure sizing /
%                  visibility, consistent with other plotting functions
%                  in this pipeline.
%   'figSaveDir', 'header', 'figSaveKeyword' : save-to-PDF options, same
%                  convention as plot_beta_with_labels.m etc.
%   'colorMap'   : [nAnimals x 3] RGB matrix aligned to animalIDs, to
%                  override the default 'lines' colormap (e.g. from
%                  get_colors_for_animalIDs.m). Default: [] (use 'lines').
%   'groupLabels': cellstr, same length as animalIDs (e.g. 'fast'/'slow'
%                  per animal), used ONLY to reorder the legend -- line
%                  colors/positions are unaffected. Default: {} (legend
%                  follows animalIDs order).
%   'groupOrder' : cellstr giving the desired group order in the legend,
%                  e.g. {'fast','slow'} lists all fast learners first.
%                  Default: {'fast','slow'}.
%
% OUTPUT
%   h : struct with fig/ax handles, per-animal line handles, and opt.
%
% EXAMPLE
%   [glmCvR2C_collect, animalIDs] = collect_cvR2_perAnimal(glmEvRezC, glmLabelC);
%   h = plot_cvR2_acrossSessions(glmCvR2C_collect, animalIDs, ...
%           'figSaveDir', figSaveDir, 'figSaveKeyword', 'cvR2acrossSessions');
%
% See also: collect_cvR2_perAnimal

% -------- parse args --------
p = inputParser;
p.addParameter('title', '', @(s) ischar(s) || isstring(s));
p.addParameter('xlabelStr', 'Session (sequential)', @(s) ischar(s) || isstring(s));
p.addParameter('ylabelStr', 'Cross-validated global R^2', @(s) ischar(s) || isstring(s));
p.addParameter('figureScaleFactor', 1, @(x) isnumeric(x) && x>0);
p.addParameter('figureWidthFactor', 1, @(x) isnumeric(x) && x>0);
p.addParameter('visible', 'on', @(s) any(strcmpi(s, {'on','off'})));
p.addParameter('figSaveDir', {}, @(x) isempty(x) || ischar(x) || isstring(x) || iscell(x));
p.addParameter('header', '', @(s) ischar(s) || isstring(s));
p.addParameter('figSaveKeyword', '', @(s) ischar(s) || isstring(s));
p.addParameter('colorMap', [], @(x) isempty(x) || (isnumeric(x) && size(x,2)==3));
p.addParameter('groupLabels', {}, @(x) isempty(x) || iscell(x));
p.addParameter('groupOrder', {'fast','slow'}, @(x) iscell(x));
p.parse(varargin{:});
opt = p.Results;

% normalize figSaveDir
figSaveDir = opt.figSaveDir;
if iscell(figSaveDir)
    if isempty(figSaveDir), figSaveDir = ""; else, figSaveDir = string(figSaveDir{1}); end
else
    figSaveDir = string(figSaveDir);
end
header         = string(opt.header);
figSaveKeyword = string(opt.figSaveKeyword);

% -------- sanity checks --------
nAnimals = numel(glmCvR2C_collect);
assert(iscell(glmCvR2C_collect), 'glmCvR2C_collect must be a cell array.');
assert(numel(animalIDs)==nAnimals, 'animalIDs must be the same length as glmCvR2C_collect.');
if ~isempty(opt.groupLabels)
    assert(numel(opt.groupLabels)==nAnimals, 'groupLabels must be the same length as animalIDs.');
end

% -------- plot --------
h = struct();
h.fig = figure('Color','w', 'Visible', opt.visible);
set(h.fig, 'Units', 'normalized');
pos = get(h.fig, 'Position');
pos(3:4) = pos(3:4) * opt.figureScaleFactor;
set(h.fig, 'Position', pos);
h.ax = axes('Parent', h.fig); hold(h.ax, 'on');

if ~isempty(opt.colorMap)
    assert(size(opt.colorMap,1)==nAnimals, ...
        'colorMap must have one row per animal (numel(animalIDs)=%d, got %d rows).', ...
        nAnimals, size(opt.colorMap,1));
    colors = opt.colorMap;
else
    colors = lines(max(nAnimals,1));
end

h.lines = gobjects(nAnimals,1);
legHandles = [];
legEntries = {};
includedIdx = [];
maxN = 0;

for i = 1:nAnimals
    vals = glmCvR2C_collect{i};
    if isempty(vals)
        continue;   % nothing to plot for this animal (all sessions missing/invalid)
    end
    x = 1:numel(vals);
    maxN = max(maxN, numel(vals));

    h.lines(i) = plot(h.ax, x, vals, '-o', ...
        'Color', colors(i,:), ...
        'MarkerFaceColor', colors(i,:), ...
        'MarkerEdgeColor', colors(i,:), ...
        'LineWidth', 1.5, 'MarkerSize', 5);

    legHandles(end+1) = h.lines(i); %#ok<AGROW>
    legEntries{end+1}  = animalIDs{i}; %#ok<AGROW>
    includedIdx(end+1) = i; %#ok<AGROW>
end

if maxN > 0
    xticks(h.ax, 1:maxN);
    xlim(h.ax, [0, maxN + 1]);   % exactly one session of margin on each side
end
xlabel(h.ax, opt.xlabelStr);
ylabel(h.ax, opt.ylabelStr);

if strlength(string(opt.title)) > 0
    title(h.ax, opt.title, 'Interpreter', 'none');
else
    title(h.ax, 'Cross-validated global R^2 across sessions', 'Interpreter', 'tex');
end

box(h.ax, 'off');
grid(h.ax, 'on');
set(h.ax, 'TickDir', 'out');

if ~isempty(legHandles)
    if ~isempty(opt.groupLabels)
        legGroupsIncluded = opt.groupLabels(includedIdx);
        groupRank = nan(1, numel(legGroupsIncluded));
        for gi = 1:numel(opt.groupOrder)
            groupRank(strcmpi(legGroupsIncluded, opt.groupOrder{gi})) = gi;
        end
        groupRank(isnan(groupRank)) = numel(opt.groupOrder) + 1;  % unlisted groups go last
        % stable sort by group rank, tie-broken by original position
        [~, sortIdx] = sortrows([groupRank(:), (1:numel(groupRank))']);
        sortIdx = sortIdx(:)';
        legHandles = legHandles(sortIdx);
        legEntries = legEntries(sortIdx);
    end
    legend(h.ax, legHandles, legEntries, 'Location', 'bestoutside');
end

r = pbaspect;
pbaspect([opt.figureWidthFactor*r(1) r(2) r(3)]);

% -------- save figure (optional) --------
if strlength(figSaveDir) > 0
    if ~isfolder(figSaveDir)
        mkdir(figSaveDir);
    end
    dateStr = char(datetime("today","Format","MMddyy"));

    parts = strings(0,1);
    if strlength(header) > 0,         parts(end+1,1) = header; end
    parts(end+1,1) = "cvR2acrossSessions";
    if strlength(figSaveKeyword) > 0, parts(end+1,1) = figSaveKeyword; end
    parts(end+1,1) = dateStr;

    figSaveName = strjoin(parts, "_");
    print(h.fig, fullfile(figSaveDir, figSaveName), '-dpdf', '-painters', '-bestfit');
end

h.opt = opt;

end % function
