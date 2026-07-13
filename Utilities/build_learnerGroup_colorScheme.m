function [colorSchemeTable, colorMapLookup] = build_learnerGroup_colorScheme(fast_learners, slow_learners, varargin)
%BUILD_LEARNERGROUP_COLORSCHEME  Build (and optionally save) a reusable animal->color map.
%
% SYNOPSIS
%   [colorSchemeTable, colorMapLookup] = build_learnerGroup_colorScheme(fast_learners, slow_learners, ...)
%
% DESCRIPTION
%   Assigns each animal a color such that:
%     - Fast learners all share one hue family (default: green).
%     - Slow learners all share a different hue family (default: purple).
%     - Within each group, individuals get a shade gradient (darker+more
%       saturated -> lighter+less saturated) at that fixed hue, so
%       individuals stay distinguishable while visibly "belonging" to
%       their group.
%
%   Colors are assigned by an animal's POSITION IN THE MASTER LIST
%   (fast_learners / slow_learners), not by its position in any particular
%   analysis's animalIDs subset. This keeps the color scheme stable across
%   figures even if a given plot is missing some animals (e.g. due to
%   failed sessions) -- pair this with get_colors_for_animalIDs.m to apply
%   the scheme to a specific animalIDs list.
%
%   Hue choice deliberately avoids blue/red, since those are already used
%   elsewhere in this pipeline's figures for Go/NoGo -- reusing them here
%   for a different variable (learning speed) would risk misreading.
%
% INPUTS
%   fast_learners : cellstr of animal IDs, e.g. {'m1044','m1045','m1092','m1094'}
%   slow_learners : cellstr of animal IDs, e.g. {'m1048','m1049','m1613','m1859','m1873'}
%
% NAME-VALUE ARGS
%   'hueFast'    : hue in [0,1] for the fast-learner group. Default: 0.36 (green).
%   'hueSlow'    : hue in [0,1] for the slow-learner group. Default: 0.80 (purple).
%   'satRange'   : [min max] saturation range across individuals within a
%                  group (first-listed animal gets satRange(2), the most
%                  saturated). Default: [0.55 0.85].
%   'valRange'   : [min max] HSV "value" (brightness) range across
%                  individuals (first-listed animal gets valRange(1), the
%                  darkest). Default: [0.45 0.85].
%   'saveDir'    : folder to save the color scheme .mat file to. If empty
%                  (default), nothing is saved.
%   'saveName'   : filename to save as. Default: 'learnerGroupColorScheme.mat'.
%
% OUTPUTS
%   colorSchemeTable : table with columns animalID, group, color (Nx3 RGB),
%                      one row per animal in fast_learners + slow_learners.
%   colorMapLookup   : containers.Map, animalID (char) -> 1x3 RGB, for fast
%                      lookup (used by get_colors_for_animalIDs.m).
%
% EXAMPLE
%   fast_learners = {'m1044','m1045','m1092','m1094'};
%   slow_learners = {'m1048','m1049','m1613','m1859','m1873'};
%   [colorSchemeTable, colorMapLookup] = build_learnerGroup_colorScheme( ...
%       fast_learners, slow_learners, ...
%       'saveDir', "Z:\Rodent Data\dualImaging_parkj\collectData", ...
%       'saveName', 'learnerGroupColorScheme.mat');
%
% See also: get_colors_for_animalIDs, plot_cvR2_acrossSessions

% -------- parse args --------
p = inputParser;
p.addParameter('hueFast', 0.36, @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
p.addParameter('hueSlow', 0.80, @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
p.addParameter('satRange', [0.55 0.85], @(x) isnumeric(x) && numel(x)==2);
p.addParameter('valRange', [0.45 0.85], @(x) isnumeric(x) && numel(x)==2);
p.addParameter('saveDir', '', @(x) ischar(x) || isstring(x));
p.addParameter('saveName', 'learnerGroupColorScheme.mat', @(x) ischar(x) || isstring(x));
p.parse(varargin{:});
opt = p.Results;

fast_learners = cellstr(fast_learners(:));
slow_learners = cellstr(slow_learners(:));

nFast = numel(fast_learners);
nSlow = numel(slow_learners);
assert(nFast>0 || nSlow>0, 'At least one of fast_learners/slow_learners must be non-empty.');

overlap = intersect(fast_learners, slow_learners);
assert(isempty(overlap), 'Animal(s) listed in both fast_learners and slow_learners: %s', strjoin(overlap, ', '));

% -------- generate per-group shade gradients --------
colorsFast = generate_group_shades_(opt.hueFast, nFast, opt.satRange, opt.valRange);
colorsSlow = generate_group_shades_(opt.hueSlow, nSlow, opt.satRange, opt.valRange);

allIDs  = [fast_learners; slow_learners];
allGrp  = [repmat({'fast'}, nFast, 1); repmat({'slow'}, nSlow, 1)];
allCols = [colorsFast; colorsSlow];

colorSchemeTable = table(allIDs, allGrp, allCols, 'VariableNames', {'animalID','group','color'});

colorMapLookup = containers.Map('KeyType','char', 'ValueType','any');
for i = 1:numel(allIDs)
    colorMapLookup(allIDs{i}) = allCols(i,:);
end

% -------- save (optional) --------
saveDir = string(opt.saveDir);
if strlength(saveDir) > 0
    if ~isfolder(saveDir)
        mkdir(saveDir);
    end
    saveName = char(opt.saveName);
    save(fullfile(saveDir, saveName), ...
        'colorSchemeTable', 'colorMapLookup', 'fast_learners', 'slow_learners', ...
        'opt');
    fprintf('Saved learner-group color scheme to %s\n', fullfile(saveDir, saveName));
end

end % function


% ===== helper: HSV gradient (dark+saturated -> light+desaturated) at fixed hue =====
function rgb = generate_group_shades_(hue, n, satRange, valRange)
if n == 0
    rgb = zeros(0,3);
    return;
end
if n == 1
    sVals = mean(satRange);
    vVals = mean(valRange);
else
    sVals = linspace(satRange(2), satRange(1), n);  % most saturated first
    vVals = linspace(valRange(1), valRange(2), n);  % darkest first
end
hsv = [repmat(hue, n, 1), sVals(:), vVals(:)];
rgb = hsv2rgb(hsv);
end
