function colorMap = get_colors_for_animalIDs(animalIDs, colorSchemeSource, varargin)
%GET_COLORS_FOR_ANIMALIDS  Look up RGB colors for a specific list of animal IDs.
%
% SYNOPSIS
%   colorMap = get_colors_for_animalIDs(animalIDs, colorSchemeSource, ...)
%
% DESCRIPTION
%   Applies a color scheme (built by build_learnerGroup_colorScheme.m,
%   loaded fresh from its saved .mat file, or passed directly in memory)
%   to a specific animalIDs list -- e.g. the animalIDs returned by
%   collect_cvR2_perAnimal.m for a given analysis. Animals not found in
%   the color scheme fall back to a neutral gray (with a warning), so a
%   new/unclassified animal doesn't error out the whole plot.
%
% INPUTS
%   animalIDs         : cellstr of animal IDs for THIS analysis, in the
%                       order you want colors returned (e.g. matching the
%                       order of glmCvR2C_collect).
%   colorSchemeSource : any of:
%                         - containers.Map (animalID -> 1x3 RGB), e.g. the
%                           colorMapLookup output of build_learnerGroup_colorScheme.m
%                         - a table with columns 'animalID' and 'color'
%                           (e.g. colorSchemeTable from the same function)
%                         - a char/string path to a .mat file containing
%                           a saved 'colorMapLookup' variable
%
% NAME-VALUE ARGS
%   'fallbackColor' : 1x3 RGB used for any animalID not found in the color
%                     scheme. Default: [0.6 0.6 0.6] (gray).
%   'warnOnFallback': logical, warn when a fallback color is used.
%                     Default: true.
%
% OUTPUT
%   colorMap : [numel(animalIDs) x 3] RGB matrix, aligned to animalIDs.
%
% EXAMPLE
%   colorMap = get_colors_for_animalIDs(animalIDs, ...
%       "Z:\Rodent Data\dualImaging_parkj\collectData\learnerGroupColorScheme.mat");
%   h = plot_cvR2_acrossSessions(glmCvR2C_collect, animalIDs, 'colorMap', colorMap);
%
% See also: build_learnerGroup_colorScheme, plot_cvR2_acrossSessions

p = inputParser;
p.addParameter('fallbackColor', [0.6 0.6 0.6], @(x) isnumeric(x) && numel(x)==3);
p.addParameter('warnOnFallback', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;

animalIDs = cellstr(animalIDs(:));

% -------- resolve colorSchemeSource into a containers.Map --------
if ischar(colorSchemeSource) || isstring(colorSchemeSource)
    matPath = char(colorSchemeSource);
    assert(exist(matPath,'file')==2, 'Color scheme file not found: %s', matPath);
    S = load(matPath, 'colorMapLookup');
    assert(isfield(S,'colorMapLookup'), '%s does not contain a colorMapLookup variable.', matPath);
    lookupMap = S.colorMapLookup;
elseif isa(colorSchemeSource, 'containers.Map')
    lookupMap = colorSchemeSource;
elseif istable(colorSchemeSource)
    assert(all(ismember({'animalID','color'}, colorSchemeSource.Properties.VariableNames)), ...
        'Table must contain animalID and color columns.');
    lookupMap = containers.Map('KeyType','char','ValueType','any');
    for i = 1:size(colorSchemeSource,1)
        lookupMap(char(colorSchemeSource.animalID{i})) = colorSchemeSource.color(i,:);
    end
else
    error('colorSchemeSource must be a containers.Map, a table, or a path to a saved .mat file.');
end

% -------- look up each requested animal --------
nAnimals = numel(animalIDs);
colorMap = nan(nAnimals, 3);
missing  = {};

for i = 1:nAnimals
    aid = animalIDs{i};
    if isKey(lookupMap, aid)
        colorMap(i,:) = lookupMap(aid);
    else
        colorMap(i,:) = opt.fallbackColor;
        missing{end+1} = aid; %#ok<AGROW>
    end
end

if ~isempty(missing) && opt.warnOnFallback
    warning('get_colors_for_animalIDs:fallback', ...
        'No color scheme entry for: %s -- using fallback gray.', strjoin(missing, ', '));
end

end % function
