function fullPath = findFolderPattern(baseFolder, pattern)
% FIND FOLDER in a given directory that matches a given pattern exactly.
%
% Inputs:
%   baseFolder - The parent directory to search within
%   pattern    - A string with '*' wildcards (e.g., 'dffG*trial_1') or
%                a cell array of substrings (e.g., {'dffG', 'trial_1'}).
%
% Output:
%   fullPath - Full path to the matched folder or empty if not found.

% List all subdirectories
dirs = dir(baseFolder);

% Filter out non-directories
dirs = dirs([dirs.isdir]);

% Remove '.' and '..' from results
dirs = dirs(~ismember({dirs.name}, {'.', '..'}));

% Extract folder names
folderNames = {dirs.name};

if ischar(pattern)
    % Convert wildcard '*' to regex '.*'
    regexPattern = strrep(pattern, '*', '.*');

    % Ensure exact number matching for digits (avoid partial matches like 'trial_11')
    regexPattern = regexprep(regexPattern, '(\d+)', '$1(?!\d)');

    % Perform regex match
    matchIdx = ~cellfun(@isempty, regexp(folderNames, ['^' regexPattern '$'], 'once'));

elseif iscell(pattern)
    % Match substrings with exact number matching
    matchIdx = cellfun(@(name) all(contains(name, pattern)) && ...
        all(cellfun(@(p) ~contains(name, [p, '0':'9']), pattern)), folderNames);

else
    error('Pattern must be a string (wildcard) or a cell array of substrings.');
end

% Return the matched folder path
if any(matchIdx)
    fullPath = fullfile(baseFolder, folderNames{matchIdx});
else
    fullPath = ''; % Return empty if no match found
end
end