function fullPaths = findFilePattern(baseFolder, pattern)
% FIND FILES in a given directory that match a given pattern exactly.
%
% Inputs:
%   baseFolder - The parent directory to search within
%   pattern    - A string with '*' wildcards (e.g., 'data_*.mat') or
%                a cell array of substrings (e.g., {'data', 'session1'}).
%
% Output:
%   fullPaths - A cell array of full paths to the matched files, or empty if not found.

% List all files (excluding folders)
files = dir(fullfile(baseFolder, '*'));  % Get all files and folders
files = files(~[files.isdir]);  % Keep only files

% Extract file names
fileNames = {files.name};

if ischar(pattern)  % Case: Wildcard string pattern
    % Convert wildcard '*' to regex '.*' (matches any characters)
    regexPattern = strrep(pattern, '*', '.*');
    
    % Ensure exact number matching (avoid partial matches like 'trial_11')
    regexPattern = regexprep(regexPattern, '(\d+)', '$1(?!\d)');
    
    % Perform regex match
    matchIdx = ~cellfun(@isempty, regexp(fileNames, ['^' regexPattern '$'], 'once'));

elseif iscell(pattern)  % Case: Cell array of substrings
    % Match only if all substrings exist AND numbers match exactly
    matchIdx = cellfun(@(name) all(contains(name, pattern)) && ...
        all(cellfun(@(p) isempty(regexp(name, [p '\d+'], 'once')), pattern)), fileNames);

else
    error('Pattern must be a string (wildcard) or a cell array of substrings.');
end

% Return full paths to all matched files
matchedFiles = fileNames(matchIdx);
if ~isempty(matchedFiles)
    fullPaths = fullfile(baseFolder, matchedFiles);
    fullPaths = natsortfiles(fullPaths); 
else
    fullPaths = {}; % Return empty if no match found
end
end
