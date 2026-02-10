function glmRezPath = find_latest_glmRez_file(filePath_mat, header)
%FIND_LATEST_GLMREZ_FILE  Find the most recent glmRez file for a session header.
%
% glmRezPath = find_latest_glmRez_file(filePath_mat, header)
%
% Looks for files of the form:
%   <header>_glmRez_MMDDYY.mat
% and returns the one with the latest date.
%
% EXCLUDES files like:
%   <header>_glmRez_region_*.mat
%
% INPUTS
%   filePath_mat : path to Matfiles directory
%   header       : session header string (e.g. 'm1045_122424')
%
% OUTPUT
%   glmRezPath   : full path to latest glmRez file
%                  '' if none found
%
% EXAMPLE
%   fp = find_latest_glmRez_file(filePath_mat, header);

% -------------------- sanity --------------------
if nargin < 2 || isempty(filePath_mat) || isempty(header)
    glmRezPath = '';
    return;
end

header = char(string(header));

% -------------------- list .mat files --------------------
d = dir(fullfile(filePath_mat, '*.mat'));
if isempty(d)
    glmRezPath = '';
    return;
end

fnames = {d.name};

% -------------------- regex: EXACT header + glmRez + MMDDYY --------------------
pat = ['^' regexptranslate('escape', header) '_glmRez_(\d{6})\.mat$'];

tok = regexp(fnames, pat, 'tokens');
isMatch = ~cellfun(@isempty, tok);

if ~any(isMatch)
    glmRezPath = '';
    return;
end

glmFiles = fnames(isMatch);
tok = tok(isMatch);

% -------------------- parse dates --------------------
mmddyy = cellfun(@(t) t{1}{1}, tok, 'UniformOutput', false);

try
    dt = datetime(mmddyy, 'InputFormat','MMddyy');
catch
    warning('find_latest_glmRez_file:DateParseFail', ...
        'Failed to parse MMDDYY dates in glmRez filenames.');
    glmRezPath = '';
    return;
end

% -------------------- pick latest --------------------
[~, idxLatest] = max(dt);
glmRezPath = fullfile(filePath_mat, glmFiles{idxLatest});

end
