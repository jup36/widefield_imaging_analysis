function glmRezPath = find_latest_glmRez_file(filePath_mat, header, keyword)
%FIND_LATEST_GLMREZ_FILE  Find latest glmRez file using header + keyword.
%
% glmRezPath = find_latest_glmRez_file(filePath_mat, header)
% glmRezPath = find_latest_glmRez_file(filePath_mat, header, keyword)
%
% Looks for .mat files containing:
%   <header>_<keyword>_<MMDDYY>.mat
%
% Example filename:
%   m1045_122424_glmRez_redCal_L10K10_070126.mat
%
% Example call:
%   glmRezPath = find_latest_glmRez_file( ...
%       filePath_mat, header, 'glmRez_redCal_L10K10');
%
% If multiple files match the same keyword, returns the one with the latest
% MMDDYY date label.
%
% INPUTS
%   filePath_mat : path to Matfiles directory
%   header       : session header string, e.g. 'm1045_122424'
%   keyword      : user keyword, e.g.
%                    'glmRez_redCal_L10K10'
%                    'glmRez_greenDA_L10K10'
%
% OUTPUT
%   glmRezPath   : full path to latest matching glmRez file
%                  '' if none found

%% -------------------- defaults / sanity --------------------

glmRezPath = '';

if nargin < 2 || isempty(filePath_mat) || isempty(header)
    return
end

if nargin < 3 || isempty(keyword)
    keyword = 'glmRez';
end

filePath_mat = char(string(filePath_mat));
header       = char(string(header));
keyword      = char(string(keyword));

if exist(filePath_mat, 'dir') ~= 7
    warning('find_latest_glmRez_file:MissingDir', ...
        'Matfiles directory does not exist:\n%s', filePath_mat);
    return
end

%% -------------------- list .mat files --------------------

d = dir(fullfile(filePath_mat, '*.mat'));

if isempty(d)
    return
end

fnames = {d.name};

%% -------------------- match header + keyword + final date --------------------
%
% This intentionally allows extra characters between keyword and the date.
%
% Matches:
%   m1045_122424_glmRez_redCal_L10K10_070126.mat
%   m1045_122424_glmRez_greenDA_L10K10_070126.mat
%
% Also matches if there is extra text:
%   m1045_122424_glmRez_redCal_L10K10_extra_070126.mat
%
% Date must be the final 6 digits before .mat.

pat = ['^' regexptranslate('escape', header) ...
       '.*' regexptranslate('escape', keyword) ...
       '.*_(\d{6})\.mat$'];

tok = regexp(fnames, pat, 'tokens', 'once');
isMatch = ~cellfun(@isempty, tok);

if ~any(isMatch)
    fprintf('No glmRez file found for header "%s" with keyword "%s" in:\n%s\n', ...
        header, keyword, filePath_mat);
    return
end

glmFiles = fnames(isMatch);
tok = tok(isMatch);

%% -------------------- parse final MMDDYY dates --------------------

mmddyy = cellfun(@(t) t{1}, tok, 'UniformOutput', false);

try
    dt = datetime(mmddyy, 'InputFormat', 'MMddyy');
catch
    warning('find_latest_glmRez_file:DateParseFail', ...
        'Failed to parse MMDDYY dates in matching glmRez filenames.');
    return
end

%% -------------------- pick latest date --------------------

[~, idxLatest] = max(dt);

glmRezPath = fullfile(filePath_mat, glmFiles{idxLatest});

fprintf('Selected latest glmRez file:\n%s\n', glmRezPath);

end