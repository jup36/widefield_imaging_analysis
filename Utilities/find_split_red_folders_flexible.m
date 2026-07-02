function redFolderC = find_split_red_folders_flexible(filePath_img, header, modeKeyword)
% Find split red folders while tolerating naming variants.
%
% Supports examples:
%   m1045_122424_baseline_1_red
%   m1613_042425_baseline1_red
%   m1045_122424_task_day4-8_1_red
%
% Inputs:
%   filePath_img  - parent image folder
%   header        - e.g. m1613_042425
%   modeKeyword   - e.g. 'baseline' or 'task'

redFolderC = {};

if exist(filePath_img, 'dir') ~= 7
    warning('Image folder does not exist:\n%s', filePath_img);
    return
end

d = dir(filePath_img);
d = d([d.isdir]);

folderNames = {d.name};
folderNames = folderNames(~ismember(folderNames, {'.', '..'}));

folderPaths = cellfun(@(x) fullfile(filePath_img, x), folderNames, 'UniformOutput', false);

% Flexible pattern:
%
% header + anything + modeKeyword + optional underscore + number + _red
%
% This matches:
%   m1613_042425_baseline_1_red
%   m1613_042425_baseline1_red
%
% Also works for task if modeKeyword = 'task':
%   m1045_122424_task_day4-8_1_red

expr = ['^' regexptranslate('escape', header) ...
        '.*' regexptranslate('escape', modeKeyword) ...
        '.*_?\d+_red$'];

matchI = cellfun(@(x) ~isempty(regexp(x, expr, 'once')), folderNames);

redFolderC = folderPaths(matchI);

% Fallback: any direct child folder ending with optional underscore-number-red
% and containing modeKeyword.
if isempty(redFolderC)
    exprFallback = ['.*' regexptranslate('escape', modeKeyword) '.*_?\d+_red$'];
    matchI = cellfun(@(x) ~isempty(regexp(lower(x), lower(exprFallback), 'once')), folderNames);
    redFolderC = folderPaths(matchI);
end

try
    redFolderC = sort_nat(redFolderC);
catch
    redFolderC = sort(redFolderC);
end

end