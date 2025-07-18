function fileList = findFilesWithPattern(baseDir, searchPattern)
%FINDFILESWITHPATTERN  Return full paths of all files that match a pattern.
%
%  fileList = findFilesWithPattern(baseDir, searchPattern)
%
%  INPUTS
%    baseDir       – directory to search (no recursion)
%    searchPattern – filename pattern with * wildcards, e.g.
%                    'm1873_041025*red*chunk*.mat'
%
%  OUTPUT
%    fileList      – cell array of full paths (char) matching the pattern.
%                    Returns {} if no file matches.
%
%  EXAMPLE
%    p = '/Volumes/buschman/Rodent Data/Wide Field Microscopy/ExampleData/Preprocessed/...
%         m1873_041025_task_day0_img_motif';
%    files = findFilesWithPattern(p,'m1873_041025*red*chunk*.mat');
%
%  NOTES
%    • Uses MATLAB's DIR, so all standard wildcard rules apply.
%    • No recursive descent; modify with **dir(fullfile(baseDir,'**',pattern))**
%      if you ever need that.
% ------------------------------------------------------------------------

arguments
    baseDir (1,1) string  {mustBeFolder}
    searchPattern (1,1) string
end

% Grab matches
info = dir(fullfile(baseDir, searchPattern));

% Build full paths
if isempty(info)
    fileList = {};
else
    fileList = fullfile({info.folder}, {info.name});
end
end
