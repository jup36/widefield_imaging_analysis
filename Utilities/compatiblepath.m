function outPath = compatiblepath(inPath)
% compatiblepath converts between Mac and Windows paths based on OS
% If running on Windows, it converts Mac path to Windows format (if needed)
% If running on Mac, it converts Windows path to Mac format (if needed)
%
% Input:
%   inPath - string, original file path (Mac or Windows)
% Output:
%   outPath - string, converted path compatible with current OS

if ispc && contains(inPath, '/Volumes/buschman/')
    outPath = ConvertMacToWinPath(inPath);
elseif ismac && contains(inPath, 'Z:\')
    outPath = ConvertWinToMacPath(inPath);
else
    % If already in the correct format or OS is unsupported, return as-is
    outPath = inPath;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function winPath = ConvertMacToWinPath(macPath)
% ConvertMacToWinPath converts a Mac-style file path to a Windows-style path
% Specific to '/Volumes/buschman/' → 'Z:\'
%
% Input:
%   macPath - string, Mac-style path
% Output:
%   winPath - string, Windows-style path

if contains(macPath, '/Volumes/buschman/')
    % Replace only the full prefix with correct Windows drive path
    winPath = strrep(macPath, '/Volumes/buschman/', 'Z:\');
    % Replace remaining slashes with backslashes
    winPath = strrep(winPath, '/', '\');
else
    warning('The input path is not a proper Mac path!');
    winPath = macPath;
end
end

function macPath = ConvertWinToMacPath(winPath)
% ConvertWinToMacPath converts a Windows-style path to a Mac-style path
% Specific to 'Z:\' → '/Volumes/buschman/'
%
% Input:
%   winPath - string, Windows-style path
% Output:
%   macPath - string, Mac-style path

if contains(winPath, 'Z:\')
    % Replace only the full prefix with correct Mac volume path
    macPath = strrep(winPath, 'Z:\', '/Volumes/buschman/');
    % Replace remaining backslashes with forward slashes
    macPath = strrep(macPath, '\', '/');
else
    warning('The input path is not a proper Windows path!');
    macPath = winPath;
end

end

end