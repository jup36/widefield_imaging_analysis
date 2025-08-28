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
