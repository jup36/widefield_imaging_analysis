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
