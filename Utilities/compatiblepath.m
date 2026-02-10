function outPath = compatiblepath(inPath)
% compatiblepath converts between Mac and Windows paths based on OS
% Robust to slash direction and fullfile() behavior

inPath = char(inPath);  % ensure char
normPath = strrep(inPath, '\', '/');  % normalize to forward slashes

if ispc
    % Mac → Windows
    if startsWith(normPath, '/Volumes/buschman/', 'IgnoreCase', true)
        outPath = ConvertMacToWinPath(normPath);
    else
        outPath = inPath;
    end

elseif ismac
    % Windows → Mac
    if startsWith(normPath, 'Z:/', 'IgnoreCase', true)
        outPath = ConvertWinToMacPath(normPath);
    else
        outPath = inPath;
    end

else
    outPath = inPath;
end
end

% -------------------------------------------------------------------------
function winPath = ConvertMacToWinPath(macPath)
% '/Volumes/buschman/...' → 'Z:\...'

winPath = strrep(macPath, '/Volumes/buschman/', 'Z:/');
winPath = strrep(winPath, '/', '\');
end

% -------------------------------------------------------------------------
function macPath = ConvertWinToMacPath(winPath)
% 'Z:\...' or 'Z:/...' → '/Volumes/buschman/...'

winPath = strrep(winPath, '\', '/');
macPath = strrep(winPath, 'Z:/', '/Volumes/buschman/');
end
