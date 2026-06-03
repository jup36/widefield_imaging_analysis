function file_bucket = ConvertWinToBucketPath(file)
% ConvertWinToBucketPath converts Windows-mounted bucket paths to Scotty paths.
%
% Example:
%   Z:\Rodent Data\... 
% becomes:
%   /jukebox/buschman/Rodent Data/...

local_bucket = "Z:\";
spock_bucket = "/jukebox/buschman/";

% Preserve input type
inputWasChar = ischar(file);

% Work internally as string scalar
file = string(file);

% Convert only if it starts with local bucket and is not already converted
if startsWith(file, local_bucket) && ~startsWith(file, spock_bucket)

    % Remove Z:\ prefix
    relative_path = erase(file, local_bucket);

    % Convert backslashes to forward slashes
    relative_path = replace(relative_path, "\", "/");

    % Concatenate as string scalar
    file_bucket = spock_bucket + relative_path;

else
    file_bucket = file;
end

% Return char if input was char
if inputWasChar
    file_bucket = char(file_bucket);
end

end