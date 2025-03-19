function file_mac = ConvertBucketToMacPath(file)

local_bucket = '/Volumes/buschman/'; 
spock_bucket = '/jukebox/buschman/'; % Corrected bucket path format

% Check if the path starts with spock_bucket and is not already converted
if startsWith(file, spock_bucket) && ~startsWith(file, local_bucket)
    % Replace the spock_bucket prefix with local_bucket
    file_mac = strrep(file, spock_bucket, local_bucket);
else
    file_mac = file; % If already converted, return as is
end

end

