function file_win = ConvertBucketToWinPath(file)

spock_bucket = '/jukebox/buschman/';
win_bucket   = 'Z:\';

% If path starts with spock bucket, replace prefix
if startsWith(file, spock_bucket)
    file_win = strrep(file, spock_bucket, win_bucket);
else
    file_win = file;
end

% Replace all forward slashes with backslashes
file_win = strrep(file_win, '/', '\');

end