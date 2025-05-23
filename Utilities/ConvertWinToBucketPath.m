function file_bucket = ConvertWinToBucketPath(file)

local_bucket = 'Z:\'; 
spock_bucket = '/jukebox/buschman/'; 

% Confirm not already converted
if ~isempty(regexp(file, regexptranslate('escape', local_bucket), 'once')) && isempty(regexp(file, regexptranslate('escape', spock_bucket), 'once'))
    file_bucket = [spock_bucket erase(file, local_bucket)];
    file_bucket = regexprep(file_bucket, '\\','/'); % convert backslashes to forward slashes
else
    file_bucket = file; 
end

end
