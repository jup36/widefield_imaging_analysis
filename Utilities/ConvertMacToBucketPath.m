function file_bucket = ConvertMacToBucketPath(file)

local_bucket = '/Volumes/buschman/'; 
spock_bucket = '\jukebox\buschman\'; 

%confirm not already converted
if ~isempty(regexp(file,local_bucket,'once')) && isempty(regexp(file,spock_bucket,'once'))
    file_bucket = [spock_bucket erase(file,local_bucket)];
    file_bucket = regexprep(file_bucket, '\','/');
else
    file_bucket = file; 
end


end