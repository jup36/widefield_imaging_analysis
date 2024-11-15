function file_bucket = ConvertMacToBucketPath(file)

gp = general_params_mac;

%confirm not already converted
if ~isempty(regexp(file,gp.local_bucket,'once')) && isempty(regexp(file,gp.spock_bucket,'once'))
    file_bucket = [gp.spock_bucket erase(file,gp.local_bucket)];
    file_bucket = regexprep(file_bucket, '\','/');
else
    file_bucket = file; 
end


end