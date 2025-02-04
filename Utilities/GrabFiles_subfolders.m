function [file_list, folder_list] = GrabFiles_subfolders(substring, searchdir)
% @Synopsis 
% Searches folder and all subfolders for any file with the defined substring.
% Outputs a list of file paths and folder paths. 

% @Inputs 
% @substring (required)
% A string containing the substring of the files to grab. Will grab any file 
% with substring anywhere in the file name; 
% Example 1: substring = '.pdf' will grab all pdfs in a folder/subfolders
% Example 2: substring = 'applesauce.mat' will grab all files that
% contain applesauce in the name applesauce irrespective of other 
% text in the file name. 
%
%
% NOTE ON Wildcards. Uses regular expressions to find matches so denote
% wildcard accordingly (e.g. \w*)
%
% @searchdir (required)
% Cell array of directories to search for files. 
%
% @Outputs: list of paths for each file and folder paths for each folder

% Validate input arguments

% Initialize output variables
file_list = [];
folder_list = [];
start = pwd;

if iscell(searchdir)
    searchdir = searchdir{1};
    warning('The input searchdir was a cell, will use the first entry!')
end


tempd = dir(searchdir);
isub = [tempd(:).isdir]; %returns logical vector
nameFolds = {tempd(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = []; %remove silly extra folders

%Loop through each subfolder in the dir.
for i = 1:length(nameFolds) %subfolder loop
    subfolder_fn = [searchdir filesep num2str(nameFolds{i})];
    [file_names,file_path,folder_names] = GrabFile(subfolder_fn,substring);
    file_list = [file_list, file_path];
    folder_list = [folder_list,folder_names];
end %Subfolderloop
cd(start)

% Remove duplicates if they occur
file_list = unique(file_list);
folder_list = unique(folder_list);

if ~isempty(file_list)
    % Reorder by the natural order of the file names (e.g. 1,2,3 not 1,10,100)
    [~, index] = natsortfiles(file_list);
    file_list = file_list(index);
end

if ~isempty(folder_list)
    % Reorder by the natural order of the folder names
    [~, index] = natsortfiles(folder_list);
    folder_list = folder_list(index);
end

end

function [in_fns, in_path, in_folder] = GrabFile(target_dir, substring)
if isempty(target_dir) 
    error('Target directory must be specified.');
end

startdir = pwd;
cd(target_dir);
files = dir(); 
matches = regexpi({files.name}, substring);
files = files(cellfun(@(x) ~isempty(x), matches));

if isempty(files)   % No files = throw a warning
    %warning('No files with "%s" in selected directory %s', substring, target_dir);
    in_fns = [];
    in_path = [];
    in_folder = [];
else
    in_fns = cell([1, length(files)]);
    in_path = cell([1, length(files)]);
    in_folder = cell([1, length(files)]);
    
    for cur_f = 1:length(files)
        in_fns{cur_f} = files(cur_f).name;
        in_folder{cur_f} = files(cur_f).folder;
        in_path{cur_f} = fullfile(files(cur_f).folder, files(cur_f).name);
    end

    in_fns = in_fns(~cellfun('isempty', in_fns));
    in_path = in_path(~cellfun('isempty', in_path));
    in_folder = in_folder(~cellfun('isempty', in_folder));

    if isempty(in_fns)
        %warning('No files with "%s" in selected directory %s', substring, target_dir);
    end
end
cd(startdir);
end
