%
filePath = '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1859_jRGECO_GRABda/m1859_051925/task'; 
filePathImg = find_keyword_containing_folder(filePath, '_img', 'recursive', false); 

[~, folder_list_dff_green] = GrabFiles_subfolders('_green_ome.tiff', filePathImg);
fileStack_green = GrabFiles_sort_trials('ome_stack.mat', 0, folder_list_dff_green(1)); 
load(fileStack_green{1}, 'stack')
wfImage_green = nanmean(stack, 3); %clear stack


[file_list_first_stack, folder_list_raw] = GrabFiles_subfolders('_green_ome.tiff', filePathImg); % use GrabFiles_sort_trials to sort both files and folders
ref_img = GetReferenceImage(file_list_first_stack{1}, 'first');
