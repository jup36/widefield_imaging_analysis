% Define paths
filePath_basisMotifs = '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1045_jRGECO_GRABda/motif_figures'; 
filePath_save = '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1045_jRGECO_GRABda/motif_figures/concat_motif_images'; 

% Get list of PNG files
file_list_png = GrabFiles_sort_trials('_Group', 0, {filePath_basisMotifs}); 

% Ensure save directory exists
if ~exist(filePath_save, 'dir')
    mkdir(filePath_save);
end

for f = 1:12
    substring = sprintf("Group%d.png", f); 
    
    % Find indices of matching files
    idx = find(cellfun(@(name) ~isempty(regexp(name, substring, 'once')), file_list_png));
    
    if length(idx) < 26
        warning('Less than 26 images found for %s. Skipping...', substring);
    end
    
    % Load all images
    images = cell(1, length(idx));
    for ff = 1:length(idx)
        images{ff} = imread(file_list_png{idx(ff)});
    end
    
    % Select 10 images from the middle of the 27
    total_images = length(images);
    mid_start = floor((total_images - 10) / 2) + 1; % Compute middle start index
    selected_images = images(mid_start:mid_start+9);
    
    % Concatenate selected images vertically
    concatenated_image = cat(2, selected_images{:});
    
    % Define output filename
    save_filename = fullfile(filePath_save, sprintf('concatenated_motif_%s.png', substring));
    
    % Save concatenated image
    imwrite(concatenated_image, save_filename);
    
    fprintf('Saved: %s\n', save_filename);
end
