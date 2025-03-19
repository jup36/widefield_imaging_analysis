mask = load('/Users/jp3025/Documents/codes/Widefield_Imaging_Analysis/Preprocessing/brainoutline.mat'); 
mask_540 = mask.brainoutline; 

% Assume your original mask is named `mask_540`
mask_64 = imresize(mask_540, [64 64], 'nearest');

% Ensure the result is still logical
mask_64 = logical(mask_64);

save(fullfile('/Users/jp3025/Documents/codes/Widefield_Imaging_Analysis/Preprocessing', 'brainoutline.mat'), 'mask_64'); 



