function rez = analyzeSVMBetaValues(filePath, betaValueC, channel)

% Define the normalization function using nanmean and nanstd
normZ = @(a) (a - nanmean(a(:))) ./ nanstd(a(:));

% Apply the normalization function to each cell in the cell array
betaValueCnorm = cellfun(normZ, betaValueC, 'UniformOutput', false);
betaValueStack = cat(3, betaValueCnorm{:}); 

% load Allen dorsalMap
load(fullfile('/Users/jp3025/Documents/codes/WidefieldAnalysis_Musall/allenDorsalMap.mat'), 'dorsalMaps', 'motorMask', 'smotorMask', 'ssMask', 'vMask', 'rsMask');

% load Allen transformation object (the alignment parameters to Allen CCF)
folder_list_imgTrial = GrabFiles_sort_trials('_img', 0, {filePath}); 
if strcmpi(channel, 'green')
    path_transParamsAllen = find_keyword_containing_files(folder_list_imgTrial{1}, 'transParamsAllen_green', 'recursive', true); 
elseif strcmpi(channel, 'red')
    path_transParamsAllen = find_keyword_containing_files(folder_list_imgTrial{1}, 'transParamsAllen_red', 'recursive', true); 
end
load(fullfile(path_transParamsAllen{1}), 'transParams')

betaValueStackAlign = alignStackToAllenKabsch(betaValueStack, dorsalMaps.dorsalMap, transParams.tformObj); % row x col x frame

% apply region masks
rez.M1 = apply2DMaskTo3DStack(betaValueStackAlign, motorMask);
rez.M2 = apply2DMaskTo3DStack(betaValueStackAlign, smotorMask);
rez.Ss = apply2DMaskTo3DStack(betaValueStackAlign, ssMask);
rez.V1 = apply2DMaskTo3DStack(betaValueStackAlign, vMask);
rez.Rs = apply2DMaskTo3DStack(betaValueStackAlign, rsMask);

end
