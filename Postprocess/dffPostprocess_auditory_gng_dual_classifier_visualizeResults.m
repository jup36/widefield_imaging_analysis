figSaveDir = '/Volumes/buschman/Rodent Data/dualImaging_parkj/collectFigure'; 

% Session 12/13/24
svm_red_121324 = load('/Volumes/buschman/Rodent Data/dualImaging_parkj/m1045_jRGECO_GRABda/m1045_121324/task/Matfiles/m1045_121324_red_svmRez.mat', 'svm');
svm_green_121324 = load('/Volumes/buschman/Rodent Data/dualImaging_parkj/m1045_jRGECO_GRABda/m1045_121324/task/Matfiles/m1045_121324_green_svmRez.mat', 'svm');
svm_red_121324 = svm_red_121324.('svm'); 
svm_green_121324 = svm_green_121324.('svm'); 

% Session 12/24/24
svm_red_122424 = load('/Volumes/buschman/Rodent Data/dualImaging_parkj/m1045_jRGECO_GRABda/m1045_122424/task/Matfiles/m1045_122424_red_svmRez.mat', 'svm');
svm_green_122424 = load('/Volumes/buschman/Rodent Data/dualImaging_parkj/m1045_jRGECO_GRABda/m1045_122424/task/Matfiles/m1045_122424_green_svmRez.mat', 'svm');
svm_red_122424 = svm_red_122424.('svm'); 
svm_green_122424 = svm_green_122424.('svm'); 

%% Visualize
timeX = -0.9:0.05:4;
%timeXI = [timeX', 1:length(timeX)]; 

% jRGECO1a Early vs Late Go-NoGo discriminability (Linear SVM classifier)
figure; hold on; 
plot(timeX, squeeze(mean(svm_red_121324.accuracy,2)));
plot(timeX, squeeze(mean(svm_red_122424.accuracy,2)));
set(gca, "TickDir", 'out')
print(fullfile(figSaveDir, 'm1045_red_gonogo_discriminability_svm_earlylate'), '-dpdf', '-bestfit', '-vector')

% GRABDA Early vs Late Go-NoGo discriminability (Linear SVM classifier)
figure; hold on; 
plot(timeX, squeeze(mean(svm_green_121324.accuracy,2)));
plot(timeX, squeeze(mean(svm_green_122424.accuracy,2)));
set(gca, "TickDir", 'out')
print(fullfile(figSaveDir, 'm1045_green_gonogo_discriminability_svm_earlylate'), '-dpdf', '-bestfit', '-vector')

% bvalues
%figure; imagesc(applyImgaussfiltSd(svm_red_122424.finalBetaValues{80}, 2)); 

%% Analyze beta values 
betaRez_svm_red_122424 = analyzeSVMBetaValues('/Volumes/buschman/Rodent Data/dualImaging_parkj/m1045_jRGECO_GRABda/m1045_122424/task', ...
    svm_red_122424.finalBetaValues, 'red'); 

plot(timeX, betaRez_svm_red_122424.M2)

betaRez_svm_green_122424 = analyzeSVMBetaValues('/Volumes/buschman/Rodent Data/dualImaging_parkj/m1045_jRGECO_GRABda/m1045_122424/task', ...
    svm_green_122424.finalBetaValues, 'green'); 

betaVal_rez = analyzeSVMBetaValues('/Volumes/buschman/Rodent Data/dualImaging_parkj/m1045_jRGECO_GRABda/m1045_121324/task', ...
    svm_red_121324.finalBetaValues, 'red'); 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rez = analyzeSVMBetaValues(filePath, betaValueC, channel)

% Define the normalization function using nanmean and nanstd
normZ = @(a) (a - nanmean(a(:))) ./ nanstd(a(:));

% Apply the normalization function to each cell in the cell array
betaValueCnorm = cellfun(normZ, betaValueC, 'UniformOutput', false);
betaValueStack = cat(3, betaValueCnorm{:}); % convert the cell array with beta coefficients to a 3D stack

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


