function [svm] = dffPostprocess_auditory_gng_dual_classfier(filePath, channel, timeWin, step)
%filePath = '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1045_jRGECO_GRABda/m1045_121324/task';

%% load trial-by-trial behavior & task data tbytDat.mat
%parentDir = fileparts(filePath);
match_header = regexp(filePath, 'm\d{1,4}_\d{6}', 'match');
header = match_header{1};

if strcmpi(channel, 'green')
    fileBeh = GrabFiles_sort_trials('_green_tbytDat_dff', 0, {fullfile(filePath, 'Matfiles')});
elseif strcmpi(channel, 'red')
    fileBeh = GrabFiles_sort_trials('_red_tbytDat_dff', 0, {fullfile(filePath, 'Matfiles')});
end
load(fullfile(fileBeh{1}), 'tbytDat')

% load the preprocessed dffs collect
if strcmpi(channel, 'green')
    load(fullfile(filePath, 'Matfiles', strcat(header, "_green_dff_smCollect.mat")), 'dffsmCell');
elseif strcmpi(channel, 'red')
    load(fullfile(filePath, 'Matfiles', strcat(header, '_red_dff_smCollect.mat')), 'dffsmCell');
end

% task event indices
taskI = taskEventIndices(tbytDat); 

%% Temporal alignment and interpolation
for t = 1:length(tbytDat)
    frT = tbytDat(t).frameTrel; % store timestamps
    [dffC{t}, dffCtF{t}] = alignToEventInterp3D(dffsmCell{t}, 0, frT, timeWin, step); 
end
svm.timePeth = dffCtF{1}; 

%% Train/Test SVM classifier
[svm.accuracy, svm.betaValues, svm.finalBetaValues] = trainDffClassifierC(dffC(taskI.goI), dffC(taskI.nogoI), 10, 5);

%% Save results
saveName = [header, '_', channel, '_svmRez.mat']; 
save(fullfile(filePath,'Matfiles', saveName), 'svm');

end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dffOut, tF] = alignToEventInterp3D(dffTs, eventTimeToAlign, frameT, timeWin, step)
% NOTE: Modified the function to return outcomes even if the timeWin goes
% out of bounds of the frame time! (2/4/2025, Junchol Park)
% e.g., step = 0.05; % 50ms

alignedFrameT = frameT + eventTimeToAlign;  
tF = timeWin(1):step:timeWin(end); 

% Define the valid range based on frameTOnTf
validTfIdx = tF >= min(alignedFrameT) & tF <= max(alignedFrameT);
tFval = tF(validTfIdx); 

% Get dimensions
[rows, cols, ~] = size(dffTs);

% Initialize output with NaNs 
dffOut = NaN(size(dffTs,1),size(dffTs,2),length(tF)); 

% Loop through each pixel
for i = 1:rows
    for j = 1:cols
        % Extract time series for this pixel
        pixelSeries = squeeze(dffTs(i, j, :));
        
        % Apply strict NaN criterion: Only interpolate if ALL time points are valid
        if all(~isnan(pixelSeries))
            dffOut(i, j, validTfIdx) = interp1(alignedFrameT, pixelSeries, tFval, 'pchip');
        end
    end
end
end

function index = taskEventIndices(tbytDat) 
%This utility function generates various task event indices from the tbytDat structure. 

index.waterI = cellfun(@(a) ~isempty(a), {tbytDat.water}); 
index.lickI = cellfun(@(a) ~isempty(a), {tbytDat.Lick}); 
index.airpuffI = cellfun(@(a) ~isempty(a), {tbytDat.airpuff}); 
index.goI = [tbytDat.rewardTrI]'==1; 
index.nogoI = [tbytDat.punishTrI]'==1; 
index.hitI = cellfun(@(a) ~isempty(a), {tbytDat.hitLicks})'; 
index.missI = cell2mat({tbytDat.rewardTrI})' & cellfun(@(a) isempty(a), {tbytDat.water})'; 
index.faI = cellfun(@(a) ~isempty(a), {tbytDat.faLicks})'; 
index.crI = cell2mat({tbytDat.punishTrI})' & cellfun(@(a) isempty(a), {tbytDat.airpuff})'; 
end