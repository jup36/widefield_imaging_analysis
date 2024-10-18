%% whereabouts
filePaths = {'/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA017/DA017_041624', ...
             '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA017/DA017_041824', ...
             '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA017/DA017_042224', ...
             '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA017/DA017_042924'}; 

% filePaths = {'/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA019/DA019_041924', ...
%              '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA019/DA019_2_042024', ...
%              '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA019/DA019_042224', ...
%              '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA019/DA019_042624'}; 

mean_stimAlignGoC = cell(length(filePaths), 2);
sem_stimAlignGoC = cell(length(filePaths), 1);

mean_hitFstLickDff = cell(length(filePaths), 1); 
sem_hitFstLickDff = cell(length(filePaths), 1); 

for ff = 1:length(filePaths)
    filePath = filePaths{ff};
    filePath_region = GrabFiles_sort_trials('_regionMask_SS', 0, {fullfile(filePath, 'Matfiles')});

    [~, header] = fileparts(filePath);
    headerParts = regexp(header, '_', 'split');  % Split the string at the underscore
    mouseIdC{ff,1} = headerParts{1};  % Part before the underscore
    dateC{ff,1} = headerParts{2};  % Part after the underscore
    figSaveDir = fullfile(fileparts(filePath), 'collectFigure');
    if exist(figSaveDir, 'dir')~=7
        mkdir(figSaveDir)
    end
    
    %load
    load(filePath_region{1}, 'rezSS'); 

    % cue-aligned go nogo
    mean_stimOnGoNogoSvmC_SS{ff, 1} = smooth2a(mean(rezSS.stimOnGoNogoSvm), 0, 4); 
    mean_stimOnGoNogoSvmC_SS{ff, 2} = rezSS.stimOnGoNogoSvmTs; 

    % cue-aligned hit FA
    mean_stimOnHitFaSvmC_SS{ff, 1} = smooth2a(mean(rezSS.stimOnHitFaSvm), 0, 1);
    mean_stimOnHitFaSvmC_SS{ff, 2} = rezSS.stimOnHitFaSvmTs; 

    % cue-aligned hit miss
    mean_stimOnHitMissSvmC_SS{ff, 1} = smooth2a(mean(rezSS.stimOnHitMissSvm), 0, 1);
    mean_stimOnHitMissSvmC_SS{ff, 2} = rezSS.stimOnHitMissSvmTs;   

    % cue-aligned hit cr
    mean_stimOnHitCrSvmC_SS{ff, 1} = smooth2a(mean(rezSS.stimOnHitCrSvm), 0, 1);
    mean_stimOnHitCrSvmC_SS{ff, 2} = rezSS.stimOnHitCrSvmTs; 

    % lick-aligned hit FA
    mean_lickHitFaSvmC_SS{ff, 1} = smooth2a(mean(rezSS.lickHitFaSvm), 0, 1); 
    mean_lickHitFaSvmC_SS{ff, 2} = rezSS.lickHitFaSvmTs; 

    fprintf("processed file #%d\n", ff)
end

%% plot cue-aligned hit FA svm
% color
learningColor = slanCM('blues', 5); 
learningColor = learningColor(2:5, :); 
learningColorC = mat2cell(learningColor, ones(size(learningColor, 1), 1), size(learningColor, 2));  
id_dateC = cellfun(@(a, b) strcat(a, '_', b), mouseIdC, dateC, 'UniformOutput', false); 

%% plot cue-aligned gng svm
plotMeanColorC(cell2mat(mean_stimOnGoNogoSvmC_SS(:, 1)), mean_stimOnGoNogoSvmC_SS{1, 2}, learningColorC, id_dateC)
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -5:0.05:5); xlim([-1 5]); ylim([0.45 0.75])
print(fullfile(figSaveDir, strcat(mouseIdC{1}, '_dff_cueOnAligned_GoNogo_SVM_SS')), '-dpdf', '-vector', '-bestfit')

%% plot cue-aligned hit fa svm
plotMeanColorC(cell2mat(mean_stimOnHitFaSvmC_SS(:, 1)), mean_stimOnHitFaSvmC_SS{1, 2}, learningColorC, id_dateC)
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -5:0.05:5); xlim([-1 5]); %ylim([-0.4 1.1])
print(fullfile(figSaveDir, strcat(mouseIdC{1}, '_dff_cueOnAligned_HitFa_SVM_SS')), '-dpdf', '-vector', '-bestfit')

%% plot cue-aligned hit miss svm
plotMeanColorC(cell2mat(mean_stimOnHitMissSvmC_SS(:, 1)), mean_stimOnHitMissSvmC_SS{1, 2}, learningColorC, id_dateC)
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -5:0.05:5); xlim([-1 5]); %ylim([-0.4 1.1])
print(fullfile(figSaveDir, strcat(mouseIdC{1}, '_dff_cueOnAligned_HitMiss_SVM_SS')), '-dpdf', '-vector', '-bestfit')

%% plot cue-aligned hit CR svm
plotMeanColorC(cell2mat(mean_stimOnHitCrSvmC_SS(:, 1)), mean_stimOnHitCrSvmC_SS{1, 2}, learningColorC, id_dateC)
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -5:0.05:5); xlim([-1 5]); %ylim([-0.4 1.1])
print(fullfile(figSaveDir, strcat(mouseIdC{1}, '_dff_cueOnAligned_HitCr_SVM_SS')), '-dpdf', '-vector', '-bestfit')

%% plot lick-aligned hit FA svm
plotMeanColorC(cell2mat(mean_lickHitFaSvmC_SS(:, 1)), mean_lickHitFaSvmC_SS{1, 2}, learningColorC, id_dateC)
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -5:0.05:5); xlim([-1 1]); %ylim([-0.4 1.1])
print(fullfile(figSaveDir, strcat(mouseIdC{1}, '_dff_lickAligned_HitFA_SVM_SS')), '-dpdf', '-vector', '-bestfit')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dffTsItp_mean, dffTsItp_sem, dffTsItp, dffTsItpTs] = getInterpolatedMeanSemDffFromDffTsCell(dffTsCell)
dffTsC = cellWithNonEmptyColumns(dffTsCell);
dffTsLinC = cellfun(@(a) linspace(-1, 1, length(a)), dffTsC(:, 2), 'UniformOutput', false);
[dffTsItp, dffTsItpTs] = temporalAlignInterp1(dffTsC(:, 1), dffTsLinC, 0.001);
[dffTsItp_mean, ~, dffTsItp_sem] = meanstdsem(cell2mat(dffTsItp));

end




