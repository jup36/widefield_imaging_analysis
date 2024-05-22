%% whereabouts
filePath1 = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA017/DA017_041624'; 
[~, header1] = fileparts(filePath1);
header1Parts = regexp(header1, '_', 'split');  % Split the string at the underscore
mouseId1 = header1Parts{1};  % Part before the underscore
date1 = header1Parts{2};  % Part after the underscore

filePath2 = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA017/DA017_042924'; 
[~, header2] = fileparts(filePath2);
header2Parts = regexp(header2, '_', 'split');  % Split the string at the underscore
mouseId2 = header2Parts{1};  % Part before the underscore
date2 = header2Parts{2};  % Part after the underscore

figSaveDir = fullfile(fileparts(filePath1), 'collectFigure'); 
if exist(figSaveDir, 'dir')~=7
    mkdir(figSaveDir)
end

%% load and prep data 
RSPath1 = GrabFiles_sort_trials('_regionMask_RS', 0, {fullfile(filePath1, 'Matfiles')});
RSPath2 = GrabFiles_sort_trials('_regionMask_RS', 0, {fullfile(filePath2, 'Matfiles')});

rezRS_early = load(RSPath1{1}, 'rezRS'); 
rezRS_early = rezRS_early.('rezRS'); 

rezRS_late = load(RSPath2{1}, 'rezRS'); 
rezRS_late = rezRS_late.('rezRS'); 

cool = colormap('cool'); % for Go/No-Go
close; 

% rez1 = load(fullfile(filePath1, 'Matfiles', strcat(header1, '_dff_evtAligned_regionMask.mat')), 'rez');
% rez1 = rez1.('rez'); 
% [~, stimOn_timepts1] = temporalAlignInterp1(rez1.stimOnDffC.rs(:, 1), rez1.stimOnDffC.rs(:, 2), 0.001);
% 
% rez2 = load(fullfile(filePath2, 'Matfiles', strcat(header2, '_dff_evtAligned_regionMask.mat')), 'rez');
% rez2 = rez2.('rez'); 
% [~, stimOn_timepts2] = temporalAlignInterp1(rez2.stimOnDffC.rs(:, 1), rez2.stimOnDffC.rs(:, 2), 0.001);

%rezRS_learn.meanStimOnDffGo = [rezRS_early.meanStimOnDffGng(1, :);   rezRS_late.meanStimOnDffGng(1, :)]; 
%rezRS_learn.semStimOnDffGo = [rezRS_early.semStimOnDffGng(1, :);   rezRS_late.semStimOnDffGng(1, :)]; 

%% generate plots
% cue-aligned go vs nogo classification dffs
assert(isequal(rezRS_early.stimOnGoNogoSvmTs, rezRS_late.stimOnGoNogoSvmTs))
plotMeanWithoutSem([mean(rezRS_early.stimOnGoNogoSvm); mean(rezRS_late.stimOnGoNogoSvm)], ... 
    rezRS_early.stimOnGoNogoSvmTs, {header1, header2})
title("early vs late rs go-nogo SVM (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(figSaveDir, strcat(mouseId1, '_', date1, '_', date2, '_dff_cueAlign_gonogo_svm_early_late_rs.pdf')), '-dpdf', '-vector', '-bestfit')

% cue-aligned hit vs FA classification dffs
assert(isequal(rezRS_early.stimOnHitFaSvmTs, rezRS_late.stimOnHitFaSvmTs))
plotMeanWithoutSem([mean(rezRS_early.stimOnHitFaSvm); mean(rezRS_late.stimOnHitFaSvm)], ... 
    rezRS_early.stimOnHitFaSvmTs, {header1, header2})
title("early vs late rs hitVsFa SVM (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(figSaveDir, strcat(mouseId1, '_', date1, '_', date2, '_dff_cueAlign_hitFA_svm_early_late_rs.pdf')), '-dpdf', '-vector', '-bestfit')

% lick-aligned hit vs FA classification dffs
assert(isequal(rezRS_early.lickHitFaSvmTs, rezRS_late.lickHitFaSvmTs))
plotMeanWithoutSem([mean(rezRS_early.lickHitFaSvm); mean(rezRS_late.lickHitFaSvm)], ... 
    rezRS_early.lickHitFaSvmTs, {header1, header2})
title("early vs late rs hitVsFa SVM (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -1:0.5:1, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 1]); %ylim([-1 1])
print(fullfile(figSaveDir, strcat(mouseId1, '_', date1, '_', date2, '_dff_lickAlign_hitFA_svm_early_late_rs.pdf')), '-dpdf', '-vector', '-bestfit')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dffTsItp_mean, dffTsItp_sem, dffTsItp, dffTsItpTs] = getInterpolatedMeanSemDffFromDffTsCell(dffTsCell)
dffTsC = cellWithNonEmptyColumns(dffTsCell);
dffTsLinC = cellfun(@(a) linspace(-1, 1, length(a)), dffTsC(:, 2), 'UniformOutput', false);
[dffTsItp, dffTsItpTs] = temporalAlignInterp1(dffTsC(:, 1), dffTsLinC, 0.001);
[dffTsItp_mean, ~, dffTsItp_sem] = meanstdsem(cell2mat(dffTsItp));

end




