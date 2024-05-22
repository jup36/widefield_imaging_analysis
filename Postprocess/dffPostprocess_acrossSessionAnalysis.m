%% whereabouts
filePath1 = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA019/DA019_041924'; 
[~, header1] = fileparts(filePath1);
header1Parts = regexp(header1, '_', 'split');  % Split the string at the underscore
mouseId1 = header1Parts{1};  % Part before the underscore
date1 = header1Parts{2};  % Part after the underscore

[~, header2] = fileparts(filePath2);
header2Parts = regexp(header2, '_', 'split');  % Split the string at the underscore
mouseId2 = header2Parts{1};  % Part before the underscore
date2 = header2Parts{2};  % Part after the underscore

filePath2 = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA019/DA019_042224'; 

figSaveDir = fullfile(fileparts(filePath1), 'collectFigure'); 
if exist(figSaveDir, 'dir')~=7
    mkdir(figSaveDir)
end

%% load and prep data 
m1Path1 = GrabFiles_sort_trials('_regionMask_m1', 0, {fullfile(filePath1, 'Matfiles')});
m1Path2 = GrabFiles_sort_trials('_regionMask_m1', 0, {fullfile(filePath2, 'Matfiles')});

rezM1_early = load(m1Path1{1}, 'rezM1'); 
rezM1_early = rezM1_early.('rezM1'); 

rezM1_late = load(m1Path2{1}, 'rezM1'); 
rezM1_late = rezM1_late.('rezM1'); 

cool = colormap('cool'); % for Go/No-Go

rez1 = load(fullfile(filePath1, 'Matfiles', strcat(header1, '_dff_evtAligned_regionMask.mat')), 'rez');
rez1 = rez1.('rez'); 
[~, stimOn_timepts1] = temporalAlignInterp1(rez1.stimOnDffC.m1(:, 1), rez1.stimOnDffC.m1(:, 2), 0.001);

rez2 = load(fullfile(filePath2, 'Matfiles', strcat(header2, '_dff_evtAligned_regionMask.mat')), 'rez');
rez2 = rez2.('rez'); 
[~, stimOn_timepts2] = temporalAlignInterp1(rez2.stimOnDffC.m1(:, 1), rez2.stimOnDffC.m1(:, 2), 0.001);

rezM1_learn.meanStimOnDffGo = [rezM1_early.meanStimOnDffGng(1, :);   rezM1_late.meanStimOnDffGng(1, :)]; 
rezM1_learn.semStimOnDffGo = [rezM1_early.semStimOnDffGng(1, :);   rezM1_late.semStimOnDffGng(1, :)]; 

%% generate plots
stimOnFig_learn_M1 = plotMeanSem(rezM1_learn.meanStimOnDffGo, rezM1_learn.semStimOnDffGo, stimOn_timepts1, {header1, header2});
title("Go Early vs Late sessions (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(figSaveDir, strcat(mouseId1, '_', date1, '_', date2, '_dff_cueOn_hit_early_late_m1')), '-dpdf', '-vector', '-bestfit')
%close(stimOnFig)

% lick rasters
[hitLickLatencyC_early, waterLatencyC_early] = hitLickRastersAuditoryGng(filePath1); 
[hitLickLatencyC_late, waterLatencyC_late] = hitLickRastersAuditoryGng(filePath2); 

hitLickLatencyFig = rasterPlotLickWaterCellGroups("Go Early vs Late sessions (cue onset at time=0)", {hitLickLatencyC_early, hitLickLatencyC_late}, ...
    {waterLatencyC_early, waterLatencyC_late}, cool, [-1 5], [2 4.5], 0.7, {header1, header2}); 
print(fullfile(figSaveDir, strcat(mouseId1, '_', date1, '_', date2, '_dff_cueOn_hit_early_late_m1_lickWaterRaster')), '-dpdf', '-vector', '-bestfit')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



