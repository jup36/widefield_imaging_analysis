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
V1Path1 = GrabFiles_sort_trials('_regionMask_V1', 0, {fullfile(filePath1, 'Matfiles')});
V1Path2 = GrabFiles_sort_trials('_regionMask_V1', 0, {fullfile(filePath2, 'Matfiles')});

rezV1_early = load(V1Path1{1}, 'rezV1'); 
rezV1_early = rezV1_early.('rezV1'); 

rezV1_late = load(V1Path2{1}, 'rezV1'); 
rezV1_late = rezV1_late.('rezV1'); 

cool = colormap('cool'); % for Go/No-Go

rez1 = load(fullfile(filePath1, 'Matfiles', strcat(header1, '_dff_evtAligned_regionMask.mat')), 'rez');
rez1 = rez1.('rez'); 
[~, stimOn_timepts1] = temporalAlignInterp1(rez1.stimOnDffC.v1(:, 1), rez1.stimOnDffC.v1(:, 2), 0.001);

rez2 = load(fullfile(filePath2, 'Matfiles', strcat(header2, '_dff_evtAligned_regionMask.mat')), 'rez');
rez2 = rez2.('rez'); 
[~, stimOn_timepts2] = temporalAlignInterp1(rez2.stimOnDffC.v1(:, 1), rez2.stimOnDffC.v1(:, 2), 0.001);

rezV1_learn.meanStimOnDffGo = [rezV1_early.meanStimOnDffGng(1, :);   rezV1_late.meanStimOnDffGng(1, :)]; 
rezV1_learn.semStimOnDffGo = [rezV1_early.semStimOnDffGng(1, :);   rezV1_late.semStimOnDffGng(1, :)]; 

%% generate plots
% cue-aligned hit dffs
stimOnFig_learn_V1 = plotMeanSem(rezV1_learn.meanStimOnDffGo, rezV1_learn.semStimOnDffGo, stimOn_timepts1, {header1, header2});
title("Go Early vs Late sessions (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(figSaveDir, strcat(mouseId1, '_', date1, '_', date2, '_dff_cueOn_hit_early_late_V1')), '-dpdf', '-vector', '-bestfit')
%close(stimOnFig)

% 1st hit lick-aligned dffs
[mHitDffFirstLick_early, sHitDffFirstLick_early, ~, hitDffFirstLick_early_ts] = getInterpolatedMeanSemDffFromDffTsCell(rez1.hitDffFirstLick.v1, [-1 1]);
[mHitDffFirstLick_late, sHitDffFirstLick_late, ~, hitDffFirstLick_late_ts] = getInterpolatedMeanSemDffFromDffTsCell(rez2.hitDffFirstLick.v1, [-1 1]);

hitLickEarlyLateFig = plotMeanSem([mHitDffFirstLick_early; mHitDffFirstLick_late], ...
    [sHitDffFirstLick_early; sHitDffFirstLick_late], ...
    hitDffFirstLick_early_ts, {header1, header2}); 
title("Lick Hit early vs late")
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -1:.5:1, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); grid on
print(fullfile(figSaveDir, strcat(mouseId1, '_', date1, '_', date2, '_dff_lick_Hit_early_late_v1.pdf')), '-dpdf', '-vector', '-bestfit')
close(hitLickEarlyLateFig)

% lick rasters
pathToLickWaterRaster = GrabFiles_sort_trials('_lickWaterRaster', 0, {figSaveDir});
if isempty(pathToLickWaterRaster)
    [hitLickLatencyC_early, waterLatencyC_early] = hitLickRastersAuditoryGng(filePath1);
    [hitLickLatencyC_late, waterLatencyC_late] = hitLickRastersAuditoryGng(filePath2);

    hitLickLatencyFig = rasterPlotLickWaterCellGroups("Go Early vs Late sessions (cue onset at time=0)", {hitLickLatencyC_early, hitLickLatencyC_late}, ...
        {waterLatencyC_early, waterLatencyC_late}, cool, [-1 5], [2 4.5], 0.7, {header1, header2});
    print(fullfile(figSaveDir, strcat(mouseId1, '_', date1, '_', date2, '_dff_cueOn_hit_early_late_V1_lickWaterRaster')), '-dpdf', '-vector', '-bestfit')
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dffTsItp_mean, dffTsItp_sem, dffTsItp, dffTsItpTs] = getInterpolatedMeanSemDffFromDffTsCell(dffTsCell)
dffTsC = cellWithNonEmptyColumns(dffTsCell);
dffTsLinC = cellfun(@(a) linspace(-1, 1, length(a)), dffTsC(:, 2), 'UniformOutput', false);
[dffTsItp, dffTsItpTs] = temporalAlignInterp1(dffTsC(:, 1), dffTsLinC, 0.001);
[dffTsItp_mean, ~, dffTsItp_sem] = meanstdsem(cell2mat(dffTsItp));

end




