%% whereabouts
filePath = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA017/DA017_042224'; 
[~, header] = fileparts(filePath);
headerParts = regexp(header, '_', 'split');  % Split the string at the underscore
mouseId = headerParts{1};  % Part before the underscore
date = headerParts{2};  % Part after the underscore
figSaveDir = fullfile(fileparts(filePath), 'collectFigure'); 
if exist(figSaveDir, 'dir')~=7
    mkdir(figSaveDir)
end

%% load and prep data 
m1Path = GrabFiles_sort_trials('_regionMask_m1', 0, {fullfile(filePath, 'Matfiles')});
m2Path = GrabFiles_sort_trials('_regionMask_m2', 0, {fullfile(filePath, 'Matfiles')});
ssPath = GrabFiles_sort_trials('_regionMask_ss', 0, {fullfile(filePath, 'Matfiles')});
rsPath = GrabFiles_sort_trials('_regionMask_rs', 0, {fullfile(filePath, 'Matfiles')});
v1Path = GrabFiles_sort_trials('_regionMask_v1', 0, {fullfile(filePath, 'Matfiles')});

load(m1Path{1}, 'rezM1'); 
load(m2Path{1}, 'rezM2'); 
load(ssPath{1}, 'rezSS'); 
load(rsPath{1}, 'rezRS'); 
load(v1Path{1}, 'rezV1'); 

% prepare color
%cool = colormap('cool'); % for Go/No-Go
m1_color = [0 255 255]./255; 
m2_color = [64 191 255]./255; 
ss_color = [128 127 255]./255; 
rs_color = [191 64 255]./255; 
v1_color = [255 0 255]./255; 

colorC = {m1_color, m2_color, ss_color, rs_color, v1_color}; 

% get temporal info for PETHs
load(fullfile(filePath, 'Matfiles', strcat(header, '_dff_evtAligned_regionMask.mat')), 'rez');
[~, stimOn_timepts] = temporalAlignInterp1(rez.stimOnDffC.m1(:, 1), rez.stimOnDffC.m1(:, 2), 0.001);

xCorrLagMat = cell2mat(rez.xcorrDffLick.m1(:, 2)); 
xCorrLag = xCorrLagMat(1, :); 

%% generate plots
% cue-aligned HIT dffs
meanStimOnDffHitFA = [rezM1.meanStimOnDffHitFA(1, :); ...
                      rezM2.meanStimOnDffHitFA(1, :); ...
                      rezSS.meanStimOnDffHitFA(1, :); ...
                      rezRS.meanStimOnDffHitFA(1, :); ...
                      rezV1.meanStimOnDffHitFA(1, :)]; 

semStimOnDffHitFA = [rezM1.semStimOnDffHitFA(1, :); ...
                     rezM2.semStimOnDffHitFA(1, :); ...
                     rezSS.semStimOnDffHitFA(1, :); ...
                     rezRS.semStimOnDffHitFA(1, :); ...
                     rezV1.semStimOnDffHitFA(1, :)]; 

h_cueAlignedHit = plotMeanSemColorC(meanStimOnDffHitFA, semStimOnDffHitFA, stimOn_timepts, colorC, {'M1', 'M2', 'SS', 'RS', 'V1'}); 
title("cue aligned hit trials across regions (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(figSaveDir, strcat(mouseId, '_', date, '_dff_cueAlign_hit_acrossRegions.pdf')), '-dpdf', '-vector', '-bestfit')

% cross-correlograms
meanXcorrDffLick = [rezM1.meanXcorrDffLick; ...
                    rezM2.meanXcorrDffLick; ...
                    rezSS.meanXcorrDffLick; ...
                    rezRS.meanXcorrDffLick; ...
                    rezV1.meanXcorrDffLick]; 

semXcorrDffLick = [rezM1.semXcorrDffLick; ...
                   rezM2.semXcorrDffLick; ...
                   rezSS.semXcorrDffLick; ...
                   rezRS.semXcorrDffLick; ...
                   rezV1.semXcorrDffLick]; 

xcorrFig = plotMeanSemColorC(meanXcorrDffLick, semXcorrDffLick, xCorrLag, colorC, {'M1', 'M2', 'SS', 'RS', 'V1'}); 
xlabel('Time (s)'); ylabel('xcorr'); set(gca, 'TickDir', 'out', 'XTick', -2000:1000:2000); grid on
print(fullfile(figSaveDir, strcat(mouseId, '_', date, '_xcorr_dff_lick_acrossRegions.pdf')), '-dpdf', '-vector', '-bestfit')
close(xcorrFig)
fprintf("Processed the xcorr between lick and dff, and saved the crosscorrelogram!\n")

% cue-aligned go dffs
meanStimOnDffGo = [rezM1.meanStimOnDffGng(1, :); ...
                      rezM2.meanStimOnDffGng(1, :); ...
                      rezSS.meanStimOnDffGng(1, :); ...
                      rezRS.meanStimOnDffGng(1, :); ...
                      rezV1.meanStimOnDffGng(1, :)]; 

semStimOnDffGo = [rezM1.semStimOnDffGng(1, :); ...
                     rezM2.semStimOnDffGng(1, :); ...
                     rezSS.semStimOnDffGng(1, :); ...
                     rezRS.semStimOnDffGng(1, :); ...
                     rezV1.semStimOnDffGng(1, :)]; 

h_cueAlignedGo = plotMeanSemColorC(meanStimOnDffGo, semStimOnDffGo, stimOn_timepts, colorC, {'M1', 'M2', 'SS', 'RS', 'V1'}); 
title("cue aligned go trials across regions (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-0.3 0.69])
print(fullfile(figSaveDir, strcat(mouseId, '_', date, '_dff_cueAlign_go_acrossRegions.pdf')), '-dpdf', '-vector', '-bestfit')

% cue-aligned no-go dffs
meanStimOnDffNoGo = [rezM1.meanStimOnDffGng(2, :); ...
                      rezM2.meanStimOnDffGng(2, :); ...
                      rezSS.meanStimOnDffGng(2, :); ...
                      rezRS.meanStimOnDffGng(2, :); ...
                      rezV1.meanStimOnDffGng(2, :)]; 

semStimOnDffNoGo = [rezM1.semStimOnDffGng(2, :); ...
                     rezM2.semStimOnDffGng(2, :); ...
                     rezSS.semStimOnDffGng(2, :); ...
                     rezRS.semStimOnDffGng(2, :); ...
                     rezV1.semStimOnDffGng(2, :)]; 

h_cueAlignedNoGo = plotMeanSemColorC(meanStimOnDffNoGo, semStimOnDffNoGo, stimOn_timepts, colorC, {'M1', 'M2', 'SS', 'RS', 'V1'}); 
title("cue aligned nogo trials across regions (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); ylim([-0.3 0.69])
print(fullfile(figSaveDir, strcat(mouseId, '_', date, '_dff_cueAlign_nogo_acrossRegions.pdf')), '-dpdf', '-vector', '-bestfit')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dffTsItp_mean, dffTsItp_sem, dffTsItp, dffTsItpTs] = getInterpolatedMeanSemDffFromDffTsCell(dffTsCell)
dffTsC = cellWithNonEmptyColumns(dffTsCell);
dffTsLinC = cellfun(@(a) linspace(-1, 1, length(a)), dffTsC(:, 2), 'UniformOutput', false);
[dffTsItp, dffTsItpTs] = temporalAlignInterp1(dffTsC(:, 1), dffTsLinC, 0.001);
[dffTsItp_mean, ~, dffTsItp_sem] = meanstdsem(cell2mat(dffTsItp));

end




