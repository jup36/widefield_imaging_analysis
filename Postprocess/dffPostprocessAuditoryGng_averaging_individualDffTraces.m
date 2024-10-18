function dffPostprocessAuditoryGng_averaging_individualDffTraces(filePath)
%filePath = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA008/DA008_101823';
% This function performs further analyses and plotting using data in
% 'rez'(outcome from 'dffPostprocessGng_averaging.m').
%

%% load rez tbytDat
[~, header] = fileparts(filePath);
load(fullfile(filePath, 'Matfiles', strcat(header, '_dff_evtAligned_regionMask.mat')), 'rez');
load(fullfile(filePath, 'Matfiles', strcat(header, '_tbytDat_parseGng')), 'tbytDat')

% user variables
m1_color = [0 240 240]./255;

if exist(fullfile(filePath, 'Figure'), 'dir')~=7
    mkdir(fullfile(filePath, 'Figure'))
end

figSavePath = fullfile(fileparts(filePath), 'collectFigure'); 

%% more indices for trial-type identification
waterI = cellfun(@(a) ~isempty(a), {tbytDat.water}); 
lickI = cellfun(@(a) ~isempty(a), {tbytDat.Lick}); 
airpuffI = cellfun(@(a) ~isempty(a), {tbytDat.airpuff}); 
goI = [tbytDat.rewardTrI]'==1; 
nogoI = [tbytDat.punishTrI]'==1; 
hitI = cellfun(@(a) ~isempty(a), {tbytDat.hitLicks})'; 
missI = cell2mat({tbytDat.rewardTrI})' & cellfun(@(a) isempty(a), {tbytDat.water})'; 
faI = cellfun(@(a) ~isempty(a), {tbytDat.faLicks})'; 
crI = cell2mat({tbytDat.punishTrI})' & cellfun(@(a) isempty(a), {tbytDat.airpuff})'; 
trialI = trialTypesFromTbytDatParseGng(tbytDat); 

%% example lick and dff traces
% exampleTrials = [2, 8, 15, 16, 21, 22, 111];
% x1 = 1;
% figure; hold on;
% for j = 1:length(exampleTrials)
% 
%     t = exampleTrials(j);
%     licks = x1+find(rez.lickOnTfBin{t, 1});
%     xend = x1+length(rez.dffsOnTfItpM1{t, 1})-1;
% 
%     plot(x1:x1+length(rez.dffsOnTfItpM1{t, 1})-1, rez.dffsOnTfItpM1{t, 1}, 'Color', m1_color, 'LineWidth', 1);
%     vertline(licks, [0 1], 'lineWidth', 1, 'alpha', 0.5);
%     vertline(xend-length(rez.dffsOnTfItpM1{t, 1})/2, ylim, 'lineWidth', 0.75, 'lineStyle', ':', 'alpha', 0.5);
% 
%     x1 = x1 + length(rez.dffsOnTfItpM1{t, 1}) + 1;
% end
% title("dff and lick traces");
% set(gca, 'TickDir', 'out');
% axis tight
% print(fullfile(filePath, 'Figure', 'dff_lickbouts_m1_example5trials.pdf'), '-dpdf', '-vector', '-bestfit')

%% Stim onset aligned activity (m1)
% align dffs take the mean and sem
% primary motor (stim onset)
% m1
[stimOnDff_m1_itp, stimOn_timepts] = temporalAlignInterp1(rez.stimOnDffC.m1(:, 1), rez.stimOnDffC.m1(:, 2), 0.001);
%figure; imagesc(smooth2a(cell2mat(stimOnDff_m1_itp(hitI)), 0, 100)); 
figure; imagesc(cell2mat(stimOnDff_m1_itp(hitI))); 
clim([-.5 .5])
set(gca, "TickDir", "out", "XTick", 0:1000:10000)
print(fullfile(figSavePath, strcat(header, "stimOnDff_m1_itp_hit")), '-dpdf', '-vector', '-bestfit')

% m1 hit lick align
hitFstLickDffm1 = cellWithNonEmptyColumns(rez.hitDffFirstLick.m1);
hitFstLickDffm1Ts = cellfun(@(a) linspace(-1, 1, length(a)), hitFstLickDffm1(:, 2), 'UniformOutput', false);
[hitFstLickDffm1Itp, hitFstLickDffm1ItpTs] = temporalAlignInterp1(hitFstLickDffm1(:, 1), hitFstLickDffm1Ts, 0.001);
figure; imagesc(cell2mat(hitFstLickDffm1Itp)); 
clim([-1 1])
set(gca, "TickDir", "out")
print(fullfile(figSavePath, strcat(header, "hitLickDff_m1_itp_hit")), '-dpdf', '-vector', '-bestfit')

% m2
[stimOnDff_m2_itp, stimOn_timepts] = temporalAlignInterp1(rez.stimOnDffC.m2(:, 1), rez.stimOnDffC.m2(:, 2), 0.001);
figure; imagesc(smooth2a(cell2mat(stimOnDff_m2_itp(hitI)), 0, 100)); 
clim([-.5 .5])
set(gca, "TickDir", "out", "XTick", 0:1000:10000)
print(fullfile(figSavePath, strcat(header, "stimOnDff_m2_itp_hit")), '-dpdf', '-vector', '-bestfit')

% m2 hit lick align
hitFstLickDffm2 = cellWithNonEmptyColumns(rez.hitDffFirstLick.m2);
hitFstLickDffm2Ts = cellfun(@(a) linspace(-1, 1, length(a)), hitFstLickDffm2(:, 2), 'UniformOutput', false);
[hitFstLickDffm2Itp, hitFstLickDffm2ItpTs] = temporalAlignInterp1(hitFstLickDffm2(:, 1), hitFstLickDffm2Ts, 0.001);
figure; imagesc(cell2mat(hitFstLickDffm2Itp)); 
clim([-1 1])
set(gca, "TickDir", "out")
print(fullfile(figSavePath, strcat(header, "hitLickDff_m2_itp_hit")), '-dpdf', '-vector', '-bestfit')

% ss
[stimOnDff_ss_itp, stimOn_timepts] = temporalAlignInterp1(rez.stimOnDffC.ss(:, 1), rez.stimOnDffC.ss(:, 2), 0.001);
figure; imagesc(smooth2a(cell2mat(stimOnDff_ss_itp(hitI)), 0, 100)); 
clim([-0.5 0.5])
set(gca, "TickDir", "out", "XTick", 0:1000:10000)
print(fullfile(figSavePath, strcat(header, "stimOnDff_ss_itp_hit")), '-dpdf', '-vector', '-bestfit')

% ss hit lick align
hitFstLickDffss = cellWithNonEmptyColumns(rez.hitDffFirstLick.ss);
hitFstLickDffssTs = cellfun(@(a) linspace(-1, 1, length(a)), hitFstLickDffss(:, 2), 'UniformOutput', false);
[hitFstLickDffssItp, hitFstLickDffssItpTs] = temporalAlignInterp1(hitFstLickDffss(:, 1), hitFstLickDffssTs, 0.001);
figure; imagesc(cell2mat(hitFstLickDffssItp)); 
clim([-1 1])
set(gca, "TickDir", "out")
print(fullfile(figSavePath, strcat(header, "hitLickDff_ss_itp_hit")), '-dpdf', '-vector', '-bestfit')

% rs
[stimOnDff_rs_itp, stimOn_timepts] = temporalAlignInterp1(rez.stimOnDffC.rs(:, 1), rez.stimOnDffC.rs(:, 2), 0.001);
figure; imagesc(smooth2a(cell2mat(stimOnDff_rs_itp(hitI)), 0, 100)); 
clim([-0.5 0.5])
set(gca, "TickDir", "out", "XTick", 0:1000:10000)
print(fullfile(figSavePath, strcat(header, "stimOnDff_rs_itp_hit")), '-dpdf', '-vector', '-bestfit')

% rs hit lick align
hitFstLickDffrs = cellWithNonEmptyColumns(rez.hitDffFirstLick.rs);
hitFstLickDffrsTs = cellfun(@(a) linspace(-1, 1, length(a)), hitFstLickDffrs(:, 2), 'UniformOutput', false);
[hitFstLickDffrsItp, hitFstLickDffrsItpTs] = temporalAlignInterp1(hitFstLickDffrs(:, 1), hitFstLickDffrsTs, 0.001);
figure; imagesc(cell2mat(hitFstLickDffrsItp)); 
clim([-1 1])
set(gca, "TickDir", "out")
print(fullfile(figSavePath, strcat(header, "hitLickDff_rs_itp_hit")), '-dpdf', '-vector', '-bestfit')

% v1
[stimOnDff_v1_itp, stimOn_timepts] = temporalAlignInterp1(rez.stimOnDffC.v1(:, 1), rez.stimOnDffC.v1(:, 2), 0.001);
figure; imagesc(smooth2a(cell2mat(stimOnDff_v1_itp(hitI)), 0, 100)); 
clim([-0.5 0.5])
set(gca, "TickDir", "out", "XTick", 0:1000:10000)
print(fullfile(figSavePath, strcat(header, "stimOnDff_v1_itp_hit")), '-dpdf', '-vector', '-bestfit')

% v1 hit lick align
hitFstLickDffv1 = cellWithNonEmptyColumns(rez.hitDffFirstLick.v1);
hitFstLickDffv1Ts = cellfun(@(a) linspace(-1, 1, length(a)), hitFstLickDffv1(:, 2), 'UniformOutput', false);
[hitFstLickDffv1Itp, hitFstLickDffv1ItpTs] = temporalAlignInterp1(hitFstLickDffv1(:, 1), hitFstLickDffv1Ts, 0.001);
figure; imagesc(cell2mat(hitFstLickDffv1Itp)); 
clim([-1 1])
set(gca, "TickDir", "out")
print(fullfile(figSavePath, strcat(header, "hitLickDff_v1_itp_hit")), '-dpdf', '-vector', '-bestfit')














