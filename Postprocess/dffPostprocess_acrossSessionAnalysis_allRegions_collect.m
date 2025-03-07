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

regions = {'m1', 'm2', 'ss', 'rs', 'v1'}; 

hitLickLatencyC = cell(1, length(filePaths)); 
waterLatencyC = cell(1, length(filePaths)); 
headerC = cell(1,length(filePaths)); 

for ff = 1:length(filePaths)
    filePath = filePaths{ff};
    filePath_region = GrabFiles_sort_trials('_regionMask_M1', 0, {fullfile(filePath, 'Matfiles')});

    [~, header] = fileparts(filePath);
    headerC{1, ff} = header; 
    headerParts = regexp(header, '_', 'split');  % Split the string at the underscore
    mouseIdC{ff,1} = headerParts{1};  % Part before the underscore
    dateC{ff,1} = headerParts{2};  % Part after the underscore
    figSaveDir = fullfile(fileparts(filePath), 'collectFigure');
    if exist(figSaveDir, 'dir')~=7
        mkdir(figSaveDir)
    end
    
    %load
    load(filePath_region{1}, 'rezM1'); 
    
    % cue-aligned go
    mean_stimAlignGoC_M1{ff, 1} = rezM1.meanStimOnDffHitFA(1, :);
    sem_stimAlignGoC_M1{ff, 1}  = rezM1.semStimOnDffHitFA(1, :);  

    % hit-aligned go
    mean_hitFstLickDff_M1{ff, 1} = rezM1.meanHitFstLickDff(1, :); 
    sem_hitFstLickDff_M1{ff, 1} = rezM1.semHitFstLickDff(1, :);
    fprintf("processed file #%d\n", ff)
    
    % lick and water rasters 
    [hitLickLatencyC{1, ff}, waterLatencyC{1, ff}] = hitLickRastersAuditoryGng(filePath);
end

%% lick raster plot
pathToLickWaterRaster = GrabFiles_sort_trials('_AcrossSession_lickWaterRaster', 0, {figSaveDir});
if isempty(pathToLickWaterRaster)
    hitLickLatencyFig = rasterPlotLickWaterCellGroups("Lick rasters across session (cue onset at time=0)", hitLickLatencyC, ...
        waterLatencyC, cool, [-1 5], [2 4.5], 0.7, headerC);
    print(fullfile(figSaveDir, strcat(mouseIdC{1}, '_dff_cueOn_Hit_AcrossSession_lickWaterRaster')), '-dpdf', '-vector', '-bestfit')
end

%% plot cue-aligned go dffs
% color
learningColor = slanCM('blues', 5); 
learningColor = learningColor(2:5, :); 
learningColorC = mat2cell(learningColor, ones(size(learningColor, 1), 1), size(learningColor, 2));  
id_dateC = cellfun(@(a, b) strcat(a, '_', b), mouseIdC, dateC, 'UniformOutput', false); 
% min-max normalization
mean_stimAlignGoC_M1_norm = cell2mat(cellfun(@(a) a./(max(a)), mean_stimAlignGoC_M1(:, 1), 'UniformOutput', false)); 
% plot
%plotMeanSemColorC(cell2mat(mean_stimAlignGoC_M1(:, 1)), cell2mat(sem_stimAlignGoC_M1), -0.9:0.001:5, learningColorC, id_dateC)
plotMeanSemColorC(mean_stimAlignGoC_M1_norm, cell2mat(sem_stimAlignGoC_M1), -0.9:0.001:5, learningColorC, id_dateC)
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -5:0.2:5); xlim([-1 5]); %ylim([-0.4 1.1])
print(fullfile(figSaveDir, strcat('DA017', '_dff_cueOnAligned_hit_across_session_M1_fullX')), '-dpdf', '-vector', '-bestfit')
xlim([1 4])
print(fullfile(figSaveDir, strcat('DA017', '_dff_cueOnAligned_hit_across_session_M1')), '-dpdf', '-vector', '-bestfit')

%% plot hit-aligned dffs
% min-max normalization
mean_hitFstLickDff_M1_norm = cell2mat(cellfun(@(a) a./(max(a)), mean_hitFstLickDff_M1(:, 1), 'UniformOutput', false)); % plot
% plot
plotMeanSemColorC(mean_hitFstLickDff_M1_norm, cell2mat(sem_hitFstLickDff_M1), -1:0.001:1, learningColorC)
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -1:0.5:1, 'TickDir', 'out', 'YTick', -5:0.2:5); xlim([-1 1]); %ylim([-0.2 1.1])

print(fullfile(figSaveDir, strcat('DA017', '_dff_hitAligned_hit_across_session_M1')), '-dpdf', '-vector', '-bestfit')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dffTsItp_mean, dffTsItp_sem, dffTsItp, dffTsItpTs] = getInterpolatedMeanSemDffFromDffTsCell(dffTsCell)
dffTsC = cellWithNonEmptyColumns(dffTsCell);
dffTsLinC = cellfun(@(a) linspace(-1, 1, length(a)), dffTsC(:, 2), 'UniformOutput', false);
[dffTsItp, dffTsItpTs] = temporalAlignInterp1(dffTsC(:, 1), dffTsLinC, 0.001);
[dffTsItp_mean, ~, dffTsItp_sem] = meanstdsem(cell2mat(dffTsItp));

end




