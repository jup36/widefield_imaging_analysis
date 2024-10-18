%% whereabouts
% filePaths = {'/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA017/DA017_041624', ...
%              '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA017/DA017_041824', ...
%              '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA017/DA017_042224', ...
%              '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA017/DA017_042924'}; 

filePaths = {'/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA019/DA019_041924', ...
             '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA019/DA019_2_042024', ...
             '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA019/DA019_042224', ...
             '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA019/DA019_042624'}; 

figSavePath = fullfile(fileparts(filePaths{1}), 'collectFigure'); 
blueShades = generateColormap([176 226 255]./255, [0 0 128]./255, 200); % for blocks

m1_color = [0 255 255]./255; 
m2_color = [64 191 255]./255;
ss_color = [128 127 255]./255;
rs_color = [191 64 255]./255; 
v1_color = [255 0 255]./255; 


for f = 1:length(filePaths)
    filePath = filePaths{f};
    filePath_region = GrabFiles_sort_trials('_regionMask.mat', 0, {fullfile(filePath, 'Matfiles')});
    filePath_tbytDat = GrabFiles_sort_trials('_tbytDat_parseGng.mat', 0, {fullfile(filePath, 'Matfiles')});

    [~, header] = fileparts(filePath);
    headerParts = regexp(header, '_', 'split');  % Split the string at the underscore
    mouseIdC{f,1} = headerParts{1};  % Part before the underscore
    dateC{f,1} = headerParts{2};  % Part after the underscore
    figSaveDir = fullfile(fileparts(filePath), 'collectFigure');
    if exist(figSaveDir, 'dir')~=7
        mkdir(figSaveDir)
    end
    
    %load
    load(filePath_region{1}, 'rez'); 
    load(filePath_tbytDat{1}, 'tbytDat')

    trI = trialTypesFromTbytDatParseGng(tbytDat);
    totTr = size(rez.xcorrDffLick.rs, 1); 

    % process m1
    xcorrDffLick_m1_hit = cell2mat(rez.xcorrDffLick.m1(trI.hitI(1:totTr), 1));
    % imagesc(xcorrDffLick_m1_hit)
    % set(gca, 'TickDir', 'out')
    % print(fullfile(figSavePath, 'tbyt_xcorrDffLick_m1_hit_DA019_042224.pdf'), '-dpdf', '-vector', '-bestfit')
    
    xcorrDffLick_m1_fa = cell2mat(rez.xcorrDffLick.m1(trI.faI(1:totTr), 1));
    denorm_m1_hit = nanmax(nanmax(xcorrDffLick_m1_hit(:, round(size(xcorrDffLick_m1_hit, 2)/2)-100:round(size(xcorrDffLick_m1_hit, 2)/2)+100))); 
    [norm_mean_xcorrDffLick_m1_hit_C{f, 1}, ~, norm_sem_xcorrDffLick_m1_hit_C{f, 1}] = meanstdsem(xcorrDffLick_m1_hit./denorm_m1_hit); 
    [mean_xcorrDffLick_m1_hit_C{f, 1}, ~, sem_xcorrDffLick_m1_hit_C{f, 1}] = meanstdsem(xcorrDffLick_m1_hit); 

    denorm_m1_fa = nanmax(nanmax(xcorrDffLick_m1_fa(:, round(size(xcorrDffLick_m1_fa, 2)/2)-100:round(size(xcorrDffLick_m1_fa, 2)/2)+100))); 
    [norm_mean_xcorrDffLick_m1_fa_C{f, 1}, ~, norm_sem_xcorrDffLick_m1_fa_C{f, 1}] = meanstdsem(xcorrDffLick_m1_fa./denorm_m1_fa); 
    [mean_xcorrDffLick_m1_fa_C{f, 1}, ~, sem_xcorrDffLick_m1_fa_C{f, 1}] = meanstdsem(xcorrDffLick_m1_fa); 

    % process m2
    xcorrDffLick_m2_hit = cell2mat(rez.xcorrDffLick.m2(trI.hitI(1:totTr), 1));
    % imagesc(xcorrDffLick_m2_hit)
    % set(gca, 'TickDir', 'out')
    % print(fullfile(figSavePath, 'tbyt_xcorrDffLick_m2_hit_DA019_042224.pdf'), '-dpdf', '-vector', '-bestfit')
    xcorrDffLick_m2_fa = cell2mat(rez.xcorrDffLick.m2(trI.faI(1:totTr), 1));
    denorm_m2_hit = nanmax(nanmax(xcorrDffLick_m2_hit(:, round(size(xcorrDffLick_m2_hit, 2)/2)-100:round(size(xcorrDffLick_m2_hit, 2)/2)+100))); 
    [norm_mean_xcorrDffLick_m2_hit_C{f, 1}, ~, norm_sem_xcorrDffLick_m2_hit_C{f, 1}] = meanstdsem(xcorrDffLick_m2_hit./denorm_m2_hit); 
    [mean_xcorrDffLick_m2_hit_C{f, 1}, ~, sem_xcorrDffLick_m2_hit_C{f, 1}] = meanstdsem(xcorrDffLick_m2_hit); 

    denorm_m2_fa = nanmax(nanmax(xcorrDffLick_m2_fa(:, round(size(xcorrDffLick_m2_fa, 2)/2)-100:round(size(xcorrDffLick_m2_fa, 2)/2)+100))); 
    [norm_mean_xcorrDffLick_m2_fa_C{f, 1}, ~, norm_sem_xcorrDffLick_m2_fa_C{f, 1}] = meanstdsem(xcorrDffLick_m2_fa./denorm_m2_fa); 
    [mean_xcorrDffLick_m2_fa_C{f, 1}, ~, sem_xcorrDffLick_m2_fa_C{f, 1}] = meanstdsem(xcorrDffLick_m2_fa); 

    % process ss
    xcorrDffLick_ss_hit = cell2mat(rez.xcorrDffLick.ss(trI.hitI(1:totTr), 1));
    xcorrDffLick_ss_fa = cell2mat(rez.xcorrDffLick.ss(trI.faI(1:totTr), 1));
    denorm_ss_hit = nanmax(nanmax(xcorrDffLick_ss_hit(:, round(size(xcorrDffLick_ss_hit, 2)/2)-100:round(size(xcorrDffLick_ss_hit, 2)/2)+100))); 
    [norm_mean_xcorrDffLick_ss_hit_C{f, 1}, ~, norm_sem_xcorrDffLick_ss_hit_C{f, 1}] = meanstdsem(xcorrDffLick_ss_hit./denorm_ss_hit); 
    [mean_xcorrDffLick_ss_hit_C{f, 1}, ~, sem_xcorrDffLick_ss_hit_C{f, 1}] = meanstdsem(xcorrDffLick_ss_hit); 

    denorm_ss_fa = nanmax(nanmax(xcorrDffLick_ss_fa(:, round(size(xcorrDffLick_ss_fa, 2)/2)-100:round(size(xcorrDffLick_ss_fa, 2)/2)+100))); 
    [norm_mean_xcorrDffLick_ss_fa_C{f, 1}, ~, norm_sem_xcorrDffLick_ss_fa_C{f, 1}] = meanstdsem(xcorrDffLick_ss_fa./denorm_ss_fa); 
    [mean_xcorrDffLick_ss_fa_C{f, 1}, ~, sem_xcorrDffLick_ss_fa_C{f, 1}] = meanstdsem(xcorrDffLick_ss_fa); 
    
    % process rs
    xcorrDffLick_rs_hit = cell2mat(rez.xcorrDffLick.rs(trI.hitI(1:totTr), 1));
    xcorrDffLick_rs_fa = cell2mat(rez.xcorrDffLick.rs(trI.faI(1:totTr), 1));
    denorm_rs_hit = nanmax(nanmax(xcorrDffLick_rs_hit(:, round(size(xcorrDffLick_rs_hit, 2)/2)-100:round(size(xcorrDffLick_rs_hit, 2)/2)+100))); 
    [norm_mean_xcorrDffLick_rs_hit_C{f, 1}, ~, norm_sem_xcorrDffLick_rs_hit_C{f, 1}] = meanstdsem(xcorrDffLick_rs_hit./denorm_rs_hit); 
    [mean_xcorrDffLick_rs_hit_C{f, 1}, ~, sem_xcorrDffLick_rs_hit_C{f, 1}] = meanstdsem(xcorrDffLick_rs_hit); 

    denorm_rs_fa = nanmax(nanmax(xcorrDffLick_rs_fa(:, round(size(xcorrDffLick_rs_fa, 2)/2)-100:round(size(xcorrDffLick_rs_fa, 2)/2)+100))); 
    [norm_mean_xcorrDffLick_rs_fa_C{f, 1}, ~, norm_sem_xcorrDffLick_rs_fa_C{f, 1}] = meanstdsem(xcorrDffLick_rs_fa./denorm_rs_fa); 
    [mean_xcorrDffLick_rs_fa_C{f, 1}, ~, sem_xcorrDffLick_rs_fa_C{f, 1}] = meanstdsem(xcorrDffLick_rs_fa); 

    % process v1
    xcorrDffLick_v1_hit = cell2mat(rez.xcorrDffLick.V1(trI.hitI(1:totTr), 1));
    xcorrDffLick_v1_fa = cell2mat(rez.xcorrDffLick.V1(trI.faI(1:totTr), 1));
    denorm_v1_hit = nanmax(nanmax(xcorrDffLick_v1_hit(:, round(size(xcorrDffLick_v1_hit, 2)/2)-100:round(size(xcorrDffLick_v1_hit, 2)/2)+100))); 
    [norm_mean_xcorrDffLick_v1_hit_C{f, 1}, ~, norm_sem_xcorrDffLick_v1_hit_C{f, 1}] = meanstdsem(xcorrDffLick_v1_hit./denorm_v1_hit); 
    [mean_xcorrDffLick_v1_hit_C{f, 1}, ~, sem_xcorrDffLick_v1_hit_C{f, 1}] = meanstdsem(xcorrDffLick_v1_hit); 

    denorm_v1_fa = nanmax(nanmax(xcorrDffLick_v1_fa(:, round(size(xcorrDffLick_v1_fa, 2)/2)-100:round(size(xcorrDffLick_v1_fa, 2)/2)+100))); 
    [norm_mean_xcorrDffLick_v1_fa_C{f, 1}, ~, norm_sem_xcorrDffLick_v1_fa_C{f, 1}] = meanstdsem(xcorrDffLick_v1_fa./denorm_v1_fa); 
    [mean_xcorrDffLick_v1_fa_C{f, 1}, ~, sem_xcorrDffLick_v1_fa_C{f, 1}] = meanstdsem(xcorrDffLick_v1_fa); 
    
    fprintf("processed file #%d\n", f)
    clearvars rez tbytDat
end

%% visualize rs xcorr
norm_mean_xcorrDffLick_rs = cell2mat(norm_mean_xcorrDffLick_rs_hit_C); 
norm_sem_xcorrDffLick_rs = cell2mat(norm_sem_xcorrDffLick_rs_hit_C); 

mean_xcorrDffLick_rs = cell2mat(mean_xcorrDffLick_rs_hit_C); 
sem_xcorrDffLick_rs = cell2mat(sem_xcorrDffLick_rs_hit_C); 

norm_xcorrFig_rs = plotMeanSemColormap(norm_mean_xcorrDffLick_rs, norm_sem_xcorrDffLick_rs, -2:0.001:2, blueShades);
xlabel('Time (s)'); ylabel('xcorr'); set(gca, 'TickDir', 'out', 'XTick', -2:1:2); grid on
print(fullfile(figSavePath, 'norm_xcorr_Gng_lick_dff_acrossSession_rs.pdf'), '-dpdf', '-vector', '-bestfit')
close(norm_xcorrFig_rs)

mean_xcorrDffLick_SI_rs = calculateXcorrSymmetryIndex(-2:0.001:2, mean_xcorrDffLick_rs); 
figure; hold on; 
plotScatterDatArrayColormap(mean_xcorrDffLick_SI_rs, blueShades);
hold off; 
%ylim([0.05 0.13])
print(fullfile(figSavePath, 'xcorr_Gng_lick_dff_xcorrSymInd_acrossSession_rs.pdf'), '-dpdf', '-vector', '-bestfit')

%% visualize m1 xcorr
norm_mean_xcorrDffLick_m1 = cell2mat(norm_mean_xcorrDffLick_m1_hit_C); 
norm_sem_xcorrDffLick_m1 = cell2mat(norm_sem_xcorrDffLick_m1_hit_C); 

mean_xcorrDffLick_m1 = cell2mat(mean_xcorrDffLick_m1_hit_C); 
sem_xcorrDffLick_m1 = cell2mat(sem_xcorrDffLick_m1_hit_C); 

norm_xcorrFig_m1 = plotMeanSemColormap(norm_mean_xcorrDffLick_m1, norm_sem_xcorrDffLick_m1, -2:0.001:2, blueShades);
xlabel('Time (s)'); ylabel('xcorr'); set(gca, 'TickDir', 'out', 'XTick', -2:1:2); grid on
print(fullfile(figSavePath, 'norm_xcorr_Gng_lick_dff_acrossSession_m1.pdf'), '-dpdf', '-vector', '-bestfit')
close(norm_xcorrFig_m1)

mean_xcorrDffLick_SI_m1 = calculateXcorrSymmetryIndex(-2:0.001:2, mean_xcorrDffLick_m1); 
figure; hold on; 
plotScatterDatArrayColormap(mean_xcorrDffLick_SI_m1, blueShades);
hold off; 
print(fullfile(figSavePath, 'xcorr_Gng_lick_dff_xcorrSymInd_acrossSession_m1.pdf'), '-dpdf', '-vector', '-bestfit')

%% visualize m2 xcorr
norm_mean_xcorrDffLick_m2 = cell2mat(norm_mean_xcorrDffLick_m2_hit_C); 
norm_sem_xcorrDffLick_m2 = cell2mat(norm_sem_xcorrDffLick_m2_hit_C); 

mean_xcorrDffLick_m2 = cell2mat(mean_xcorrDffLick_m2_hit_C); 
sem_xcorrDffLick_m2 = cell2mat(sem_xcorrDffLick_m2_hit_C); 

norm_xcorrFig_m2 = plotMeanSemColormap(norm_mean_xcorrDffLick_m2, norm_sem_xcorrDffLick_m2, -2:0.001:2, blueShades);
xlabel('Time (s)'); ylabel('xcorr'); set(gca, 'TickDir', 'out', 'XTick', -2:1:2); grid on
print(fullfile(figSavePath, 'norm_xcorr_Gng_lick_dff_acrossSession_m2.pdf'), '-dpdf', '-vector', '-bestfit')
close(norm_xcorrFig_m2)

mean_xcorrDffLick_SI_m2 = calculateXcorrSymmetryIndex(-2:0.001:2, mean_xcorrDffLick_m2); 
figure; hold on; 
plotScatterDatArrayColormap(mean_xcorrDffLick_SI_m2, blueShades);
hold off; 
print(fullfile(figSavePath, 'xcorr_Gng_lick_dff_xcorrSymInd_acrossSession_m2.pdf'), '-dpdf', '-vector', '-bestfit')

%% visualize ss xcorr
norm_mean_xcorrDffLick_ss = cell2mat(norm_mean_xcorrDffLick_ss_hit_C); 
norm_sem_xcorrDffLick_ss = cell2mat(norm_sem_xcorrDffLick_ss_hit_C); 

mean_xcorrDffLick_ss = cell2mat(mean_xcorrDffLick_ss_hit_C); 
sem_xcorrDffLick_ss = cell2mat(sem_xcorrDffLick_ss_hit_C); 

norm_xcorrFig_ss = plotMeanSemColormap(norm_mean_xcorrDffLick_ss, norm_sem_xcorrDffLick_ss, -2:0.001:2, blueShades);
xlabel('Time (s)'); ylabel('xcorr'); set(gca, 'TickDir', 'out', 'XTick', -2:1:2); grid on
print(fullfile(figSavePath, 'norm_xcorr_Gng_lick_dff_acrossSession_ss.pdf'), '-dpdf', '-vector', '-bestfit')
close(norm_xcorrFig_ss)

mean_xcorrDffLick_SI_ss = calculateXcorrSymmetryIndex(-2:0.001:2, mean_xcorrDffLick_ss); 
figure; hold on; 
plotScatterDatArrayColormap(mean_xcorrDffLick_SI_ss, blueShades);
hold off; 
print(fullfile(figSavePath, 'xcorr_Gng_lick_dff_xcorrSymInd_acrossSession_ss.pdf'), '-dpdf', '-vector', '-bestfit')

%% visualize v1 xcorr
norm_mean_xcorrDffLick_v1 = cell2mat(norm_mean_xcorrDffLick_v1_hit_C); 
norm_sem_xcorrDffLick_v1 = cell2mat(norm_sem_xcorrDffLick_v1_hit_C); 

mean_xcorrDffLick_v1 = cell2mat(mean_xcorrDffLick_v1_hit_C); 
sem_xcorrDffLick_v1 = cell2mat(sem_xcorrDffLick_v1_hit_C); 

norm_xcorrFig_v1 = plotMeanSemColormap(norm_mean_xcorrDffLick_v1, norm_sem_xcorrDffLick_v1, -2:0.001:2, blueShades);
xlabel('Time (s)'); ylabel('xcorr'); set(gca, 'TickDir', 'out', 'XTick', -2:1:2); grid on
print(fullfile(figSavePath, 'norm_xcorr_Gng_lick_dff_acrossSession_v1.pdf'), '-dpdf', '-vector', '-bestfit')
close(norm_xcorrFig_v1)

mean_xcorrDffLick_SI_v1 = calculateXcorrSymmetryIndex(-2:0.001:2, mean_xcorrDffLick_v1); 
figure; hold on; 
plotScatterDatArrayColormap(mean_xcorrDffLick_SI_v1, blueShades);
hold off; 
print(fullfile(figSavePath, 'xcorr_Gng_lick_dff_xcorrSymInd_acrossSession_v1.pdf'), '-dpdf', '-vector', '-bestfit')

%%
% mean_xcorrDffLick_SI_m1 = calculateXcorrSymmetryIndex(-2:0.001:2, mean_xcorrDffLick_m1); 
% mean_xcorrDffLick_SI_m2 = calculateXcorrSymmetryIndex(-2:0.001:2, mean_xcorrDffLick_m2); 
% mean_xcorrDffLick_SI_ss = calculateXcorrSymmetryIndex(-2:0.001:2, mean_xcorrDffLick_ss); 
% mean_xcorrDffLick_SI_rs = calculateXcorrSymmetryIndex(-2:0.001:2, mean_xcorrDffLick_rs); 
% mean_xcorrDffLick_SI_v1 = calculateXcorrSymmetryIndex(-2:0.001:2, mean_xcorrDffLick_v1); 

figure; hold on; 
plotScatterDatArrayColormap(mean_xcorrDffLick_SI_m1, repmat(m1_color, 4, 1));
plotScatterDatArrayColormap(mean_xcorrDffLick_SI_m2, repmat(m2_color, 4, 1));
plotScatterDatArrayColormap(mean_xcorrDffLick_SI_ss, repmat(ss_color, 4, 1));
plotScatterDatArrayColormap(mean_xcorrDffLick_SI_rs, repmat(rs_color, 4, 1));
plotScatterDatArrayColormap(mean_xcorrDffLick_SI_v1, repmat(v1_color, 4, 1));
hold off; 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dffTsItp_mean, dffTsItp_sem, dffTsItp, dffTsItpTs] = getInterpolatedMeanSemDffFromDffTsCell(dffTsCell)
dffTsC = cellWithNonEmptyColumns(dffTsCell);
dffTsLinC = cellfun(@(a) linspace(-1, 1, length(a)), dffTsC(:, 2), 'UniformOutput', false);
[dffTsItp, dffTsItpTs] = temporalAlignInterp1(dffTsC(:, 1), dffTsLinC, 0.001);
[dffTsItp_mean, ~, dffTsItp_sem] = meanstdsem(cell2mat(dffTsItp));

end




