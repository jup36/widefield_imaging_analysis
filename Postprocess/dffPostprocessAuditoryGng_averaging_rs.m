function dffPostprocessAuditoryGng_averaging_rs(filePath)
%filePath = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA008/DA008_101823';
% This function performs further analyses and plotting using data in
% 'rez'(outcome from 'dffPostprocessGng_averaging.m').
%

%% load rez tbytDat
[~, header] = fileparts(filePath);
load(fullfile(filePath, 'Matfiles', strcat(header, '_dff_evtAligned_regionMask.mat')), 'rez');
load(fullfile(filePath, 'Matfiles', strcat(header, '_tbytDat_parseGng')), 'tbytDat')

% user variables
rs_color = [64 191 255]./255; 

if exist(fullfile(filePath, 'Figure'), 'dir')~=7
    mkdir(fullfile(filePath, 'Figure'))
end

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

%% cross correlogram lick bouts and dff
xCorrLagMat = cell2mat(rez.xcorrDffLick.rs(:, 2)); 
xCorrLag = xCorrLagMat(1, :); 

[rezRS.meanXcorrDffLick, ~, rezRS.semXcorrDffLick] = meanstdsem(cell2mat(rez.xcorrDffLick.rs(:, 1)));
xcorrFig = plotMeanSemColor(rezRS.meanXcorrDffLick, rezRS.semXcorrDffLick, xCorrLag, rs_color, {'xcorr rs'});
xlabel('Time (s)'); ylabel('xcorr'); set(gca, 'TickDir', 'out', 'XTick', -2000:1000:2000); grid on
print(fullfile(filePath, 'Figure', 'xcorr_Gng_lick_dff_rs.pdf'), '-dpdf', '-vector', '-bestfit')
close(xcorrFig)
fprintf("Processed the xcorr between lick and dff, and saved the crosscorrelogram!\n")

%% example lick and dff traces
% exampleTrials = [2, 8, 15, 16, 21, 22, 111];
% x1 = 1;
% figure; hold on;
% for j = 1:length(exampleTrials)
% 
%     t = exampleTrials(j);
%     licks = x1+find(rez.lickOnTfBin{t, 1});
%     xend = x1+length(rez.dffsOnTfItpRS{t, 1})-1;
% 
%     plot(x1:x1+length(rez.dffsOnTfItpRS{t, 1})-1, rez.dffsOnTfItpRS{t, 1}, 'Color', rs_color, 'LineWidth', 1);
%     vertline(licks, [0 1], 'lineWidth', 1, 'alpha', 0.5);
%     vertline(xend-length(rez.dffsOnTfItpRS{t, 1})/2, ylim, 'lineWidth', 0.75, 'lineStyle', ':', 'alpha', 0.5);
% 
%     x1 = x1 + length(rez.dffsOnTfItpRS{t, 1}) + 1;
% end
% title("dff and lick traces");
% set(gca, 'TickDir', 'out');
% axis tight
% print(fullfile(filePath, 'Figure', 'dff_lickbouts_rs_example5trials.pdf'), '-dpdf', '-vector', '-bestfit')

%% Stim onset aligned activity (rs)
% align dffs take the mean and sem
% primary motor (stim onset)
[stimOnDff_rs_itp, stimOn_timepts] = temporalAlignInterp1(rez.stimOnDffC.rs(:, 1), rez.stimOnDffC.rs(:, 2), 0.001);
[rezRS.meanStimOnDff, ~, rezRS.semStimOnDff] = meanstdsem(cell2mat(stimOnDff_rs_itp));
[stimOnDff_rs_go_itp, ~] = temporalAlignInterp1(rez.stimOnDffC.rs(goI, 1), rez.stimOnDffC.rs(goI, 2), 0.001);
[stimOnDff_rs_nogo_itp, ~] = temporalAlignInterp1(rez.stimOnDffC.rs(nogoI, 1), rez.stimOnDffC.rs(nogoI, 2), 0.001);

[stimOnDff_rs_hit_itp, ~] = temporalAlignInterp1(rez.stimOnDffC.rs(hitI, 1), rez.stimOnDffC.rs(hitI, 2), 0.001);
[stimOnDff_rs_fa_itp, ~] = temporalAlignInterp1(rez.stimOnDffC.rs(faI, 1), rez.stimOnDffC.rs(faI, 2), 0.001);
[stimOnDff_rs_miss_itp, ~] = temporalAlignInterp1(rez.stimOnDffC.rs(missI, 1), rez.stimOnDffC.rs(missI, 2), 0.001);
[stimOnDff_rs_cr_itp, ~] = temporalAlignInterp1(rez.stimOnDffC.rs(crI, 1), rez.stimOnDffC.rs(crI, 2), 0.001);

% primary motor (stim onset) go no-go
[rezRS.meanStimOnDffGng, rezRS.semStimOnDffGng] = trialGroupMeanSem(stimOnDff_rs_itp, {[tbytDat.rewardTrI], [tbytDat.punishTrI]});
stimOnFig = plotMeanSem(rezRS.meanStimOnDffGng, rezRS.semStimOnDffGng, stimOn_timepts, {'Go', 'NoGo'});
title("Go Nogo rs (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_cueOn_Gng_rs.pdf'), '-dpdf', '-vector', '-bestfit')
close(stimOnFig)
fprintf("Saved tone(stim)-aligned dff plot!\n")

% primary motor (stim onset) hit vs cr
[rezRS.meanStimOnDffHitCr, rezRS.semStimOnDffHitCr] = trialGroupMeanSem(stimOnDff_rs_itp, {hitI, crI});
stimOnHitCrFig = plotMeanSem(rezRS.meanStimOnDffHitCr, rezRS.semStimOnDffHitCr, stimOn_timepts, {'Hit', 'CR'});
title("Hit vs Cr rs (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_cueOn_HitCr_rs.pdf'), '-dpdf', '-vector', '-bestfit')
close(stimOnHitCrFig)
fprintf("Saved tone(stim)-aligned dff plot comparing Hit vs Cr!\n")

% primary motor (stim onset) hit vs FA
[rezRS.meanStimOnDffHitFA, rezRS.semStimOnDffHitFA] = trialGroupMeanSem(stimOnDff_rs_itp, {hitI, faI});
stimOnHitFaFig = plotMeanSem(rezRS.meanStimOnDffHitFA, rezRS.semStimOnDffHitFA, stimOn_timepts, {'Hit', 'FA'});
title("Hit vs FA rs (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_cueOn_HitFA_rs.pdf'), '-dpdf', '-vector', '-bestfit')
close(stimOnHitFaFig)
fprintf("Saved tone(stim)-aligned dff plot comparing Hit vs FA!\n")

% primary motor (stim onset) Miss vs CR
[rezRS.meanStimOnDffMissCR, rezRS.semStimOnDffMissCR] = trialGroupMeanSem(stimOnDff_rs_itp, {missI, crI});
stimOnMissCrFig = plotMeanSem(rezRS.meanStimOnDffMissCR, rezRS.semStimOnDffMissCR, stimOn_timepts, {'Miss', 'CR'});
title("Miss vs CR rs (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_cueOn_MissCR_rs.pdf'), '-dpdf', '-vector', '-bestfit')
close(stimOnMissCrFig)
fprintf("Saved tone(stim)-aligned dff plot comparing Miss vs CR!\n")

% primary motor (stim onset) FA vs CR
[rezRS.meanStimOnDffFaCR, rezRS.semStimOnDffFaCR] = trialGroupMeanSem(stimOnDff_rs_itp, {faI, crI});
stimOnFaCrFig = plotMeanSem(rezRS.meanStimOnDffFaCR, rezRS.semStimOnDffFaCR, stimOn_timepts, {'FA', 'CR'});
title("FA vs CR rs (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_cueOn_FaCR_rs.pdf'), '-dpdf', '-vector', '-bestfit')
close(stimOnFaCrFig)
fprintf("Saved tone(stim)-aligned dff plot comparing Miss vs CR!\n")

% rs (stim onset) go vs no-go (train classifier)
[rezRS.stimOnGoNogoSvm, rezRS.stimOnGoNogoSvmTs, rezRS.stimOnGoNogoNb, rezRS.stimOnGoNogoNbTs] = trainDffClassifier(stimOnDff_rs_go_itp, stimOnDff_rs_nogo_itp, stimOn_timepts, 50, 50, 10); 

% rs (stim onset) hit vs CR (train classifier)
[rezRS.stimOnHitCrSvm, rezRS.stimOnHitCrSvmTs, rezRS.stimOnHitCrNb, rezRS.stimOnHitCrNbTs] = trainDffClassifier(stimOnDff_rs_hit_itp, stimOnDff_rs_cr_itp, stimOn_timepts, 50, 50, 10); 

% rs (stim onset) hit vs FA (train classifier)
[rezRS.stimOnHitFaSvm, rezRS.stimOnHitFaSvmTs, rezRS.stimOnHitFaNb, rezRS.stimOnHitFaNbTs] = trainDffClassifier(stimOnDff_rs_hit_itp, stimOnDff_rs_fa_itp, stimOn_timepts, 50, 50, 10); 

% rs (stim onset) hit vs Miss (train classifier)
[rezRS.stimOnHitMissSvm, rezRS.stimOnHitMissSvmTs, rezRS.stimOnHitMissNb, rezRS.stimOnHitMissNbTs] = trainDffClassifier(stimOnDff_rs_hit_itp, stimOnDff_rs_miss_itp, stimOn_timepts, 50, 50, 10); 

% rs (stim onset) hit vs Miss (train classifier)
[rezRS.stimOnMissCrSvm, rezRS.stimOnMissCrSvmTs, rezRS.stimOnMissCrNb, rezRS.stimOnMissCrNbTs] = trainDffClassifier(stimOnDff_rs_miss_itp, stimOnDff_rs_cr_itp, stimOn_timepts, 50, 50, 10); 

fprintf("Completed cueOn-aligned dff and classification of trial types!\n")

%% iti licks (rs) Technically there's no ITI licks with the retractable spout
% itiLickDffrs = flatCell(rez.itiLickDff.rs); % unnest the cell arrays
% itiLickDffrsTs = cellfun(@(a) linspace(-1, 1, length(a)), itiLickDffrs(:, 2), 'UniformOutput', false);
% [itiLickDffrsItp, itiLickDffrsItpTs] = temporalAlignInterp1(itiLickDffrs (:, 1), itiLickDffrsTs);
% [rez.meanItiLickDff.rs, ~, rez.semItiLickDff.rs] = meanstdsem(cell2mat(itiLickDffrsItp));
% 
% plotMeanSemColor(rez.meanItiLickDff.rs, rez.semItiLickDff.rs, itiLickDffrsItpTs, rs_color, {'rs'});
% title("lick during ITI");
% xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); ylim([-.3 .3])

%% postStim licks (rs)
postStimLickDffrs = flatCell(rez.postStimLickDff.rs); % unnest the cell arrays
postStimLickDffrsTs = cellfun(@(a) linspace(-1, 1, length(a)), postStimLickDffrs(:, 2), 'UniformOutput', false);
[postStimLickDffrsItp, postStimLickDffrsItpTs] = temporalAlignInterp1(postStimLickDffrs(:, 1), postStimLickDffrsTs, 0.001);
[rezRS.meanPostStimLickDff, ~, rezRS.semPostStimLickDff] = meanstdsem(cell2mat(postStimLickDffrsItp));

postStimLickFig = plotMeanSemColor(rezRS.meanPostStimLickDff, rezRS.semPostStimLickDff, postStimLickDffrsItpTs, rs_color, {'rs'});
title("postStim licks");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); %ylim([-.3 .3])
print(fullfile(filePath, 'Figure', 'dff_postStimLick_Gng_rs.pdf'), '-dpdf', '-vector', '-bestfit')
close(postStimLickFig)
fprintf("Saved postStimLick-aligned dff plot!\n")

%% hit first licks (rs)
hitFstLickDffrs = cellWithNonEmptyColumns(rez.hitDffFirstLick.rs);
hitFstLickDffrsTs = cellfun(@(a) linspace(-1, 1, length(a)), hitFstLickDffrs(:, 2), 'UniformOutput', false);
[hitFstLickDffrsItp, hitFstLickDffrsItpTs] = temporalAlignInterp1(hitFstLickDffrs(:, 1), hitFstLickDffrsTs, 0.001);
[rezRS.meanHitFstLickDff, ~, rezRS.semHitFstLickDff] = meanstdsem(cell2mat(hitFstLickDffrsItp));

hitFstLickFig = plotMeanSemColor(rezRS.meanHitFstLickDff, rezRS.semHitFstLickDff, hitFstLickDffrsItpTs, rs_color, {'rs'});
title("Lick (first) Hit trial");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_hitFirstLick_Gng_rs.pdf'), '-dpdf', '-vector', '-bestfit')
close(hitFstLickFig)
fprintf("Saved hitFirstLick-aligned dff plot!\n")

%% first water (rs)
fstWaterDffrs = cellWithNonEmptyColumns(rez.firstWater.rs);
fstWaterDffrsTs = cellfun(@(a) linspace(-1, 1, length(a)), fstWaterDffrs(:, 2), 'UniformOutput', false);
[fstWaterDffrsItp, fstWaterDffrsItpTs] = temporalAlignInterp1(fstWaterDffrs (:, 1), fstWaterDffrsTs, 0.001);
[rezRS.meanFstWaterDff, ~, rezRS.semFstWaterDff] = meanstdsem(cell2mat(fstWaterDffrsItp));

firstWaterFig = plotMeanSemColor(rezRS.meanFstWaterDff, rezRS.semFstWaterDff, fstWaterDffrsItpTs, rs_color, {'rs'});
title("Water (first)");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); grid on;
print(fullfile(filePath, 'Figure', 'dff_firstWater_rs.pdf'), '-dpdf', '-vector', '-bestfit')
close(firstWaterFig)
fprintf("Saved firstWater-aligned dff plot!\n")

%% first air (rs)
if isfield(rez, 'firstAirpuff')
    fstAirDffrs = cellWithNonEmptyColumns(rez.firstAirpuff.rs);
    fstAirDffrsTs = cellfun(@(a) linspace(-1, 1, length(a)), fstAirDffrs(:, 2), 'UniformOutput', false);
    [fstAirDffrsItp, fstAirDffrsItpTs] = temporalAlignInterp1(fstAirDffrs (:, 1), fstAirDffrsTs, 0.001);
    [rezRS.meanFstAirDff, ~, rezRS.semFstAirDff] = meanstdsem(cell2mat(fstAirDffrsItp));

    fstAirFig = plotMeanSemColor(rezRS.meanFstAirDff, rezRS.semFstAirDff, fstAirDffrsItpTs, rs_color, {'rs'});
    title("Air (first)");
    xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); %ylim([-0.8 0.8])

    print(fullfile(filePath, 'Figure', 'dff_firstAirpuff_rs.pdf'), '-dpdf', '-vector', '-bestfit')
    close(fstAirFig)
    fprintf("Saved firstAir-aligned dff plot!\n")
end

%% FA first licks (rs)
faFstLickDffrs = cellWithNonEmptyColumns(rez.faDffFirstLick.rs);
faFstLickDffrsTs = cellfun(@(a) linspace(-1, 1, length(a)), faFstLickDffrs(:, 2), 'UniformOutput', false);
[faFstLickDffrsItp, faFstLickDffrsItpTs] = temporalAlignInterp1(faFstLickDffrs (:, 1), faFstLickDffrsTs, 0.001);
[rezRS.meanFaFstLickDff, ~, rezRS.semFaFstLickDff] = meanstdsem(cell2mat(faFstLickDffrsItp));

faFstLickFig = plotMeanSem(rezRS.meanFaFstLickDff, rezRS.semFaFstLickDff, faFstLickDffrsItpTs, {'rs'});
title("Lick (first) Fa trial");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); 
print(fullfile(filePath, 'Figure', 'dff_faFstLickFig_rs.pdf'), '-dpdf', '-vector', '-bestfit')
close(faFstLickFig)

%% Compare licks during Hit vs FA vs postStim (rs)
% Without normalization
hitFaPostStimFig = plotMeanSem([rezRS.meanHitFstLickDff; rezRS.meanFaFstLickDff; rezRS.meanPostStimLickDff], ...
    [rezRS.semHitFstLickDff; rezRS.semFaFstLickDff; rezRS.semPostStimLickDff], ...
    hitFstLickDffrsItpTs, {'Hit', 'FA', 'postStim'}); 
title("Lick Hit vs FA vs postStim")
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); grid on
print(fullfile(filePath, 'Figure', 'dff_lick_Hit_Fa_postStim_rs.pdf'), '-dpdf', '-vector', '-bestfit')
close(hitFaPostStimFig)

% With normalization
rezRS.meanHitFstLickDffBaseSub = baseSubNorm(rezRS.meanHitFstLickDff, hitFstLickDffrsItpTs, [-1 -0.5]);
rezRS.meanFaFstLickDffBaseSub = baseSubNorm(rezRS.meanFaFstLickDff, faFstLickDffrsItpTs, [-1 -0.5]);
rezRS.meanPostStimLickDffBaseSub = baseSubNorm(rezRS.meanPostStimLickDff, postStimLickDffrsItpTs, [-1 -0.5]);

hitFaPostStimFig_baseSub = plotMeanSem([rezRS.meanHitFstLickDffBaseSub; rezRS.meanFaFstLickDffBaseSub; rezRS.meanPostStimLickDffBaseSub], ...
    [rezRS.semHitFstLickDff; rezRS.semFaFstLickDff; rezRS.semPostStimLickDff], ...
    hitFstLickDffrsItpTs, {'Hit', 'FA', 'postStim'}); 
title("Lick Hit vs FA vs ITI vs postStim (baseSub)")
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); grid on
print(fullfile(filePath, 'Figure', 'dff_lick_Hit_Fa_postStim_baseSub_rs.pdf'), '-dpdf', '-vector', '-bestfit')
close(hitFaPostStimFig_baseSub)

%% train svm to classify licks at different task epochs
% rs (stim onset) go vs no-go (train classifier)
[rezRS.lickHitFaSvm, rezRS.lickHitFaSvmTs, rezRS.lickHitFaNb, rezRS.lickHitFaNbTs, X_hit, X_fa] = trainDffClassifier(hitFstLickDffrsItp, faFstLickDffrsItp, hitFstLickDffrsItpTs, 50, 50, 10); 

% imagesc the hit and fa trials dff
X_hit_fig = imagescWithTimeInfo(smooth2a(X_hit, 0, 1), rezRS.lickHitFaSvmTs); 
set(gca, 'XTick', -2:.5:4, 'TickDir', 'out');
clim([-0.7 0.7])
print(fullfile(filePath, 'Figure', 'dff_hit_timeBin_by_trial_rs.pdf'), '-dpdf', '-vector', '-bestfit')
close(X_hit_fig)

X_fa_fig = imagescWithTimeInfo(smooth2a(X_fa, 0, 1), rezRS.lickHitFaSvmTs); 
set(gca, 'XTick', -2:.5:4, 'TickDir', 'out');
clim([-0.7 0.7])
print(fullfile(filePath, 'Figure', 'dff_fa_timeBin_by_trial_rs.pdf'), '-dpdf', '-vector', '-bestfit')
close(X_fa_fig)

%% Miss CueOn (rs)
% primary motor (stim onset) go no-go
missCueI = [tbytDat.rewardTrI] & ~waterI & ~lickI;  
[rezRS.meanStimOnMissDff, rezRS.semStimOnMissDff] = trialGroupMeanSem(stimOnDff_rs_itp, {missCueI});
stimOnMissFig = plotMeanSem(rezRS.meanStimOnMissDff, rezRS.semStimOnMissDff, stimOn_timepts, {'rs miss'});
title("Go Miss rs (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_cueOn_goMiss_rs.pdf'), '-dpdf', '-vector', '-bestfit')
close(stimOnMissFig)
fprintf("Saved tone(stim)-aligned dff of Miss trials!\n")

%% CR CueOn (rs)
% primary motor (stim onset) go no-go
crCueI = [tbytDat.punishTrI] & ~airpuffI & ~lickI;  
[rezRS.meanStimOnCrDff, rezRS.semStimOnCrDff] = trialGroupMeanSem(stimOnDff_rs_itp, {crCueI});
stimOnCrFig = plotMeanSem(rezRS.meanStimOnCrDff, rezRS.semStimOnCrDff, stimOn_timepts, {'rs correct rejection'});
title("NoGo CR rs (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_cueOn_nogoCR_rs.pdf'), '-dpdf', '-vector', '-bestfit')
close(stimOnCrFig)
fprintf("Saved tone(stim)-aligned dff of Correct Rejection trials!\n")

% %% plot stim onset PETHs together
% plotMeanSem([rez.mStimOnDffrs; rez.mStimOnDffRS; rez.mStimOnDffS1; rez.mStimOnDffRs; rez.mStimOnDffV], ...
%     [rez.semStimOnDffrs; rez.semStimOnDffRS; rez.semStimOnDffS1; rez.semStimOnDffRs; rez.semStimOnDffV], ...
%     stimOn_timepts, ...
%     {'rs', 'RS', 'S1', 'RS', 'V'});
% xlabel('Time (s)')
% ylabel('DFF')
% set(gca, 'XTick', -2:1:4)
% xlim([-2 2])
% print(fullfile(filePath, 'Figure', 'dff_cueOn_Collect.pdf'), '-dpdf', '-vector', '-bestfit')

%% plot classification accuracy
% plotMeanSem(smooth2a([rez.logitGngCueOnrs; rez.logitGngCueOnRS; rez.logitGngCueOnS1; rez.logitGngCueOnRs; rez.logitGngCueOnV], 0, 2), ...
%     zeros(5, length(stimOn_timepts)), stimOn_timepts, {'rs', 'RS', 'S1', 'RS', 'V'});
% title("Cross-validated classification accuracy")
% xlabel('Time (s)')
% ylabel('Classification Accuracy')
% set(gca, 'XTick', -2:1:4, 'YTick', -1:0.1:1, 'TickDir', 'out')
% xlim([-2 2])
% ylim([0.45 0.85])
% print(fullfile(filePath, 'Figure', 'logisticRegression_Gng_cueOn_Collect.pdf'), '-dpdf', '-vector', '-bestfit')

% baseline dff image
% rez.meanBaseDffImgMotor = nanmean(cell2mat(reshape(rwdDffC_motor(:, 3), [1, 1, length(tbytDat)])), 3);
% imageFrameWithNaNsEdgeMap(rez.meanBaseDffImgMotor, [-.3 .3], dorsalMaps.edgeMap, [0, 0, 1], 0.3)
% title('mean baseline rs')
% set(gca, 'XTickLabel', {})
% set(gca, 'YTickLabel', {})
% print(fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA003/DA003_083023/Figure', 'dff_base_rs_Img.pdf'), '-dpdf', '-vector', '-bestfit')
% 
% % reward dff image
% rez.meanRwdDffImgMotor = nanmean(cell2mat(reshape(rwdDffC_motor(:, 4), [1, 1, length(tbytDat)])), 3);
% imageFrameWithNaNsEdgeMap(rez.meanRwdDffImgMotor, [-.3 .3], dorsalMaps.edgeMap, [0, 0, 1], 0.3)
% print(fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA003/DA003_083023/Figure', 'dff_reward_rs_Img.pdf'), '-dpdf', '-vector', '-bestfit')
% title('mean reward rs')
% set(gca, 'XTickLabel', {})
% set(gca, 'YTickLabel', {})
% print(fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA003/DA003_083023/Figure', 'dff_reward_rs_Img.pdf'), '-dpdf', '-vector', '-bestfit')

%showAllenEdgesOnTransformedWF(rez.meanBaseDffImgMotor, dorsalMaps, 0, 1)

%imageFrameWithNaNs(rez.meanRwdDffImgMotor, [-.3 .3])


%% save rez
save(fullfile(filePath, 'Matfiles', strcat(header, '_dff_evtAligned_regionMask_rs.mat')), 'rezRS');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [evtAlignedDffTs, evtAlignedTs] = alignToEvent(dffTs, eventTimeToAlign, frameT, timeWin)
        win = timeWin + eventTimeToAlign;

        if min(win) >= min(frameT) && max(win) <= max(frameT)
            frameI = frameT >= min(win) & frameT <= max(win);
            evtAlignedDffTs = dffTs(frameI);
            evtAlignedDffTs = evtAlignedDffTs(:).'; % make it a row vector
            evtAlignedTs = frameT(frameI);
            evtAlignedTs = evtAlignedTs(:).'; % make it a row vector
        else
            warning("The peri-event window goes out of bound!")
            evtAlignedDffTs = [];
            evtAlignedTs = [];
        end
    end


    function [imagStackMaskedMean, imgStackMasked] = apply2DMaskTo3DStack(imgStack, mask2D)
        % Replicate the 2D mask to match the 3D stack dimensions
        replicatedMask = repmat(mask2D, [1, 1, size(imgStack, 3)]);

        % Apply the mask: replace values in the stack where the mask is zero (or false) with NaN
        imgStack(~replicatedMask) = NaN;

        imgStackMasked = imgStack;

        imagStackMaskedMean = squeeze(nanmean(nanmean(imgStackMasked, 1), 2));

    end


    function [mRez, sRez] = trialGroupMeanSem(tbytPsthC, groupC)
        % 'tbytPsthC' contains trial-by-trial data in each cell (it is assumed that data are temporally aligned across trials already).
        % 'groupC' contains logical for each grouping, logicals are assumed to be of same lengths as the number of trials.
        % 'mRez' returns the group mean for each group logical per row.
        % 'sRez' returns the group sem for each group logical per row.

        % tbytPsthC = stimOnDff_rs_itp;
        % groupC = {[tbytDat.rewardTrI], [tbytDat.punishTrI]};

        % sanity check 1: all trial numbers must match!
        lenT = length(tbytPsthC);
        trialIC = cell2mat(cellfun(@length, groupC, 'UniformOutput', false));
        if length(unique([lenT, trialIC]))~=1
            error("The length of trials and length(s) of group logics do not match!")
        end

        % sanity check 2:
        if length(unique(cell2mat(cellfun(@length, tbytPsthC, 'UniformOutput', false))))~=1
            error("Some trials have different lengths, e.g., temporally misaligned!")
        end

        nGroups = length(groupC);

        mRez = [];
        sRez = [];

        for gg = 1:nGroups
            gDat = cell2mat(tbytPsthC(logical(groupC{gg}), 1));
            [gMean, ~, gSem] = meanstdsem(gDat);
            mRez = [mRez; gMean];
            sRez = [sRez; gSem];
        end
    end


    function accuracy = logireg(X, Y)
        % % X: observations-by-feature
        % %   e.g., cell2mat(stimOnDff_rs_itp);
        % % Y: label

        Y = Y(:);

        % Determine the number of samples for each class in the test set
        numClass0 = sum(Y == 0);
        numClass1 = sum(Y == 1);
        testSetSize = min(numClass0, numClass1) / 2; % Or any other criterion

        % Randomly sample for the test set
        indices0 = find(Y == 0);
        indices1 = find(Y == 1);
        testIndices0 = randsample(indices0, testSetSize);
        testIndices1 = randsample(indices1, testSetSize);

        % Combine test samples and create the test set
        testInd = [testIndices0; testIndices1];
        testInd = testInd(randperm(length(testInd)));

        % Use the remaining data as the training set
        trainInd = setdiff(1:length(Y), testInd);

        accuracy = nan(1, size(X, 2));
        % Training
        for jj = 1:size(X, 2) % iterate features
            model = fitglm(X(trainInd, jj), Y(trainInd), 'Distribution', 'binomial');
            % Testing
            Y_pred = predict(model, X(testInd,jj)) > 0.5; % Thresholding at 0.5
            accuracy(1, jj) = sum(Y_pred == Y(testInd)) / length(Y_pred);
        end

    end

    function accuracy = multiClass_svm_peth(Xs, y, resample)
        % Xs = [cell2mat(hitFstLickDffSsItp); cell2mat(faFstLickDffSsItp); cell2mat(postStimLickDffSsItp); cell2mat(itiLickDffSsItp)];
        % y = [ones(size(hitFstLickDffSsItp, 1), 1)*1; ones(size(faFstLickDffSsItp, 1), 1)*2; ones(size(postStimLickDffSsItp, 1), 1)*3; ones(size(itiLickDffSsItp, 1), 1)*4];
        % resample = 10;


        accuracy = zeros(resample, size(Xs, 2));

        % Find the minimum class size
        minSize = min(histcounts(y));

        for jj = 1:size(Xs, 2)

            X = Xs(:, jj);

            for rs = 1:resample
                % Create new variables for balanced data
                X_balanced = [];
                y_balanced = [];

                % Balance each class
                uniqueClasses = unique(y);
                for i = 1:length(uniqueClasses)
                    classIndices = find(y == uniqueClasses(i));
                    % Randomly select 'minSize' samples from each class
                    randIndices = randsample(classIndices, minSize);
                    X_balanced = [X_balanced; X(randIndices, :)];
                    y_balanced = [y_balanced; y(randIndices)];
                end

                % Split Data into Training and Testing
                cv = cvpartition(size(X_balanced,1),'HoldOut',0.3);
                idx = cv.test;

                % Separate to training and test data
                XTrain = X_balanced(~idx,:);
                YTrain = y_balanced(~idx,:);
                XTest  = X_balanced(idx,:);
                YTest  = y_balanced(idx,:);

                % Step 4: Train SVM Classifier
                SVMModel = fitcecoc(XTrain, YTrain);
                %SVMModel = fitcnb(XTrain, YTrain);

                % Step 5: Test the Classifier
                YPred = predict(SVMModel, XTest);

                % Evaluate performance
                accuracy(rs, jj) = sum(YPred == YTest) / length(YTest);
                fprintf('Finished resample #%d of bin #%d with accuracy: %.2f%%\n', rs, jj, accuracy(rs, jj) * 100);
            end
        end

    end



end


