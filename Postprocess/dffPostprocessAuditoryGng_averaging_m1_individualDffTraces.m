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
xCorrLagMat = cell2mat(rez.xcorrDffLick.m1(:, 2)); 
xCorrLag = xCorrLagMat(1, :); 

[rezM1.meanXcorrDffLick, ~, rezM1.semXcorrDffLick] = meanstdsem(cell2mat(rez.xcorrDffLick.m1(:, 1)));
xcorrFig = plotMeanSemColor(rezM1.meanXcorrDffLick, rezM1.semXcorrDffLick, xCorrLag, m1_color, {'xcorr m1'});
xlabel('Time (s)'); ylabel('xcorr'); set(gca, 'TickDir', 'out', 'XTick', -2000:1000:2000); grid on
print(fullfile(filePath, 'Figure', 'xcorr_Gng_lick_dff_m1.pdf'), '-dpdf', '-vector', '-bestfit')
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
[stimOnDff_m1_itp, stimOn_timepts] = temporalAlignInterp1(rez.stimOnDffC.m1(:, 1), rez.stimOnDffC.m1(:, 2), 0.001);
[rezM1.meanStimOnDff, ~, rezM1.semStimOnDff] = meanstdsem(cell2mat(stimOnDff_m1_itp));
[stimOnDff_m1_go_itp, ~] = temporalAlignInterp1(rez.stimOnDffC.m1(goI, 1), rez.stimOnDffC.m1(goI, 2), 0.001);
[stimOnDff_m1_nogo_itp, ~] = temporalAlignInterp1(rez.stimOnDffC.m1(nogoI, 1), rez.stimOnDffC.m1(nogoI, 2), 0.001);

[stimOnDff_m1_hit_itp, ~] = temporalAlignInterp1(rez.stimOnDffC.m1(hitI, 1), rez.stimOnDffC.m1(hitI, 2), 0.001);
[stimOnDff_m1_fa_itp, ~] = temporalAlignInterp1(rez.stimOnDffC.m1(faI, 1), rez.stimOnDffC.m1(faI, 2), 0.001);
[stimOnDff_m1_miss_itp, ~] = temporalAlignInterp1(rez.stimOnDffC.m1(missI, 1), rez.stimOnDffC.m1(missI, 2), 0.001);
[stimOnDff_m1_cr_itp, ~] = temporalAlignInterp1(rez.stimOnDffC.m1(crI, 1), rez.stimOnDffC.m1(crI, 2), 0.001);

% primary motor (stim onset) go no-go
[rezM1.meanStimOnDffGng, rezM1.semStimOnDffGng] = trialGroupMeanSem(stimOnDff_m1_itp, {[tbytDat.rewardTrI], [tbytDat.punishTrI]});
stimOnFig = plotMeanSem(rezM1.meanStimOnDffGng, rezM1.semStimOnDffGng, stimOn_timepts, {'Go', 'NoGo'});
title("Go Nogo m1 (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_cueOn_Gng_m1.pdf'), '-dpdf', '-vector', '-bestfit')
close(stimOnFig)
fprintf("Saved tone(stim)-aligned dff plot!\n")

% primary motor (stim onset) hit vs cr
[rezM1.meanStimOnDffHitCr, rezM1.semStimOnDffHitCr] = trialGroupMeanSem(stimOnDff_m1_itp, {hitI, crI});
stimOnHitCrFig = plotMeanSem(rezM1.meanStimOnDffHitCr, rezM1.semStimOnDffHitCr, stimOn_timepts, {'Hit', 'CR'});
title("Hit vs Cr m1 (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_cueOn_HitCr_m1.pdf'), '-dpdf', '-vector', '-bestfit')
close(stimOnHitCrFig)
fprintf("Saved tone(stim)-aligned dff plot comparing Hit vs Cr!\n")

% primary motor (stim onset) hit vs FA
[rezM1.meanStimOnDffHitFA, rezM1.semStimOnDffHitFA] = trialGroupMeanSem(stimOnDff_m1_itp, {hitI, faI});
stimOnHitFaFig = plotMeanSem(rezM1.meanStimOnDffHitFA, rezM1.semStimOnDffHitFA, stimOn_timepts, {'Hit', 'FA'});
title("Hit vs FA m1 (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_cueOn_HitFA_m1.pdf'), '-dpdf', '-vector', '-bestfit')
close(stimOnHitFaFig)
fprintf("Saved tone(stim)-aligned dff plot comparing Hit vs FA!\n")

% primary motor (stim onset) Miss vs CR
[rezM1.meanStimOnDffMissCR, rezM1.semStimOnDffMissCR] = trialGroupMeanSem(stimOnDff_m1_itp, {missI, crI});
stimOnMissCrFig = plotMeanSem(rezM1.meanStimOnDffMissCR, rezM1.semStimOnDffMissCR, stimOn_timepts, {'Miss', 'CR'});
title("Miss vs CR m1 (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_cueOn_MissCR_m1.pdf'), '-dpdf', '-vector', '-bestfit')
close(stimOnMissCrFig)
fprintf("Saved tone(stim)-aligned dff plot comparing Miss vs CR!\n")

% primary motor (stim onset) FA vs CR
[rezM1.meanStimOnDffFaCR, rezM1.semStimOnDffFaCR] = trialGroupMeanSem(stimOnDff_m1_itp, {faI, crI});
stimOnFaCrFig = plotMeanSem(rezM1.meanStimOnDffFaCR, rezM1.semStimOnDffFaCR, stimOn_timepts, {'FA', 'CR'});
title("FA vs CR m1 (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_cueOn_FaCR_m1.pdf'), '-dpdf', '-vector', '-bestfit')
close(stimOnFaCrFig)
fprintf("Saved tone(stim)-aligned dff plot comparing Miss vs CR!\n")

% m1 (stim onset) go vs no-go (train classifier)
[rezM1.stimOnGoNogoSvm, rezM1.stimOnGoNogoSvmTs, rezM1.stimOnGoNogoNb, rezM1.stimOnGoNogoNbTs] = trainDffClassifier(stimOnDff_m1_go_itp, stimOnDff_m1_nogo_itp, stimOn_timepts, 50, 50, 10); 
% figure; plot(rezM1.stimOnGoNogoSvmTs, mean(rezM1.stimOnGoNogoSvm))

% m1 (stim onset) hit vs CR (train classifier)
[rezM1.stimOnHitCrSvm, rezM1.stimOnHitCrSvmTs, rezM1.stimOnHitCrNb, rezM1.stimOnHitCrNbTs] = trainDffClassifier(stimOnDff_m1_hit_itp, stimOnDff_m1_cr_itp, stimOn_timepts, 50, 50, 10); 
% figure; plot(rezM1.stimOnHitCrSvmTs, mean(rezM1.stimOnHitCrSvm))

% m1 (stim onset) hit vs FA (train classifier)
[rezM1.stimOnHitFaSvm, rezM1.stimOnHitFaSvmTs, rezM1.stimOnHitFaNb, rezM1.stimOnHitFaNbTs] = trainDffClassifier(stimOnDff_m1_hit_itp, stimOnDff_m1_fa_itp, stimOn_timepts, 50, 50, 10); 
% figure; plot(rezM1.stimOnHitFaSvmTs, mean(rezM1.stimOnHitFaSvm))

% m1 (stim onset) hit vs Miss (train classifier)
[rezM1.stimOnHitMissSvm, rezM1.stimOnHitMissSvmTs, rezM1.stimOnHitMissNb, rezM1.stimOnHitMissNbTs] = trainDffClassifier(stimOnDff_m1_hit_itp, stimOnDff_m1_miss_itp, stimOn_timepts, 50, 50, 10); 
% figure; plot(rezM1.stimOnHitMissSvmTs, mean(rezM1.stimOnHitMissSvm))

% m1 (stim onset) Miss vs CR (train classifier)
[rezM1.stimOnMissCrSvm, rezM1.stimOnMissCrSvmTs, rezM1.stimOnMissCrNb, rezM1.stimOnMissCrNbTs] = trainDffClassifier(stimOnDff_m1_miss_itp, stimOnDff_m1_cr_itp, stimOn_timepts, 50, 50, 10); 
% figure; plot(rezM1.stimOnMissCrSvmTs, mean(rezM1.stimOnMissCrSvm))

fprintf("Completed cueOn-aligned dff and classification of trial types!\n")

%% iti licks (m1) Technically there's no ITI licks with the retractable spout
% itiLickDffm1 = flatCell(rez.itiLickDff.m1); % unnest the cell arrays
% itiLickDffm1Ts = cellfun(@(a) linspace(-1, 1, length(a)), itiLickDffm1(:, 2), 'UniformOutput', false);
% [itiLickDffm1Itp, itiLickDffm1ItpTs] = temporalAlignInterp1(itiLickDffm1 (:, 1), itiLickDffm1Ts);
% [rez.meanItiLickDff.m1, ~, rez.semItiLickDff.m1] = meanstdsem(cell2mat(itiLickDffm1Itp));
% 
% plotMeanSemColor(rez.meanItiLickDff.m1, rez.semItiLickDff.m1, itiLickDffm1ItpTs, m1_color, {'m1'});
% title("lick during ITI");
% xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); ylim([-.3 .3])

%% postStim licks (m1)
postStimLickDffm1 = flatCell(rez.postStimLickDff.m1); % unnest the cell arrays
postStimLickDffm1Ts = cellfun(@(a) linspace(-1, 1, length(a)), postStimLickDffm1(:, 2), 'UniformOutput', false);
[postStimLickDffm1Itp, postStimLickDffm1ItpTs] = temporalAlignInterp1(postStimLickDffm1(:, 1), postStimLickDffm1Ts, 0.001);
[rezM1.meanPostStimLickDff, ~, rezM1.semPostStimLickDff] = meanstdsem(cell2mat(postStimLickDffm1Itp));

postStimLickFig = plotMeanSemColor(rezM1.meanPostStimLickDff, rezM1.semPostStimLickDff, postStimLickDffm1ItpTs, m1_color, {'m1'});
title("postStim licks");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); %ylim([-.3 .3])
print(fullfile(filePath, 'Figure', 'dff_postStimLick_Gng_m1.pdf'), '-dpdf', '-vector', '-bestfit')
close(postStimLickFig)
fprintf("Saved postStimLick-aligned dff plot!\n")

%% hit firstlicks (m1)
hitFstLickDffm1 = cellWithNonEmptyColumns(rez.hitDffFirstLick.m1);
hitFstLickDffm1Ts = cellfun(@(a) linspace(-1, 1, length(a)), hitFstLickDffm1(:, 2), 'UniformOutput', false);
[hitFstLickDffm1Itp, hitFstLickDffm1ItpTs] = temporalAlignInterp1(hitFstLickDffm1(:, 1), hitFstLickDffm1Ts, 0.001);
[rezM1.meanHitFstLickDff, ~, rezM1.semHitFstLickDff] = meanstdsem(cell2mat(hitFstLickDffm1Itp));

hitFstLickFig = plotMeanSemColor(rezM1.meanHitFstLickDff, rezM1.semHitFstLickDff, hitFstLickDffm1ItpTs, m1_color, {'m1'});
title("Lick (first) Hit trial");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_hitFirstLick_Gng_m1.pdf'), '-dpdf', '-vector', '-bestfit')
close(hitFstLickFig)
fprintf("Saved hitFirstLick-aligned dff plot!\n")

%% first water (m1)
fstWaterDffm1 = cellWithNonEmptyColumns(rez.firstWater.m1);
fstWaterDffm1Ts = cellfun(@(a) linspace(-1, 1, length(a)), fstWaterDffm1(:, 2), 'UniformOutput', false);
[fstWaterDffm1Itp, fstWaterDffm1ItpTs] = temporalAlignInterp1(fstWaterDffm1 (:, 1), fstWaterDffm1Ts, 0.001);
[rezM1.meanFstWaterDff, ~, rezM1.semFstWaterDff] = meanstdsem(cell2mat(fstWaterDffm1Itp));

firstWaterFig = plotMeanSemColor(rezM1.meanFstWaterDff, rezM1.semFstWaterDff, fstWaterDffm1ItpTs, m1_color, {'m1'});
title("Water (first)");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); grid on;
print(fullfile(filePath, 'Figure', 'dff_firstWater_m1.pdf'), '-dpdf', '-vector', '-bestfit')
close(firstWaterFig)
fprintf("Saved firstWater-aligned dff plot!\n")

%% first air (m1)
if isfield(rez, 'firstAirpuff')
    fstAirDffm1 = cellWithNonEmptyColumns(rez.firstAirpuff.m1);
    fstAirDffm1Ts = cellfun(@(a) linspace(-1, 1, length(a)), fstAirDffm1(:, 2), 'UniformOutput', false);
    [fstAirDffm1Itp, fstAirDffm1ItpTs] = temporalAlignInterp1(fstAirDffm1 (:, 1), fstAirDffm1Ts, 0.001);
    [rezM1.meanFstAirDff, ~, rezM1.semFstAirDff] = meanstdsem(cell2mat(fstAirDffm1Itp));

    fstAirFig = plotMeanSemColor(rezM1.meanFstAirDff, rezM1.semFstAirDff, fstAirDffm1ItpTs, m1_color, {'m1'});
    title("Air (first)");
    xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); %ylim([-0.8 0.8])

    print(fullfile(filePath, 'Figure', 'dff_firstAirpuff_m1.pdf'), '-dpdf', '-vector', '-bestfit')
    close(fstAirFig)
    fprintf("Saved firstAir-aligned dff plot!\n")
end

%% FA first licks (m1)
faFstLickDffm1 = cellWithNonEmptyColumns(rez.faDffFirstLick.m1);
faFstLickDffm1Ts = cellfun(@(a) linspace(-1, 1, length(a)), faFstLickDffm1(:, 2), 'UniformOutput', false);
[faFstLickDffm1Itp, faFstLickDffm1ItpTs] = temporalAlignInterp1(faFstLickDffm1 (:, 1), faFstLickDffm1Ts, 0.001);
[rezM1.meanFaFstLickDff, ~, rezM1.semFaFstLickDff] = meanstdsem(cell2mat(faFstLickDffm1Itp));

faFstLickFig = plotMeanSem(rezM1.meanFaFstLickDff, rezM1.semFaFstLickDff, faFstLickDffm1ItpTs, {'m1'});
title("Lick (first) Fa trial");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); 
print(fullfile(filePath, 'Figure', 'dff_faFstLickFig_m1.pdf'), '-dpdf', '-vector', '-bestfit')
close(faFstLickFig)

%% Compare licks during Hit vs FA vs postStim (m1)
% Without normalization
hitFaPostStimFig = plotMeanSem([rezM1.meanHitFstLickDff; rezM1.meanFaFstLickDff; rezM1.meanPostStimLickDff], ...
    [rezM1.semHitFstLickDff; rezM1.semFaFstLickDff; rezM1.semPostStimLickDff], ...
    hitFstLickDffm1ItpTs, {'Hit', 'FA', 'postStim'}); 
title("Lick Hit vs FA vs postStim")
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); grid on
print(fullfile(filePath, 'Figure', 'dff_lick_Hit_Fa_postStim_m1.pdf'), '-dpdf', '-vector', '-bestfit')
close(hitFaPostStimFig)

% With normalization
rezM1.meanHitFstLickDffBaseSub = baseSubNorm(rezM1.meanHitFstLickDff, hitFstLickDffm1ItpTs, [-1 -0.5]);
rezM1.meanFaFstLickDffBaseSub = baseSubNorm(rezM1.meanFaFstLickDff, faFstLickDffm1ItpTs, [-1 -0.5]);
rezM1.meanPostStimLickDffBaseSub = baseSubNorm(rezM1.meanPostStimLickDff, postStimLickDffm1ItpTs, [-1 -0.5]);

hitFaPostStimFig_baseSub = plotMeanSem([rezM1.meanHitFstLickDffBaseSub; rezM1.meanFaFstLickDffBaseSub; rezM1.meanPostStimLickDffBaseSub], ...
    [rezM1.semHitFstLickDff; rezM1.semFaFstLickDff; rezM1.semPostStimLickDff], ...
    hitFstLickDffm1ItpTs, {'Hit', 'FA', 'postStim'}); 
title("Lick Hit vs FA vs ITI vs postStim (baseSub)")
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); grid on
print(fullfile(filePath, 'Figure', 'dff_lick_Hit_Fa_postStim_baseSub_m1.pdf'), '-dpdf', '-vector', '-bestfit')
close(hitFaPostStimFig_baseSub)

%% train svm to classify licks at different task epochs
% m1 (stim onset) go vs no-go (train classifier)
[rezM1.lickHitFaSvm, rezM1.lickHitFaSvmTs, rezM1.lickHitFaNb, rezM1.lickHitFaNbTs, X_hit, X_fa] = trainDffClassifier(hitFstLickDffm1Itp, faFstLickDffm1Itp, hitFstLickDffm1ItpTs, 50, 50, 10); 
% figure; plot(rezM1.lickHitFaSvmTs, mean(rezM1.lickHitFaSvm))


% imagesc the hit and fa trials dff
X_hit_fig = imagescWithTimeInfo(smooth2a(X_hit, 0, 1), rezM1.lickHitFaSvmTs); 
set(gca, 'XTick', -2:.5:4, 'TickDir', 'out');
clim([-0.5 0.5])
print(fullfile(filePath, 'Figure', 'dff_hit_timeBin_by_trial_m1.pdf'), '-dpdf', '-vector', '-bestfit')
close(X_hit_fig)

X_fa_fig = imagescWithTimeInfo(smooth2a(X_fa, 0, 1), rezM1.lickHitFaSvmTs); 
set(gca, 'XTick', -2:.5:4, 'TickDir', 'out');
clim([-0.7 0.7])
print(fullfile(filePath, 'Figure', 'dff_fa_timeBin_by_trial_m1.pdf'), '-dpdf', '-vector', '-bestfit')
close(X_fa_fig)

%% Miss CueOn (m1)
% primary motor (stim onset) go no-go
missCueI = [tbytDat.rewardTrI] & ~waterI & ~lickI;  
[rezM1.meanStimOnMissDff, rezM1.semStimOnMissDff] = trialGroupMeanSem(stimOnDff_m1_itp, {missCueI});
stimOnMissFig = plotMeanSem(rezM1.meanStimOnMissDff, rezM1.semStimOnMissDff, stimOn_timepts, {'m1 miss'});
title("Go Miss m1 (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_cueOn_goMiss_m1.pdf'), '-dpdf', '-vector', '-bestfit')
close(stimOnMissFig)
fprintf("Saved tone(stim)-aligned dff of Miss trials!\n")

%% CR CueOn (m1)
% primary motor (stim onset) go no-go
crCueI = [tbytDat.punishTrI] & ~airpuffI & ~lickI;  
[rezM1.meanStimOnCrDff, rezM1.semStimOnCrDff] = trialGroupMeanSem(stimOnDff_m1_itp, {crCueI});
stimOnCrFig = plotMeanSem(rezM1.meanStimOnCrDff, rezM1.semStimOnCrDff, stimOn_timepts, {'m1 correct rejection'});
title("NoGo CR m1 (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_cueOn_nogoCR_m1.pdf'), '-dpdf', '-vector', '-bestfit')
close(stimOnCrFig)
fprintf("Saved tone(stim)-aligned dff of Correct Rejection trials!\n")

% %% plot stim onset PETHs together
% plotMeanSem([rez.mStimOnDffm1; rez.mStimOnDffM1; rez.mStimOnDffS1; rez.mStimOnDffRs; rez.mStimOnDffV], ...
%     [rez.semStimOnDffm1; rez.semStimOnDffM1; rez.semStimOnDffS1; rez.semStimOnDffRs; rez.semStimOnDffV], ...
%     stimOn_timepts, ...
%     {'m1', 'M1', 'S1', 'M1', 'V'});
% xlabel('Time (s)')
% ylabel('DFF')
% set(gca, 'XTick', -2:1:4)
% xlim([-2 2])
% print(fullfile(filePath, 'Figure', 'dff_cueOn_Collect.pdf'), '-dpdf', '-vector', '-bestfit')

%% plot classification accuracy
% plotMeanSem(smooth2a([rez.logitGngCueOnm1; rez.logitGngCueOnM1; rez.logitGngCueOnS1; rez.logitGngCueOnRs; rez.logitGngCueOnV], 0, 2), ...
%     zeros(5, length(stimOn_timepts)), stimOn_timepts, {'m1', 'M1', 'S1', 'M1', 'V'});
% title("Cross-validated classification accuracy")
% xlabel('Time (s)')
% ylabel('Classification Accuracy')
% set(gca, 'XTick', -2:1:4, 'YTick', -1:0.1:1, 'TickDir', 'out')
% xlim([-2 2])
% ylim([0.45 0.85])
% print(fullfile(filePath, 'Figure', 'logisticRegression_Gng_cueOn_Collect.pdf'), '-dpdf', '-vector', '-bestfit')

% baseline dff image
% rez.meanBaseDffImgMotor = nanmean(cell2mat(reshape(rwdDffC_motor(:, 3), [1, 1, length(tbytDat)])), 3);
% imageFrameWithNaNsEdgeMap(rez.meanBaseDffImgMotor, [-.3 .3], dom1alMaps.edgeMap, [0, 0, 1], 0.3)
% title('mean baseline m1')
% set(gca, 'XTickLabel', {})
% set(gca, 'YTickLabel', {})
% print(fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA003/DA003_083023/Figure', 'dff_base_m1_Img.pdf'), '-dpdf', '-vector', '-bestfit')
% 
% % reward dff image
% rez.meanRwdDffImgMotor = nanmean(cell2mat(reshape(rwdDffC_motor(:, 4), [1, 1, length(tbytDat)])), 3);
% imageFrameWithNaNsEdgeMap(rez.meanRwdDffImgMotor, [-.3 .3], dom1alMaps.edgeMap, [0, 0, 1], 0.3)
% print(fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA003/DA003_083023/Figure', 'dff_reward_m1_Img.pdf'), '-dpdf', '-vector', '-bestfit')
% title('mean reward m1')
% set(gca, 'XTickLabel', {})
% set(gca, 'YTickLabel', {})
% print(fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA003/DA003_083023/Figure', 'dff_reward_m1_Img.pdf'), '-dpdf', '-vector', '-bestfit')

%showAllenEdgesOnTransformedWF(rez.meanBaseDffImgMotor, dom1alMaps, 0, 1)

%imageFrameWithNaNs(rez.meanRwdDffImgMotor, [-.3 .3])


%% save rez
save(fullfile(filePath, 'Matfiles', strcat(header, '_dff_evtAligned_regionMask_m1.mat')), 'rezM1');

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

        % tbytPsthC = stimOnDff_m1_itp;
        % groupC = {[tbytDat.rewardTrI], [tbytDat.punishTrI]};

        % sanity check 1: all trial numbem1 must match!
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
        % %   e.g., cell2mat(stimOnDff_m1_itp);
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

            for m1 = 1:resample
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
                accuracy(m1, jj) = sum(YPred == YTest) / length(YTest);
                fprintf('Finished resample #%d of bin #%d with accuracy: %.2f%%\n', m1, jj, accuracy(m1, jj) * 100);
            end
        end

    end



end


