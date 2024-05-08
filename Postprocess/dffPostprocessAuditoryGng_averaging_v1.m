function dffPostprocessAuditoryGng_averaging_v1(filePath)
%filePath = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA008/DA008_101823';
% This function performs further analyses and plotting using data in
% 'rez'(outcome from 'dffPostprocessGng_averaging.m').
%

%% load rez tbytDat
[~, header] = fileparts(filePath);
load(fullfile(filePath, 'Matfiles', strcat(header, '_dff_evtAligned_regionMask.mat')), 'rez');
load(fullfile(filePath, 'Matfiles', strcat(header, '_tbytDat_parseGng')), 'tbytDat')

% user variables
v1_color = [0 240 240]./255;

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
xCorrLagMat = cell2mat(rez.xcorrDffLick.V1(:, 2)); 
xCorrLag = xCorrLagMat(1, :); 

[rezV1.meanXcorrDffLick, ~, rezV1.semXcorrDffLick] = meanstdsem(cell2mat(rez.xcorrDffLick.V1(:, 1)));
xcorrFig = plotMeanSemColor(rezV1.meanXcorrDffLick, rezV1.semXcorrDffLick, xCorrLag, v1_color, {'xcorr v1'});
xlabel('Time (s)'); ylabel('xcorr'); set(gca, 'TickDir', 'out', 'XTick', -2000:1000:2000); grid on
print(fullfile(filePath, 'Figure', 'xcorr_Gng_lick_dff_v1.pdf'), '-dpdf', '-vector', '-bestfit')
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
%     xend = x1+length(rez.dffsOnTfItpV1{t, 1})-1;
% 
%     plot(x1:x1+length(rez.dffsOnTfItpV1{t, 1})-1, rez.dffsOnTfItpV1{t, 1}, 'Color', v1_color, 'LineWidth', 1);
%     vertline(licks, [0 1], 'lineWidth', 1, 'alpha', 0.5);
%     vertline(xend-length(rez.dffsOnTfItpV1{t, 1})/2, ylim, 'lineWidth', 0.75, 'lineStyle', ':', 'alpha', 0.5);
% 
%     x1 = x1 + length(rez.dffsOnTfItpV1{t, 1}) + 1;
% end
% title("dff and lick traces");
% set(gca, 'TickDir', 'out');
% axis tight
% print(fullfile(filePath, 'Figure', 'dff_lickbouts_v1_example5trials.pdf'), '-dpdf', '-vector', '-bestfit')

%% Stim onset aligned activity (v1)
% align dffs take the mean and sem
% primary motor (stim onset)
[stimOnDff_v1_itp, stimOn_timepts] = temporalAlignInterp1(rez.stimOnDffC.v1(:, 1), rez.stimOnDffC.v1(:, 2), 0.001);
[rezV1.meanStimOnDff, ~, rezV1.semStimOnDff] = meanstdsem(cell2mat(stimOnDff_v1_itp));
[stimOnDff_v1_go_itp, ~] = temporalAlignInterp1(rez.stimOnDffC.v1(goI, 1), rez.stimOnDffC.v1(goI, 2), 0.001);
[stimOnDff_v1_nogo_itp, ~] = temporalAlignInterp1(rez.stimOnDffC.v1(nogoI, 1), rez.stimOnDffC.v1(nogoI, 2), 0.001);

[stimOnDff_v1_hit_itp, ~] = temporalAlignInterp1(rez.stimOnDffC.v1(hitI, 1), rez.stimOnDffC.v1(hitI, 2), 0.001);
[stimOnDff_v1_fa_itp, ~] = temporalAlignInterp1(rez.stimOnDffC.v1(faI, 1), rez.stimOnDffC.v1(faI, 2), 0.001);
[stimOnDff_v1_miss_itp, ~] = temporalAlignInterp1(rez.stimOnDffC.v1(missI, 1), rez.stimOnDffC.v1(missI, 2), 0.001);
[stimOnDff_v1_cr_itp, ~] = temporalAlignInterp1(rez.stimOnDffC.v1(crI, 1), rez.stimOnDffC.v1(crI, 2), 0.001);

% primary motor (stim onset) go no-go
[rezV1.meanStimOnDffGng, rezV1.semStimOnDffGng] = trialGroupMeanSem(stimOnDff_v1_itp, {[tbytDat.rewardTrI], [tbytDat.punishTrI]});
stimOnFig = plotMeanSem(rezV1.meanStimOnDffGng, rezV1.semStimOnDffGng, stimOn_timepts, {'Go', 'NoGo'});
title("Go Nogo v1 (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_cueOn_Gng_v1.pdf'), '-dpdf', '-vector', '-bestfit')
close(stimOnFig)
fprintf("Saved tone(stim)-aligned dff plot!\n")

% primary motor (stim onset) hit vs cr
[rezV1.meanStimOnDffHitCr, rezV1.semStimOnDffHitCr] = trialGroupMeanSem(stimOnDff_v1_itp, {hitI, crI});
stimOnHitCrFig = plotMeanSem(rezV1.meanStimOnDffHitCr, rezV1.semStimOnDffHitCr, stimOn_timepts, {'Hit', 'CR'});
title("Hit vs Cr v1 (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_cueOn_HitCr_v1.pdf'), '-dpdf', '-vector', '-bestfit')
close(stimOnHitCrFig)
fprintf("Saved tone(stim)-aligned dff plot comparing Hit vs Cr!\n")

% primary motor (stim onset) hit vs FA
[rezV1.meanStimOnDffHitFA, rezV1.semStimOnDffHitFA] = trialGroupMeanSem(stimOnDff_v1_itp, {hitI, faI});
stimOnHitFaFig = plotMeanSem(rezV1.meanStimOnDffHitFA, rezV1.semStimOnDffHitFA, stimOn_timepts, {'Hit', 'FA'});
title("Hit vs FA v1 (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_cueOn_HitFA_v1.pdf'), '-dpdf', '-vector', '-bestfit')
close(stimOnHitFaFig)
fprintf("Saved tone(stim)-aligned dff plot comparing Hit vs FA!\n")

% primary motor (stim onset) Miss vs CR
[rezV1.meanStimOnDffMissCR, rezV1.semStimOnDffMissCR] = trialGroupMeanSem(stimOnDff_v1_itp, {missI, crI});
stimOnMissCrFig = plotMeanSem(rezV1.meanStimOnDffMissCR, rezV1.semStimOnDffMissCR, stimOn_timepts, {'Miss', 'CR'});
title("Miss vs CR v1 (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_cueOn_MissCR_v1.pdf'), '-dpdf', '-vector', '-bestfit')
close(stimOnMissCrFig)
fprintf("Saved tone(stim)-aligned dff plot comparing Miss vs CR!\n")

% primary motor (stim onset) FA vs CR
[rezV1.meanStimOnDffFaCR, rezV1.semStimOnDffFaCR] = trialGroupMeanSem(stimOnDff_v1_itp, {faI, crI});
stimOnFaCrFig = plotMeanSem(rezV1.meanStimOnDffFaCR, rezV1.semStimOnDffFaCR, stimOn_timepts, {'FA', 'CR'});
title("FA vs CR v1 (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_cueOn_FaCR_v1.pdf'), '-dpdf', '-vector', '-bestfit')
close(stimOnFaCrFig)
fprintf("Saved tone(stim)-aligned dff plot comparing Miss vs CR!\n")

% v1 (stim onset) go vs no-go (train classifier)
[rezV1.stimOnGoNogoSvm, rezV1.stimOnGoNogoSvmTs, rezV1.stimOnGoNogoNb, rezV1.stimOnGoNogoNbTs] = trainDffClassifier(stimOnDff_v1_go_itp, stimOnDff_v1_nogo_itp, stimOn_timepts, 50, 50, 10); 
% figure; plot(rezV1.stimOnGoNogoSvmTs, mean(rezV1.stimOnGoNogoSvm))

% v1 (stim onset) hit vs CR (train classifier)
[rezV1.stimOnHitCrSvm, rezV1.stimOnHitCrSvmTs, rezV1.stimOnHitCrNb, rezV1.stimOnHitCrNbTs] = trainDffClassifier(stimOnDff_v1_hit_itp, stimOnDff_v1_cr_itp, stimOn_timepts, 50, 50, 10); 
% figure; plot(rezV1.stimOnHitCrSvmTs, mean(rezV1.stimOnHitCrSvm))

% v1 (stim onset) hit vs FA (train classifier)
[rezV1.stimOnHitFaSvm, rezV1.stimOnHitFaSvmTs, rezV1.stimOnHitFaNb, rezV1.stimOnHitFaNbTs] = trainDffClassifier(stimOnDff_v1_hit_itp, stimOnDff_v1_fa_itp, stimOn_timepts, 50, 50, 10); 
% figure; plot(rezV1.stimOnHitFaSvmTs, mean(rezV1.stimOnHitFaSvm))

% v1 (stim onset) hit vs Miss (train classifier)
[rezV1.stimOnHitMissSvm, rezV1.stimOnHitMissSvmTs, rezV1.stimOnHitMissNb, rezV1.stimOnHitMissNbTs] = trainDffClassifier(stimOnDff_v1_hit_itp, stimOnDff_v1_miss_itp, stimOn_timepts, 50, 50, 10); 
% figure; plot(rezV1.stimOnHitMissSvmTs, mean(rezV1.stimOnHitMissSvm))

% v1 (stim onset) Miss vs CR (train classifier)
[rezV1.stimOnMissCrSvm, rezV1.stimOnMissCrSvmTs, rezV1.stimOnMissCrNb, rezV1.stimOnMissCrNbTs] = trainDffClassifier(stimOnDff_v1_miss_itp, stimOnDff_v1_cr_itp, stimOn_timepts, 50, 50, 10); 
% figure; plot(rezV1.stimOnMissCrSvmTs, mean(rezV1.stimOnMissCrSvm))

fprintf("Completed cueOn-aligned dff and classification of trial types!\n")

%% iti licks (v1) Technically there's no ITI licks with the retractable spout
% itiLickDffv1 = flatCell(rez.itiLickDff.v1); % unnest the cell arrays
% itiLickDffv1Ts = cellfun(@(a) linspace(-1, 1, length(a)), itiLickDffv1(:, 2), 'UniformOutput', false);
% [itiLickDffv1Itp, itiLickDffv1ItpTs] = temporalAlignInterp1(itiLickDffv1 (:, 1), itiLickDffv1Ts);
% [rez.meanItiLickDff.v1, ~, rez.semItiLickDff.v1] = meanstdsem(cell2mat(itiLickDffv1Itp));
% 
% plotMeanSemColor(rez.meanItiLickDff.v1, rez.semItiLickDff.v1, itiLickDffv1ItpTs, v1_color, {'v1'});
% title("lick during ITI");
% xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); ylim([-.3 .3])

%% postStim licks (v1)
postStimLickDffv1 = flatCell(rez.postStimLickDff.v1); % unnest the cell arrays
postStimLickDffv1Ts = cellfun(@(a) linspace(-1, 1, length(a)), postStimLickDffv1(:, 2), 'UniformOutput', false);
[postStimLickDffv1Itp, postStimLickDffv1ItpTs] = temporalAlignInterp1(postStimLickDffv1(:, 1), postStimLickDffv1Ts, 0.001);
[rezV1.meanPostStimLickDff, ~, rezV1.semPostStimLickDff] = meanstdsem(cell2mat(postStimLickDffv1Itp));

postStimLickFig = plotMeanSemColor(rezV1.meanPostStimLickDff, rezV1.semPostStimLickDff, postStimLickDffv1ItpTs, v1_color, {'v1'});
title("postStim licks");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); %ylim([-.3 .3])
print(fullfile(filePath, 'Figure', 'dff_postStimLick_Gng_v1.pdf'), '-dpdf', '-vector', '-bestfit')
close(postStimLickFig)
fprintf("Saved postStimLick-aligned dff plot!\n")

%% hit firstlicks (v1)
hitFstLickDffv1 = cellWithNonEmptyColumns(rez.hitDffFirstLick.v1);
hitFstLickDffv1Ts = cellfun(@(a) linspace(-1, 1, length(a)), hitFstLickDffv1(:, 2), 'UniformOutput', false);
[hitFstLickDffv1Itp, hitFstLickDffv1ItpTs] = temporalAlignInterp1(hitFstLickDffv1(:, 1), hitFstLickDffv1Ts, 0.001);
[rezV1.meanHitFstLickDff, ~, rezV1.semHitFstLickDff] = meanstdsem(cell2mat(hitFstLickDffv1Itp));

hitFstLickFig = plotMeanSemColor(rezV1.meanHitFstLickDff, rezV1.semHitFstLickDff, hitFstLickDffv1ItpTs, v1_color, {'v1'});
title("Lick (first) Hit trial");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_hitFirstLick_Gng_v1.pdf'), '-dpdf', '-vector', '-bestfit')
close(hitFstLickFig)
fprintf("Saved hitFirstLick-aligned dff plot!\n")

%% first water (v1)
fstWaterDffv1 = cellWithNonEmptyColumns(rez.firstWater.v1);
fstWaterDffv1Ts = cellfun(@(a) linspace(-1, 1, length(a)), fstWaterDffv1(:, 2), 'UniformOutput', false);
[fstWaterDffv1Itp, fstWaterDffv1ItpTs] = temporalAlignInterp1(fstWaterDffv1 (:, 1), fstWaterDffv1Ts, 0.001);
[rezV1.meanFstWaterDff, ~, rezV1.semFstWaterDff] = meanstdsem(cell2mat(fstWaterDffv1Itp));

firstWaterFig = plotMeanSemColor(rezV1.meanFstWaterDff, rezV1.semFstWaterDff, fstWaterDffv1ItpTs, v1_color, {'v1'});
title("Water (first)");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); grid on;
print(fullfile(filePath, 'Figure', 'dff_firstWater_v1.pdf'), '-dpdf', '-vector', '-bestfit')
close(firstWaterFig)
fprintf("Saved firstWater-aligned dff plot!\n")

%% first air (v1)
if isfield(rez, 'firstAirpuff')
    fstAirDffv1 = cellWithNonEmptyColumns(rez.firstAirpuff.v1);
    fstAirDffv1Ts = cellfun(@(a) linspace(-1, 1, length(a)), fstAirDffv1(:, 2), 'UniformOutput', false);
    [fstAirDffv1Itp, fstAirDffv1ItpTs] = temporalAlignInterp1(fstAirDffv1 (:, 1), fstAirDffv1Ts, 0.001);
    [rezV1.meanFstAirDff, ~, rezV1.semFstAirDff] = meanstdsem(cell2mat(fstAirDffv1Itp));

    fstAirFig = plotMeanSemColor(rezV1.meanFstAirDff, rezV1.semFstAirDff, fstAirDffv1ItpTs, v1_color, {'v1'});
    title("Air (first)");
    xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); %ylim([-0.8 0.8])

    print(fullfile(filePath, 'Figure', 'dff_firstAirpuff_v1.pdf'), '-dpdf', '-vector', '-bestfit')
    close(fstAirFig)
    fprintf("Saved firstAir-aligned dff plot!\n")
end

%% FA first licks (v1)
faFstLickDffv1 = cellWithNonEmptyColumns(rez.faDffFirstLick.v1);
faFstLickDffv1Ts = cellfun(@(a) linspace(-1, 1, length(a)), faFstLickDffv1(:, 2), 'UniformOutput', false);
[faFstLickDffv1Itp, faFstLickDffv1ItpTs] = temporalAlignInterp1(faFstLickDffv1 (:, 1), faFstLickDffv1Ts, 0.001);
[rezV1.meanFaFstLickDff, ~, rezV1.semFaFstLickDff] = meanstdsem(cell2mat(faFstLickDffv1Itp));

faFstLickFig = plotMeanSem(rezV1.meanFaFstLickDff, rezV1.semFaFstLickDff, faFstLickDffv1ItpTs, {'v1'});
title("Lick (first) Fa trial");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); 
print(fullfile(filePath, 'Figure', 'dff_faFstLickFig_v1.pdf'), '-dpdf', '-vector', '-bestfit')
close(faFstLickFig)

%% Compare licks during Hit vs FA vs postStim (v1)
% Without normalization
hitFaPostStimFig = plotMeanSem([rezV1.meanHitFstLickDff; rezV1.meanFaFstLickDff; rezV1.meanPostStimLickDff], ...
    [rezV1.semHitFstLickDff; rezV1.semFaFstLickDff; rezV1.semPostStimLickDff], ...
    hitFstLickDffv1ItpTs, {'Hit', 'FA', 'postStim'}); 
title("Lick Hit vs FA vs postStim")
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); grid on
print(fullfile(filePath, 'Figure', 'dff_lick_Hit_Fa_postStim_v1.pdf'), '-dpdf', '-vector', '-bestfit')
close(hitFaPostStimFig)

% With normalization
rezV1.meanHitFstLickDffBaseSub = baseSubNorm(rezV1.meanHitFstLickDff, hitFstLickDffv1ItpTs, [-1 -0.5]);
rezV1.meanFaFstLickDffBaseSub = baseSubNorm(rezV1.meanFaFstLickDff, faFstLickDffv1ItpTs, [-1 -0.5]);
rezV1.meanPostStimLickDffBaseSub = baseSubNorm(rezV1.meanPostStimLickDff, postStimLickDffv1ItpTs, [-1 -0.5]);

hitFaPostStimFig_baseSub = plotMeanSem([rezV1.meanHitFstLickDffBaseSub; rezV1.meanFaFstLickDffBaseSub; rezV1.meanPostStimLickDffBaseSub], ...
    [rezV1.semHitFstLickDff; rezV1.semFaFstLickDff; rezV1.semPostStimLickDff], ...
    hitFstLickDffv1ItpTs, {'Hit', 'FA', 'postStim'}); 
title("Lick Hit vs FA vs ITI vs postStim (baseSub)")
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); grid on
print(fullfile(filePath, 'Figure', 'dff_lick_Hit_Fa_postStim_baseSub_v1.pdf'), '-dpdf', '-vector', '-bestfit')
close(hitFaPostStimFig_baseSub)

%% train svm to classify licks at different task epochs
% v1 hit lick vs FA lick (train classifier)
[rezV1.lickHitFaSvm, rezV1.lickHitFaSvmTs, rezV1.lickHitFaNb, rezV1.lickHitFaNbTs, X_hit, X_fa] = trainDffClassifier(hitFstLickDffv1Itp, faFstLickDffv1Itp, hitFstLickDffv1ItpTs, 50, 50, 10); 
% figure; plot(rezV1.lickHitFaSvmTs, mean(rezV1.lickHitFaSvm))

% imagesc the hit and fa trials dff
X_hit_fig = imagescWithTimeInfo(smooth2a(X_hit, 0, 1), rezV1.lickHitFaSvmTs); 
set(gca, 'XTick', -2:.5:4, 'TickDir', 'out');
clim([-0.5 0.5])
print(fullfile(filePath, 'Figure', 'dff_hit_timeBin_by_trial_v1.pdf'), '-dpdf', '-vector', '-bestfit')
close(X_hit_fig)

X_fa_fig = imagescWithTimeInfo(smooth2a(X_fa, 0, 1), rezV1.lickHitFaSvmTs); 
set(gca, 'XTick', -2:.5:4, 'TickDir', 'out');
clim([-0.7 0.7])
print(fullfile(filePath, 'Figure', 'dff_fa_timeBin_by_trial_v1.pdf'), '-dpdf', '-vector', '-bestfit')
close(X_fa_fig)

%% Miss CueOn (v1)
% primary motor (stim onset) go no-go
missCueI = [tbytDat.rewardTrI] & ~waterI & ~lickI;  
[rezV1.meanStimOnMissDff, rezV1.semStimOnMissDff] = trialGroupMeanSem(stimOnDff_v1_itp, {missCueI});
stimOnMissFig = plotMeanSem(rezV1.meanStimOnMissDff, rezV1.semStimOnMissDff, stimOn_timepts, {'v1 miss'});
title("Go Miss v1 (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_cueOn_goMiss_v1.pdf'), '-dpdf', '-vector', '-bestfit')
close(stimOnMissFig)
fprintf("Saved tone(stim)-aligned dff of Miss trials!\n")

%% CR CueOn (v1)
% primary motor (stim onset) go no-go
crCueI = [tbytDat.punishTrI] & ~airpuffI & ~lickI;  
[rezV1.meanStimOnCrDff, rezV1.semStimOnCrDff] = trialGroupMeanSem(stimOnDff_v1_itp, {crCueI});
stimOnCrFig = plotMeanSem(rezV1.meanStimOnCrDff, rezV1.semStimOnCrDff, stimOn_timepts, {'v1 correct rejection'});
title("NoGo CR v1 (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_cueOn_nogoCR_v1.pdf'), '-dpdf', '-vector', '-bestfit')
close(stimOnCrFig)
fprintf("Saved tone(stim)-aligned dff of Correct Rejection trials!\n")

% %% plot stim onset PETHs together
% plotMeanSem([rez.mStimOnDffv1; rez.mStimOnDffV1; rez.mStimOnDffS1; rez.mStimOnDffRs; rez.mStimOnDffV], ...
%     [rez.semStimOnDffv1; rez.semStimOnDffV1; rez.semStimOnDffS1; rez.semStimOnDffRs; rez.semStimOnDffV], ...
%     stimOn_timepts, ...
%     {'v1', 'V1', 'S1', 'V1', 'V'});
% xlabel('Time (s)')
% ylabel('DFF')
% set(gca, 'XTick', -2:1:4)
% xlim([-2 2])
% print(fullfile(filePath, 'Figure', 'dff_cueOn_Collect.pdf'), '-dpdf', '-vector', '-bestfit')

%% plot classification accuracy
% plotMeanSem(smooth2a([rez.logitGngCueOnv1; rez.logitGngCueOnV1; rez.logitGngCueOnS1; rez.logitGngCueOnRs; rez.logitGngCueOnV], 0, 2), ...
%     zeros(5, length(stimOn_timepts)), stimOn_timepts, {'v1', 'V1', 'S1', 'V1', 'V'});
% title("Cross-validated classification accuracy")
% xlabel('Time (s)')
% ylabel('Classification Accuracy')
% set(gca, 'XTick', -2:1:4, 'YTick', -1:0.1:1, 'TickDir', 'out')
% xlim([-2 2])
% ylim([0.45 0.85])
% print(fullfile(filePath, 'Figure', 'logisticRegression_Gng_cueOn_Collect.pdf'), '-dpdf', '-vector', '-bestfit')

% baseline dff image
% rez.meanBaseDffImgMotor = nanmean(cell2mat(reshape(rwdDffC_motor(:, 3), [1, 1, length(tbytDat)])), 3);
% imageFrameWithNaNsEdgeMap(rez.meanBaseDffImgMotor, [-.3 .3], dov1alMaps.edgeMap, [0, 0, 1], 0.3)
% title('mean baseline v1')
% set(gca, 'XTickLabel', {})
% set(gca, 'YTickLabel', {})
% print(fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA003/DA003_083023/Figure', 'dff_base_v1_Img.pdf'), '-dpdf', '-vector', '-bestfit')
% 
% % reward dff image
% rez.meanRwdDffImgMotor = nanmean(cell2mat(reshape(rwdDffC_motor(:, 4), [1, 1, length(tbytDat)])), 3);
% imageFrameWithNaNsEdgeMap(rez.meanRwdDffImgMotor, [-.3 .3], dov1alMaps.edgeMap, [0, 0, 1], 0.3)
% print(fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA003/DA003_083023/Figure', 'dff_reward_v1_Img.pdf'), '-dpdf', '-vector', '-bestfit')
% title('mean reward v1')
% set(gca, 'XTickLabel', {})
% set(gca, 'YTickLabel', {})
% print(fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA003/DA003_083023/Figure', 'dff_reward_v1_Img.pdf'), '-dpdf', '-vector', '-bestfit')

%showAllenEdgesOnTransformedWF(rez.meanBaseDffImgMotor, dov1alMaps, 0, 1)

%imageFrameWithNaNs(rez.meanRwdDffImgMotor, [-.3 .3])


%% save rez
save(fullfile(filePath, 'Matfiles', strcat(header, '_dff_evtAligned_regionMask_v1.mat')), 'rezV1');

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

        % tbytPsthC = stimOnDff_v1_itp;
        % groupC = {[tbytDat.rewardTrI], [tbytDat.punishTrI]};

        % sanity check 1: all trial numbev1 must match!
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
        % %   e.g., cell2mat(stimOnDff_v1_itp);
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

            for v1 = 1:resample
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
                accuracy(v1, jj) = sum(YPred == YTest) / length(YTest);
                fprintf('Finished resample #%d of bin #%d with accuracy: %.2f%%\n', v1, jj, accuracy(v1, jj) * 100);
            end
        end

    end



end


