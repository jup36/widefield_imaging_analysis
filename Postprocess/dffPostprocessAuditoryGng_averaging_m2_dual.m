function dffPostprocessAuditoryGng_averaging_m2_dual(filePath, channel)
%filePath = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA008/DA008_101823';
% This function performs further analyses and plotting using data in
% 'rez'(outcome from 'dffPostprocessGng_averaging.m').
%

%% load rez tbytDat
match_header = regexp(filePath, 'm\d{1,4}_\d{6}', 'match');
header = match_header{1};
load(fullfile(filePath, 'Matfiles', strcat(header, '_', channel, '_dff_evtAligned_regionMask.mat')), 'rez');
load(fullfile(filePath, 'Matfiles', strcat(header, '_tbytDat_parseGng')), 'tbytDat')

% user variables
m2_color = [64 191 255]./255; 

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
xCorrLagMat = cell2mat(rez.xcorrDffLick.m2(:, 2)); 
xCorrLag = xCorrLagMat(1, :); 

[rezM2.meanXcorrDffLick, ~, rezM2.semXcorrDffLick] = meanstdsem(cell2mat(rez.xcorrDffLick.m2(:, 1)));
xcorrFig = plotMeanSemColor(rezM2.meanXcorrDffLick, rezM2.semXcorrDffLick, xCorrLag, m2_color, {'xcorr m2'});
xlabel('Time (s)'); ylabel('xcorr'); set(gca, 'TickDir', 'out', 'XTick', -2000:1000:2000); grid on
print(fullfile(filePath, 'Figure', strcat('xcorr_Gng_lick_dff_m2_', channel, '.pdf')), '-dpdf', '-vector', '-bestfit')
close(xcorrFig)
fprintf("Processed the xcorr between lick and dff, and saved the crosscorrelogram!\n")

%% example lick and dff traces
% exampleTrials = [14];
% x1 = 1;
% figure; hold on;
% for j = 1:length(exampleTrials)
% 
%     t = exampleTrials(j);
%     licks = x1+find(rez.lickOnTfBin{t, 1});
%     xend = x1+length(rez.dffsOnTfItpM2{t, 1})-1;
% 
%     plot(x1:x1+length(rez.dffsOnTfItpM2{t, 1})-1, rez.dffsOnTfItpM2{t, 1}, 'Color', m2_color, 'LineWidth', 1);
%     vertline(licks, [0 1], 'lineWidth', 1, 'alpha', 0.5);
%     vertline(xend-length(rez.dffsOnTfItpM2{t, 1})/2, ylim, 'lineWidth', 0.75, 'lineStyle', ':', 'alpha', 0.5);
% 
%     x1 = x1 + length(rez.dffsOnTfItpM2{t, 1}) + 1;
% end
% title("dff and lick traces");
% set(gca, 'TickDir', 'out');
% axis tight
% print(fullfile(filePath, 'Figure', 'dff_lickbouts_m2_example5trials.pdf'), '-dpdf', '-vector', '-bestfit')

%% Stim onset aligned activity (m2)
% align dffs take the mean and sem
% primary motor (stim onset)
[stimOnDff_m2_itp, stimOn_timepts] = temporalAlignInterp1(rez.stimOnDffC.m2(:, 1), rez.stimOnDffC.m2(:, 2), 0.001);
[rezM2.meanStimOnDff, ~, rezM2.semStimOnDff] = meanstdsem(cell2mat(stimOnDff_m2_itp));
[stimOnDff_m2_go_itp, ~] = temporalAlignInterp1(rez.stimOnDffC.m2(goI, 1), rez.stimOnDffC.m2(goI, 2), 0.001);
[stimOnDff_m2_nogo_itp, ~] = temporalAlignInterp1(rez.stimOnDffC.m2(nogoI, 1), rez.stimOnDffC.m2(nogoI, 2), 0.001);

[stimOnDff_m2_hit_itp, ~] = temporalAlignInterp1(rez.stimOnDffC.m2(hitI, 1), rez.stimOnDffC.m2(hitI, 2), 0.001);
[stimOnDff_m2_fa_itp, ~] = temporalAlignInterp1(rez.stimOnDffC.m2(faI, 1), rez.stimOnDffC.m2(faI, 2), 0.001);
[stimOnDff_m2_miss_itp, ~] = temporalAlignInterp1(rez.stimOnDffC.m2(missI, 1), rez.stimOnDffC.m2(missI, 2), 0.001);
[stimOnDff_m2_cr_itp, ~] = temporalAlignInterp1(rez.stimOnDffC.m2(crI, 1), rez.stimOnDffC.m2(crI, 2), 0.001);

% primary motor (stim onset) go no-go
[rezM2.meanStimOnDffGng, rezM2.semStimOnDffGng] = trialGroupMeanSem(stimOnDff_m2_itp, {[tbytDat.rewardTrI], [tbytDat.punishTrI]});
stimOnFig = plotMeanSem(rezM2.meanStimOnDffGng, rezM2.semStimOnDffGng, stimOn_timepts, {'Go', 'NoGo'});
title("Go Nogo m2 (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', strcat('dff_cueOn_Gng_m2_', channel, '.pdf')), '-dpdf', '-vector', '-bestfit')
close(stimOnFig)
fprintf("Saved tone(stim)-aligned dff plot!\n")

% primary motor (stim onset) hit vs cr
[rezM2.meanStimOnDffHitCr, rezM2.semStimOnDffHitCr] = trialGroupMeanSem(stimOnDff_m2_itp, {hitI, crI});
stimOnHitCrFig = plotMeanSem(rezM2.meanStimOnDffHitCr, rezM2.semStimOnDffHitCr, stimOn_timepts, {'Hit', 'CR'});
title("Hit vs Cr m2 (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', strcat('dff_cueOn_HitCr_m2_', channel, '.pdf')), '-dpdf', '-vector', '-bestfit')
close(stimOnHitCrFig)
fprintf("Saved tone(stim)-aligned dff plot comparing Hit vs Cr!\n")

% primary motor (stim onset) hit vs FA
[rezM2.meanStimOnDffHitFA, rezM2.semStimOnDffHitFA] = trialGroupMeanSem(stimOnDff_m2_itp, {hitI, faI});
stimOnHitFaFig = plotMeanSem(rezM2.meanStimOnDffHitFA, rezM2.semStimOnDffHitFA, stimOn_timepts, {'Hit', 'FA'});
title("Hit vs FA m2 (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', strcat('dff_cueOn_HitFA_m2_', channel, '.pdf')), '-dpdf', '-vector', '-bestfit')
close(stimOnHitFaFig)
fprintf("Saved tone(stim)-aligned dff plot comparing Hit vs FA!\n")

% primary motor (stim onset) Miss vs CR
[rezM2.meanStimOnDffMissCR, rezM2.semStimOnDffMissCR] = trialGroupMeanSem(stimOnDff_m2_itp, {missI, crI});
stimOnMissCrFig = plotMeanSem(rezM2.meanStimOnDffMissCR, rezM2.semStimOnDffMissCR, stimOn_timepts, {'Miss', 'CR'});
title("Miss vs CR m2 (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', strcat('dff_cueOn_MissCR_m2_', channel, '.pdf')), '-dpdf', '-vector', '-bestfit')
close(stimOnMissCrFig)
fprintf("Saved tone(stim)-aligned dff plot comparing Miss vs CR!\n")

% primary motor (stim onset) FA vs CR
[rezM2.meanStimOnDffFaCR, rezM2.semStimOnDffFaCR] = trialGroupMeanSem(stimOnDff_m2_itp, {faI, crI});
stimOnFaCrFig = plotMeanSem(rezM2.meanStimOnDffFaCR, rezM2.semStimOnDffFaCR, stimOn_timepts, {'FA', 'CR'});
title("FA vs CR m2 (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', strcat('dff_cueOn_FaCR_m2_', channel, '.pdf')), '-dpdf', '-vector', '-bestfit')
close(stimOnFaCrFig)
fprintf("Saved tone(stim)-aligned dff plot comparing Miss vs CR!\n")

% m2 (stim onset) go vs no-go (train classifier)
[rezM2.stimOnGoNogoSvm, rezM2.stimOnGoNogoSvmTs, rezM2.stimOnGoNogoNb, rezM2.stimOnGoNogoNbTs] = trainDffClassifier(stimOnDff_m2_go_itp, stimOnDff_m2_nogo_itp, stimOn_timepts, 50, 50, 10); 
% figure; plot(rezM2.stimOnGoNogoSvmTs, mean(rezM2.stimOnGoNogoSvm))

% m2 (stim onset) hit vs CR (train classifier)
[rezM2.stimOnHitCrSvm, rezM2.stimOnHitCrSvmTs, rezM2.stimOnHitCrNb, rezM2.stimOnHitCrNbTs] = trainDffClassifier(stimOnDff_m2_hit_itp, stimOnDff_m2_cr_itp, stimOn_timepts, 50, 50, 10); 
% figure; plot(rezM2.stimOnHitCrSvmTs, mean(rezM2.stimOnHitCrSvm))

% m2 (stim onset) hit vs FA (train classifier)
[rezM2.stimOnHitFaSvm, rezM2.stimOnHitFaSvmTs, rezM2.stimOnHitFaNb, rezM2.stimOnHitFaNbTs] = trainDffClassifier(stimOnDff_m2_hit_itp, stimOnDff_m2_fa_itp, stimOn_timepts, 50, 50, 10); 
% figure; plot(rezM2.stimOnHitFaSvmTs, mean(rezM2.stimOnHitFaSvm))

% m2 (stim onset) hit vs Miss (train classifier)
[rezM2.stimOnHitMissSvm, rezM2.stimOnHitMissSvmTs, rezM2.stimOnHitMissNb, rezM2.stimOnHitMissNbTs] = trainDffClassifier(stimOnDff_m2_hit_itp, stimOnDff_m2_miss_itp, stimOn_timepts, 50, 50, 10); 
% figure; plot(rezM2.stimOnHitMissSvmTs, mean(rezM2.stimOnHitMissSvm))

% m2 (stim onset) Miss vs CR (train classifier)
[rezM2.stimOnMissCrSvm, rezM2.stimOnMissCrSvmTs, rezM2.stimOnMissCrNb, rezM2.stimOnMissCrNbTs] = trainDffClassifier(stimOnDff_m2_miss_itp, stimOnDff_m2_cr_itp, stimOn_timepts, 50, 50, 10); 
% figure; plot(rezM2.stimOnMissCrSvmTs, mean(rezM2.stimOnMissCrSvm))

fprintf("Completed cueOn-aligned dff and classification of trial types!\n")

%% iti licks (m2) Technically there's no ITI licks with the retractable spout
% itiLickDffm2 = flatCell(rez.itiLickDff.m2); % unnest the cell arrays
% itiLickDffm2Ts = cellfun(@(a) linspace(-1, 1, length(a)), itiLickDffm2(:, 2), 'UniformOutput', false);
% [itiLickDffm2Itp, itiLickDffm2ItpTs] = temporalAlignInterp1(itiLickDffm2 (:, 1), itiLickDffm2Ts);
% [rez.meanItiLickDff.m2, ~, rez.semItiLickDff.m2] = meanstdsem(cell2mat(itiLickDffm2Itp));
% 
% plotMeanSemColor(rez.meanItiLickDff.m2, rez.semItiLickDff.m2, itiLickDffm2ItpTs, m2_color, {'m2'});
% title("lick during ITI");
% xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); ylim([-.3 .3])

%% postStim licks (m2)
postStimLickDffm2 = flatCell(rez.postStimLickDff.m2); % unnest the cell arrays
postStimLickDffm2Ts = cellfun(@(a) linspace(-1, 1, length(a)), postStimLickDffm2(:, 2), 'UniformOutput', false);
[postStimLickDffm2Itp, postStimLickDffm2ItpTs] = temporalAlignInterp1(postStimLickDffm2(:, 1), postStimLickDffm2Ts, 0.001);
[rezM2.meanPostStimLickDff, ~, rezM2.semPostStimLickDff] = meanstdsem(cell2mat(postStimLickDffm2Itp));

postStimLickFig = plotMeanSemColor(rezM2.meanPostStimLickDff, rezM2.semPostStimLickDff, postStimLickDffm2ItpTs, m2_color, {'m2'});
title("postStim licks");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); %ylim([-.3 .3])
print(fullfile(filePath, 'Figure', strcat('dff_postStimLick_Gng_m2_', channel, '.pdf')), '-dpdf', '-vector', '-bestfit')
close(postStimLickFig)
fprintf("Saved postStimLick-aligned dff plot!\n")

%% hit firstlicks (m2)
hitFstLickDffm2 = cellWithNonEmptyColumns(rez.hitDffFirstLick.m2);
hitFstLickDffm2Ts = cellfun(@(a) linspace(-1, 1, length(a)), hitFstLickDffm2(:, 2), 'UniformOutput', false);
[hitFstLickDffm2Itp, hitFstLickDffm2ItpTs] = temporalAlignInterp1(hitFstLickDffm2(:, 1), hitFstLickDffm2Ts, 0.001);
[rezM2.meanHitFstLickDff, ~, rezM2.semHitFstLickDff] = meanstdsem(cell2mat(hitFstLickDffm2Itp));

hitFstLickFig = plotMeanSemColor(rezM2.meanHitFstLickDff, rezM2.semHitFstLickDff, hitFstLickDffm2ItpTs, m2_color, {'m2'});
title("Lick (first) Hit trial");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', strcat('dff_hitFirstLick_Gng_m2_', channel, '.pdf')), '-dpdf', '-vector', '-bestfit')
close(hitFstLickFig)
fprintf("Saved hitFirstLick-aligned dff plot!\n")

%% first water (m2)
fstWaterDffm2 = cellWithNonEmptyColumns(rez.firstWater.m2);
fstWaterDffm2Ts = cellfun(@(a) linspace(-1, 1, length(a)), fstWaterDffm2(:, 2), 'UniformOutput', false);
[fstWaterDffm2Itp, fstWaterDffm2ItpTs] = temporalAlignInterp1(fstWaterDffm2 (:, 1), fstWaterDffm2Ts, 0.001);
[rezM2.meanFstWaterDff, ~, rezM2.semFstWaterDff] = meanstdsem(cell2mat(fstWaterDffm2Itp));

firstWaterFig = plotMeanSemColor(rezM2.meanFstWaterDff, rezM2.semFstWaterDff, fstWaterDffm2ItpTs, m2_color, {'m2'});
title("Water (first)");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); grid on;
print(fullfile(filePath, 'Figure', strcat('dff_firstWater_m2_', channel, '.pdf')), '-dpdf', '-vector', '-bestfit')
close(firstWaterFig)
fprintf("Saved firstWater-aligned dff plot!\n")

%% first air (m2)
if isfield(rez, 'firstAirpuff')
    fstAirDffm2 = cellWithNonEmptyColumns(rez.firstAirpuff.m2);
    fstAirDffm2Ts = cellfun(@(a) linspace(-1, 1, length(a)), fstAirDffm2(:, 2), 'UniformOutput', false);
    [fstAirDffm2Itp, fstAirDffm2ItpTs] = temporalAlignInterp1(fstAirDffm2 (:, 1), fstAirDffm2Ts, 0.001);
    [rezM2.meanFstAirDff, ~, rezM2.semFstAirDff] = meanstdsem(cell2mat(fstAirDffm2Itp));

    fstAirFig = plotMeanSemColor(rezM2.meanFstAirDff, rezM2.semFstAirDff, fstAirDffm2ItpTs, m2_color, {'m2'});
    title("Air (first)");
    xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); %ylim([-0.8 0.8])

    print(fullfile(filePath, 'Figure', strcat('dff_firstAirpuff_m2_', channel, '.pdf')), '-dpdf', '-vector', '-bestfit')
    close(fstAirFig)
    fprintf("Saved firstAir-aligned dff plot!\n")
end

%% FA first licks (m2)
faFstLickDffm2 = cellWithNonEmptyColumns(rez.faDffFirstLick.m2);
faFstLickDffm2Ts = cellfun(@(a) linspace(-1, 1, length(a)), faFstLickDffm2(:, 2), 'UniformOutput', false);
[faFstLickDffm2Itp, faFstLickDffm2ItpTs] = temporalAlignInterp1(faFstLickDffm2 (:, 1), faFstLickDffm2Ts, 0.001);
[rezM2.meanFaFstLickDff, ~, rezM2.semFaFstLickDff] = meanstdsem(cell2mat(faFstLickDffm2Itp));

faFstLickFig = plotMeanSem(rezM2.meanFaFstLickDff, rezM2.semFaFstLickDff, faFstLickDffm2ItpTs, {'m2'});
title("Lick (first) Fa trial");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); 
print(fullfile(filePath, 'Figure', strcat('dff_faFstLickFig_m2_', channel, '.pdf')), '-dpdf', '-vector', '-bestfit')
close(faFstLickFig)

%% Compare licks during Hit vs FA vs postStim (m2)
% Without normalization
hitFaPostStimFig = plotMeanSem([rezM2.meanHitFstLickDff; rezM2.meanFaFstLickDff; rezM2.meanPostStimLickDff], ...
    [rezM2.semHitFstLickDff; rezM2.semFaFstLickDff; rezM2.semPostStimLickDff], ...
    hitFstLickDffm2ItpTs, {'Hit', 'FA', 'postStim'}); 
title("Lick Hit vs FA vs postStim")
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); grid on
print(fullfile(filePath, 'Figure', strcat('dff_lick_Hit_Fa_postStim_m2_', channel, '.pdf')), '-dpdf', '-vector', '-bestfit')
close(hitFaPostStimFig)

% With normalization
rezM2.meanHitFstLickDffBaseSub = baseSubNorm(rezM2.meanHitFstLickDff, hitFstLickDffm2ItpTs, [-1 -0.5]);
rezM2.meanFaFstLickDffBaseSub = baseSubNorm(rezM2.meanFaFstLickDff, faFstLickDffm2ItpTs, [-1 -0.5]);
rezM2.meanPostStimLickDffBaseSub = baseSubNorm(rezM2.meanPostStimLickDff, postStimLickDffm2ItpTs, [-1 -0.5]);

hitFaPostStimFig_baseSub = plotMeanSem([rezM2.meanHitFstLickDffBaseSub; rezM2.meanFaFstLickDffBaseSub; rezM2.meanPostStimLickDffBaseSub], ...
    [rezM2.semHitFstLickDff; rezM2.semFaFstLickDff; rezM2.semPostStimLickDff], ...
    hitFstLickDffm2ItpTs, {'Hit', 'FA', 'postStim'}); 
title("Lick Hit vs FA vs ITI vs postStim (baseSub)")
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); grid on
print(fullfile(filePath, 'Figure', strcat('dff_lick_Hit_Fa_postStim_baseSub_m2_', channel, '.pdf')), '-dpdf', '-vector', '-bestfit')
close(hitFaPostStimFig_baseSub)

%% train svm to classify licks at different task epochs
% m2 (stim onset) go vs no-go (train classifier)
[rezM2.lickHitFaSvm, rezM2.lickHitFaSvmTs, rezM2.lickHitFaNb, rezM2.lickHitFaNbTs, X_hit, X_fa] = trainDffClassifier(hitFstLickDffm2Itp, faFstLickDffm2Itp, hitFstLickDffm2ItpTs, 50, 50, 10); 
% figure; plot(rezM2.lickHitFaSvmTs, mean(rezM2.lickHitFaSvm))


% imagesc the hit and fa trials dff
X_hit_fig = imagescWithTimeInfo(smooth2a(X_hit, 0, 1), rezM2.lickHitFaSvmTs); 
set(gca, 'XTick', -2:.5:4, 'TickDir', 'out');
clim([-0.5 0.5])
print(fullfile(filePath, 'Figure', strcat('dff_hit_timeBin_by_trial_m2_', channel', '.pdf')), '-dpdf', '-vector', '-bestfit')
close(X_hit_fig)

X_fa_fig = imagescWithTimeInfo(smooth2a(X_fa, 0, 1), rezM2.lickHitFaSvmTs); 
set(gca, 'XTick', -2:.5:4, 'TickDir', 'out');
clim([-0.7 0.7])
print(fullfile(filePath, 'Figure', strcat('dff_fa_timeBin_by_trial_m2_', channel, '.pdf')), '-dpdf', '-vector', '-bestfit')
close(X_fa_fig)

%% Miss CueOn (m2)
% primary motor (stim onset) go no-go
missCueI = [tbytDat.rewardTrI] & ~waterI & ~lickI;  
[rezM2.meanStimOnMissDff, rezM2.semStimOnMissDff] = trialGroupMeanSem(stimOnDff_m2_itp, {missCueI});
stimOnMissFig = plotMeanSem(rezM2.meanStimOnMissDff, rezM2.semStimOnMissDff, stimOn_timepts, {'m2 miss'});
title("Go Miss m2 (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', strcat('dff_cueOn_goMiss_m2_', channel, '.pdf')), '-dpdf', '-vector', '-bestfit')
close(stimOnMissFig)
fprintf("Saved tone(stim)-aligned dff of Miss trials!\n")

%% CR CueOn (m2)
% primary motor (stim onset) go no-go
crCueI = [tbytDat.punishTrI] & ~airpuffI & ~lickI;  
[rezM2.meanStimOnCrDff, rezM2.semStimOnCrDff] = trialGroupMeanSem(stimOnDff_m2_itp, {crCueI});
stimOnCrFig = plotMeanSem(rezM2.meanStimOnCrDff, rezM2.semStimOnCrDff, stimOn_timepts, {'m2 correct rejection'});
title("NoGo CR m2 (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', strcat('dff_cueOn_nogoCR_m2_', channel, '.pdf'), '-dpdf', '-vector', '-bestfit')
close(stimOnCrFig)
fprintf("Saved tone(stim)-aligned dff of Correct Rejection trials!\n")

% %% plot stim onset PETHs together
% plotMeanSem([rez.mStimOnDffm2; rez.mStimOnDffM2; rez.mStimOnDffS1; rez.mStimOnDffRs; rez.mStimOnDffV], ...
%     [rez.semStimOnDffm2; rez.semStimOnDffM2; rez.semStimOnDffS1; rez.semStimOnDffRs; rez.semStimOnDffV], ...
%     stimOn_timepts, ...
%     {'m2', 'M2', 'S1', 'M2', 'V'});
% xlabel('Time (s)')
% ylabel('DFF')
% set(gca, 'XTick', -2:1:4)
% xlim([-2 2])
% print(fullfile(filePath, 'Figure', 'dff_cueOn_Collect.pdf'), '-dpdf', '-vector', '-bestfit')

%% plot classification accuracy
% plotMeanSem(smooth2a([rez.logitGngCueOnm2; rez.logitGngCueOnM2; rez.logitGngCueOnS1; rez.logitGngCueOnRs; rez.logitGngCueOnV], 0, 2), ...
%     zeros(5, length(stimOn_timepts)), stimOn_timepts, {'m2', 'M2', 'S1', 'M2', 'V'});
% title("Cross-validated classification accuracy")
% xlabel('Time (s)')
% ylabel('Classification Accuracy')
% set(gca, 'XTick', -2:1:4, 'YTick', -1:0.1:1, 'TickDir', 'out')
% xlim([-2 2])
% ylim([0.45 0.85])
% print(fullfile(filePath, 'Figure', 'logisticRegression_Gng_cueOn_Collect.pdf'), '-dpdf', '-vector', '-bestfit')

% baseline dff image
% rez.meanBaseDffImgMotor = nanmean(cell2mat(reshape(rwdDffC_motor(:, 3), [1, 1, length(tbytDat)])), 3);
% imageFrameWithNaNsEdgeMap(rez.meanBaseDffImgMotor, [-.3 .3], dom2alMaps.edgeMap, [0, 0, 1], 0.3)
% title('mean baseline m2')
% set(gca, 'XTickLabel', {})
% set(gca, 'YTickLabel', {})
% print(fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA003/DA003_083023/Figure', 'dff_base_m2_Img.pdf'), '-dpdf', '-vector', '-bestfit')
% 
% % reward dff image
% rez.meanRwdDffImgMotor = nanmean(cell2mat(reshape(rwdDffC_motor(:, 4), [1, 1, length(tbytDat)])), 3);
% imageFrameWithNaNsEdgeMap(rez.meanRwdDffImgMotor, [-.3 .3], dom2alMaps.edgeMap, [0, 0, 1], 0.3)
% print(fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA003/DA003_083023/Figure', 'dff_reward_m2_Img.pdf'), '-dpdf', '-vector', '-bestfit')
% title('mean reward m2')
% set(gca, 'XTickLabel', {})
% set(gca, 'YTickLabel', {})
% print(fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA003/DA003_083023/Figure', 'dff_reward_m2_Img.pdf'), '-dpdf', '-vector', '-bestfit')

%showAllenEdgesOnTransformedWF(rez.meanBaseDffImgMotor, dom2alMaps, 0, 1)

%imageFrameWithNaNs(rez.meanRwdDffImgMotor, [-.3 .3])


%% save rez
save(fullfile(filePath, 'Matfiles', strcat(header, '_dff_evtAligned_regionMask_m2_', channel, '.mat')), 'rezM2');

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

        % tbytPsthC = stimOnDff_m2_itp;
        % groupC = {[tbytDat.rewardTrI], [tbytDat.punishTrI]};

        % sanity check 1: all trial numbem2 must match!
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
        % %   e.g., cell2mat(stimOnDff_m2_itp);
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

            for m2 = 1:resample
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
                accuracy(m2, jj) = sum(YPred == YTest) / length(YTest);
                fprintf('Finished resample #%d of bin #%d with accuracy: %.2f%%\n', m2, jj, accuracy(m2, jj) * 100);
            end
        end

    end



end


