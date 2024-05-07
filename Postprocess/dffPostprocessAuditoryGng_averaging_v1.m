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

%% cross correlogram lick bouts and dff
xCorrLagMat = cell2mat(rez.xcorrDffLick.V1(:, 2)); 
xCorrLag = xCorrLagMat(1, :); 

[rezV1.meanXcorrDffLick, ~, rezV1.semXcorrDffLick] = meanstdsem(cell2mat(rez.xcorrDffLick.V1(:, 1)));
xcorrFig = plotMeanSemColor(rezV1.meanXcorrDffLick, rezV1.semXcorrDffLick, xCorrLag, v1_color, {'xcorr v1'});
xlabel('Time (s)'); ylabel('xcorr'); set(gca, 'TickDir', 'out', 'XTick', -2000:1000:2000); grid on
print(fullfile(filePath, 'Figure', 'dff_cueOn_Gng_v1.pdf'), '-dpdf', '-vector', '-bestfit')
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
%print(fullfile(filePath, 'Figure', 'dff_lickbouts_v1_example5trials.pdf'), '-dpdf', '-vector', '-bestfit')

%% Stim onset aligned activity (V1)
% align dffs take the mean and sem
% primary motor (stim onset)
[stimOnDff_v1_itp, stimOn_timepts] = temporalAlignInterp1(rez.stimOnDffC.v1(:, 1), rez.stimOnDffC.v1(:, 2), 0.001);
[rezV1.meanStimOnDff, ~, rezV1.semStimOnDff] = meanstdsem(cell2mat(stimOnDff_v1_itp));
[stimOnDff_v1_go_itp, stimOn_timepts_go] = temporalAlignInterp1(rez.stimOnDffC.v1(goI, 1), rez.stimOnDffC.v1(goI, 2), 0.001);
[stimOnDff_v1_nogo_itp, stimOn_timepts_nogo] = temporalAlignInterp1(rez.stimOnDffC.v1(nogoI, 1), rez.stimOnDffC.v1(nogoI, 2), 0.001);

% primary motor (stim onset) go no-go
[rezV1.meanStimOnDffGng, rezV1.semStimOnDffGng] = trialGroupMeanSem(stimOnDff_v1_itp, {[tbytDat.rewardTrI], [tbytDat.punishTrI]});
stimOnFig = plotMeanSem(rezV1.meanStimOnDffGng, rezV1.semStimOnDffGng, stimOn_timepts, {'Go', 'NoGo'});
title("Go Nogo V1 (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_cueOn_Gng_v1.pdf'), '-dpdf', '-vector', '-bestfit')
close(stimOnFig)
fprintf("Saved tone(stim)-aligned dff plot!\n")

% primary motor (stim onset) go no-go (logistic regression)
[X_stimOnGo, ~] = binAvg1msSpkCountMat(cell2mat(stimOnDff_v1_go_itp), 50, 50); 
[X_stimOnNoGo, ~] = binAvg1msSpkCountMat(cell2mat(stimOnDff_v1_nogo_itp), 50, 50); 

X_stimOn_bins = stimOn_timepts(1:50:end); 

if length(X_stimOn_bins)-size(X_stimOnGo, 2)==1
    X_stimOn_bins = X_stimOn_bins(1:end-1)+0.025; 
end

assert(size(X_stimOnGo, 2)==size(X_stimOnNoGo, 2))

% train and test svm
Xs = [X_stimOnGo; X_stimOnNoGo];
y = [ones(size(X_stimOnGo, 1), 1)*1; ones(size(X_stimOnNoGo, 1), 1)*2];
rezV1.stimOnGoNogoSvm = multiClass_svm_peth(Xs, y, 10);
rezV1.stimOnGoNogoSvmTs = X_stimOn_bins; 
rezV1.stimOnGoNogoNb = multiClass_naiveBayesClassifier_peth(Xs, y, 10);
rezV1.stimOnGoNogoNbTs = X_stimOn_bins; 

fprintf("Completed cueOn-aligned dff and classification of go vs nogo trials!\n")

%% iti licks (V1) Technically there's no ITI licks with the retractable spout
% itiLickDffV1 = flatCell(rez.itiLickDff.v1); % unnest the cell arrays
% itiLickDffV1Ts = cellfun(@(a) linspace(-1, 1, length(a)), itiLickDffV1(:, 2), 'UniformOutput', false);
% [itiLickDffV1Itp, itiLickDffV1ItpTs] = temporalAlignInterp1(itiLickDffV1 (:, 1), itiLickDffV1Ts);
% [rez.meanItiLickDff.v1, ~, rez.semItiLickDff.v1] = meanstdsem(cell2mat(itiLickDffV1Itp));
% 
% plotMeanSemColor(rez.meanItiLickDff.v1, rez.semItiLickDff.v1, itiLickDffV1ItpTs, v1_color, {'V1'});
% title("lick during ITI");
% xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); ylim([-.3 .3])

%% postStim licks (V1)
postStimLickDffV1 = flatCell(rez.postStimLickDff.v1); % unnest the cell arrays
postStimLickDffV1Ts = cellfun(@(a) linspace(-1, 1, length(a)), postStimLickDffV1(:, 2), 'UniformOutput', false);
[postStimLickDffV1Itp, postStimLickDffV1ItpTs] = temporalAlignInterp1(postStimLickDffV1(:, 1), postStimLickDffV1Ts, 0.001);
[rezV1.meanPostStimLickDff, ~, rezV1.semPostStimLickDff] = meanstdsem(cell2mat(postStimLickDffV1Itp));

postStimLickFig = plotMeanSemColor(rezV1.meanPostStimLickDff, rezV1.semPostStimLickDff, postStimLickDffV1ItpTs, v1_color, {'V1'});
title("lick during ITI");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); %ylim([-.3 .3])
print(fullfile(filePath, 'Figure', 'dff_postStimLick_Gng_v1.pdf'), '-dpdf', '-vector', '-bestfit')
close(postStimLickFig)
fprintf("Saved postStimLick-aligned dff plot!\n")

%% hit first licks (V1)
hitFstLickDffV1 = cellWithNonEmptyColumns(rez.hitDffFirstLick.v1);
hitFstLickDffV1Ts = cellfun(@(a) linspace(-1, 1, length(a)), hitFstLickDffV1(:, 2), 'UniformOutput', false);
[hitFstLickDffV1Itp, hitFstLickDffV1ItpTs] = temporalAlignInterp1(hitFstLickDffV1(:, 1), hitFstLickDffV1Ts, 0.001);
[rezV1.meanHitFstLickDff, ~, rezV1.semHitFstLickDff] = meanstdsem(cell2mat(hitFstLickDffV1Itp));

hitFstLickFig = plotMeanSemColor(rezV1.meanHitFstLickDff, rezV1.semHitFstLickDff, hitFstLickDffV1ItpTs, v1_color, {'V1'});
title("Lick (first) Hit trial");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_hitFirstLick_Gng_v1.pdf'), '-dpdf', '-vector', '-bestfit')
close(hitFstLickFig)
fprintf("Saved hitFirstLick-aligned dff plot!\n")

%% first water (V1)
fstWaterDffV1 = cellWithNonEmptyColumns(rez.firstWater.v1);
fstWaterDffV1Ts = cellfun(@(a) linspace(-1, 1, length(a)), fstWaterDffV1(:, 2), 'UniformOutput', false);
[fstWaterDffV1Itp, fstWaterDffV1ItpTs] = temporalAlignInterp1(fstWaterDffV1 (:, 1), fstWaterDffV1Ts, 0.001);
[rezV1.meanFstWaterDff, ~, rezV1.semFstWaterDff] = meanstdsem(cell2mat(fstWaterDffV1Itp));

firstWaterFig = plotMeanSemColor(rezV1.meanFstWaterDff, rezV1.semFstWaterDff, fstWaterDffV1ItpTs, v1_color, {'V1'});
title("Water (first)");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); grid on;
print(fullfile(filePath, 'Figure', 'dff_firstWater_v1.pdf'), '-dpdf', '-vector', '-bestfit')
close(firstWaterFig)
fprintf("Saved firstWater-aligned dff plot!\n")

%% first air (V1)
if isfield(rez, 'firstAirpuff')
    fstAirDffV1 = cellWithNonEmptyColumns(rez.firstAirpuff.v1);
    fstAirDffV1Ts = cellfun(@(a) linspace(-1, 1, length(a)), fstAirDffV1(:, 2), 'UniformOutput', false);
    [fstAirDffV1Itp, fstAirDffV1ItpTs] = temporalAlignInterp1(fstAirDffV1 (:, 1), fstAirDffV1Ts);
    [rezV1.meanFstAirDff, ~, rezV1.semFstAirDff] = meanstdsem(cell2mat(fstAirDffV1Itp));

    fstAirFig = plotMeanSemColor(rezV1.meanFstAirDff, rezV1.semFstAirDff, fstAirDffV1ItpTs, v1_color, {'V1'});
    title("Air (first)");
    xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); ylim([-0.8 0.8])

    print(fullfile(filePath, 'Figure', 'dff_firstAirpuff_v1.pdf'), '-dpdf', '-vector', '-bestfit')
    close(fstAirFig)
    fprintf("Saved firstAir-aligned dff plot!\n")
end

%% FA first licks (V1)
faFstLickDffV1 = cellWithNonEmptyColumns(rez.faDffFirstLick.v1);
faFstLickDffV1Ts = cellfun(@(a) linspace(-1, 1, length(a)), faFstLickDffV1(:, 2), 'UniformOutput', false);
[faFstLickDffV1Itp, faFstLickDffV1ItpTs] = temporalAlignInterp1(faFstLickDffV1 (:, 1), faFstLickDffV1Ts, 0.001);
[rezV1.meanFaFstLickDff, ~, rezV1.semFaFstLickDff] = meanstdsem(cell2mat(faFstLickDffV1Itp));

faFstLickFig = plotMeanSem(rezV1.meanFaFstLickDff, rezV1.semFaFstLickDff, faFstLickDffV1ItpTs, {'V1'});
title("Lick (first) Fa trial");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); 
print(fullfile(filePath, 'Figure', 'dff_faFstLickFig_v1.pdf'), '-dpdf', '-vector', '-bestfit')
close(faFstLickFig)

%% Compare licks during Hit vs FA vs postStim (V1)
% Without normalization
hitFaPostStimFig = plotMeanSem([rezV1.meanHitFstLickDff; rezV1.meanFaFstLickDff; rezV1.meanPostStimLickDff], ...
    [rezV1.semHitFstLickDff; rezV1.semFaFstLickDff; rezV1.semPostStimLickDff], ...
    hitFstLickDffV1ItpTs, {'Hit', 'FA', 'postStim'}); 
title("Lick Hit vs FA vs postStim")
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); grid on
print(fullfile(filePath, 'Figure', 'dff_lick_Hit_Fa_postStim_v1.pdf'), '-dpdf', '-vector', '-bestfit')
close(hitFaPostStimFig)

% With normalization
rezV1.meanHitFstLickDffBaseSub = baseSubNorm(rezV1.meanHitFstLickDff, hitFstLickDffV1ItpTs, [-1 -0.5]);
rezV1.meanFaFstLickDffBaseSub = baseSubNorm(rezV1.meanFaFstLickDff, faFstLickDffV1ItpTs, [-1 -0.5]);
rezV1.meanPostStimLickDffBaseSub = baseSubNorm(rezV1.meanPostStimLickDff, postStimLickDffV1ItpTs, [-1 -0.5]);

hitFaPostStimFig_baseSub = plotMeanSem([rezV1.meanHitFstLickDffBaseSub; rezV1.meanFaFstLickDffBaseSub; rezV1.meanPostStimLickDffBaseSub], ...
    [rezV1.semHitFstLickDff; rezV1.semFaFstLickDff; rezV1.semPostStimLickDff], ...
    hitFstLickDffV1ItpTs, {'Hit', 'FA', 'postStim'}); 
title("Lick Hit vs FA vs ITI vs postStim (baseSub)")
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); grid on
print(fullfile(filePath, 'Figure', 'dff_lick_Hit_Fa_postStim_baseSub_v1.pdf'), '-dpdf', '-vector', '-bestfit')
close(hitFaPostStimFig_baseSub)

%% train svm to classify licks at different task epochs
[X_hit, ~] = binAvg1msSpkCountMat(cell2mat(hitFstLickDffV1Itp), 50, 50); 
[X_fa, ~] = binAvg1msSpkCountMat(cell2mat(faFstLickDffV1Itp), 50, 50); 

X_hit_bins = hitFstLickDffV1ItpTs(1:50:end); 
if length(X_hit_bins)-size(X_hit, 2)==1
    X_hit_bins = X_hit_bins(1:end-1)+0.025; 
end

% imagesc the hit and fa trials dff
X_hit_fig = imagescWithTimeInfo(smooth2a(X_hit, 0, 1), X_hit_bins); 
set(gca, 'XTick', -2:.5:4, 'TickDir', 'out');
clim([-0.5 0.5])
print(fullfile(filePath, 'Figure', 'dff_hit_timeBin_by_trial_v1.pdf'), '-dpdf', '-vector', '-bestfit')
close(X_hit_fig)

X_fa_fig = imagescWithTimeInfo(smooth2a(X_fa, 0, 1), X_hit_bins); 
set(gca, 'XTick', -2:.5:4, 'TickDir', 'out');
clim([-0.5 0.5])
print(fullfile(filePath, 'Figure', 'dff_fa_timeBin_by_trial_v1.pdf'), '-dpdf', '-vector', '-bestfit')
close(X_fa_fig)

% train and test svm
Xs = [X_hit; X_fa];
y = [ones(size(X_hit, 1), 1)*1; ones(size(X_fa, 1), 1)*2];
rezV1.lickHitFaSvm = multiClass_svm_peth(Xs, y, 10);
rezV1.lickHitFaSvmTs = X_hit_bins; 
rezV1.lickHitFaNb = multiClass_naiveBayesClassifier_peth(Xs, y, 10);
rezV1.lickHitFaNbTs = X_hit_bins; 

%% Miss CueOn (V1)
% primary motor (stim onset) go no-go
missCueI = [tbytDat.rewardTrI] & ~waterI & ~lickI;  
[rezV1.meanStimOnMissDff, rezV1.semStimOnMissDff] = trialGroupMeanSem(stimOnDff_v1_itp, {missCueI});
stimOnMissFig = plotMeanSem(rezV1.meanStimOnMissDff, rezV1.semStimOnMissDff, stimOn_timepts, {'V1 miss'});
title("Go Miss V1 (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_cueOn_goMiss_v1.pdf'), '-dpdf', '-vector', '-bestfit')
close(stimOnMissFig)
fprintf("Saved tone(stim)-aligned dff of Miss trials!\n")

%% CR CueOn (V1)
% primary motor (stim onset) go no-go
crCueI = [tbytDat.punishTrI] & ~airpuffI & ~lickI;  
[rezV1.meanStimOnCrDff, rezV1.semStimOnCrDff] = trialGroupMeanSem(stimOnDff_v1_itp, {crCueI});
stimOnCrFig = plotMeanSem(rezV1.meanStimOnCrDff, rezV1.semStimOnCrDff, stimOn_timepts, {'V1 correct rejection'});
title("NoGo CR V1 (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_cueOn_nogoCR_v1.pdf'), '-dpdf', '-vector', '-bestfit')
close(stimOnCrFig)
fprintf("Saved tone(stim)-aligned dff of Correct Rejection trials!\n")

% %% plot stim onset PETHs together
% plotMeanSem([rez.mStimOnDffV1; rez.mStimOnDffM2; rez.mStimOnDffS1; rez.mStimOnDffRs; rez.mStimOnDffV], ...
%     [rez.semStimOnDffV1; rez.semStimOnDffM2; rez.semStimOnDffS1; rez.semStimOnDffRs; rez.semStimOnDffV], ...
%     stimOn_timepts, ...
%     {'V1', 'M2', 'S1', 'RS', 'V'});
% xlabel('Time (s)')
% ylabel('DFF')
% set(gca, 'XTick', -2:1:4)
% xlim([-2 2])
% print(fullfile(filePath, 'Figure', 'dff_cueOn_Collect.pdf'), '-dpdf', '-vector', '-bestfit')

%% plot classification accuracy
% plotMeanSem(smooth2a([rez.logitGngCueOnV1; rez.logitGngCueOnM2; rez.logitGngCueOnS1; rez.logitGngCueOnRs; rez.logitGngCueOnV], 0, 2), ...
%     zeros(5, length(stimOn_timepts)), stimOn_timepts, {'V1', 'M2', 'S1', 'RS', 'V'});
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
% title('mean baseline V1')
% set(gca, 'XTickLabel', {})
% set(gca, 'YTickLabel', {})
% print(fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA003/DA003_083023/Figure', 'dff_base_V1_Img.pdf'), '-dpdf', '-vector', '-bestfit')
% 
% % reward dff image
% rez.meanRwdDffImgMotor = nanmean(cell2mat(reshape(rwdDffC_motor(:, 4), [1, 1, length(tbytDat)])), 3);
% imageFrameWithNaNsEdgeMap(rez.meanRwdDffImgMotor, [-.3 .3], dorsalMaps.edgeMap, [0, 0, 1], 0.3)
% print(fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA003/DA003_083023/Figure', 'dff_reward_V1_Img.pdf'), '-dpdf', '-vector', '-bestfit')
% title('mean reward V1')
% set(gca, 'XTickLabel', {})
% set(gca, 'YTickLabel', {})
% print(fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA003/DA003_083023/Figure', 'dff_reward_V1_Img.pdf'), '-dpdf', '-vector', '-bestfit')

%showAllenEdgesOnTransformedWF(rez.meanBaseDffImgMotor, dorsalMaps, 0, 1)

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


