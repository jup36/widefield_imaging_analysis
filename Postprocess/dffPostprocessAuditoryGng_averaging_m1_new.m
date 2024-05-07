function dffPostprocessAuditoryGng_averaging_m1(filePath)
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

%% cross correlogram lick bouts and dff
xCorrLagMat = cell2mat(rez.xcorrDffLick.m1(:, 2)); 
xCorrLag = xCorrLagMat(1, :); 

[rezM1.meanXcorrDffLick, ~, rezM1.semXcorrDffLick] = meanstdsem(cell2mat(rez.xcorrDffLick.m1(:, 1)));
xcorrFig = plotMeanSemColor(rezM1.meanXcorrDffLick, rezM1.semXcorrDffLick, xCorrLag, m1_color, {'xcorr m1'});
xlabel('Time (s)'); ylabel('xcorr'); set(gca, 'TickDir', 'out', 'XTick', -2000:1000:2000); grid on
print(fullfile(filePath, 'Figure', 'dff_cueOn_Gng_m1.pdf'), '-dpdf', '-vector', '-bestfit')
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
%print(fullfile(filePath, 'Figure', 'dff_lickbouts_m1_example5trials.pdf'), '-dpdf', '-vector', '-bestfit')

%% Stim onset aligned activity (M1)
% align dffs take the mean and sem
% primary motor (stim onset)
[stimOnDff_m1_itp, timepts] = temporalAlignInterp1(rez.stimOnDffC.m1(:, 1), rez.stimOnDffC.m1(:, 2), 0.001);
[rezM1.meanStimOnDff, ~, rezM1.semStimOnDff] = meanstdsem(cell2mat(stimOnDff_m1_itp));

% primary motor (stim onset) go no-go
[rezM1.meanStimOnDffGng, rezM1.semStimOnDffGng] = trialGroupMeanSem(stimOnDff_m1_itp, {[tbytDat.rewardTrI], [tbytDat.punishTrI]});
stimOnFig = plotMeanSem(rezM1.meanStimOnDffGng, rezM1.semStimOnDffGng, timepts, {'Go', 'NoGo'});
title("Go Nogo M1 (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_cueOn_Gng_m1.pdf'), '-dpdf', '-vector', '-bestfit')
close(stimOnFig)
fprintf("Saved tone(stim)-aligned dff plot!\n")

% primary motor (stim onset) go no-go (logistic regression)
rezM1.logitGngCueOn = logireg(cell2mat(stimOnDff_m1_itp), [tbytDat.rewardTrI]); % logistic regression
fprintf("Completed logistic regression on cueOn-aligned dff!\n")

%% iti licks (M1) Technically there's no ITI licks with the retractable spout
% itiLickDffM1 = flatCell(rez.itiLickDff.m1); % unnest the cell arrays
% itiLickDffM1Ts = cellfun(@(a) linspace(-1, 1, length(a)), itiLickDffM1(:, 2), 'UniformOutput', false);
% [itiLickDffM1Itp, itiLickDffM1ItpTs] = temporalAlignInterp1(itiLickDffM1 (:, 1), itiLickDffM1Ts);
% [rez.meanItiLickDff.m1, ~, rez.semItiLickDff.m1] = meanstdsem(cell2mat(itiLickDffM1Itp));
% 
% plotMeanSemColor(rez.meanItiLickDff.m1, rez.semItiLickDff.m1, itiLickDffM1ItpTs, m1_color, {'M1'});
% title("lick during ITI");
% xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); ylim([-.3 .3])

%% postStim licks (M1)
postStimLickDffM1 = flatCell(rez.postStimLickDff.m1); % unnest the cell arrays
postStimLickDffM1Ts = cellfun(@(a) linspace(-1, 1, length(a)), postStimLickDffM1(:, 2), 'UniformOutput', false);
[postStimLickDffM1Itp, postStimLickDffM1ItpTs] = temporalAlignInterp1(postStimLickDffM1(:, 1), postStimLickDffM1Ts, 0.001);
[rezM1.meanPostStimLickDff, ~, rezM1.semPostStimLickDff] = meanstdsem(cell2mat(postStimLickDffM1Itp));

postStimLickFig = plotMeanSemColor(rezM1.meanPostStimLickDff, rezM1.semPostStimLickDff, postStimLickDffM1ItpTs, m1_color, {'M1'});
title("lick during ITI");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); %ylim([-.3 .3])
print(fullfile(filePath, 'Figure', 'dff_postStimLick_Gng_m1.pdf'), '-dpdf', '-vector', '-bestfit')
close(postStimLickFig)
fprintf("Saved postStimLick-aligned dff plot!\n")

%% hit first licks (M1)
hitFstLickDffM1 = cellWithNonEmptyColumns(rez.hitDffFirstLick.m1);
hitFstLickDffM1Ts = cellfun(@(a) linspace(-1, 1, length(a)), hitFstLickDffM1(:, 2), 'UniformOutput', false);
[hitFstLickDffM1Itp, hitFstLickDffM1ItpTs] = temporalAlignInterp1(hitFstLickDffM1(:, 1), hitFstLickDffM1Ts, 0.001);
[rezM1.meanHitFstLickDff, ~, rezM1.semHitFstLickDff] = meanstdsem(cell2mat(hitFstLickDffM1Itp));

hitFstLickFig = plotMeanSemColor(rezM1.meanHitFstLickDff, rezM1.semHitFstLickDff, hitFstLickDffM1ItpTs, m1_color, {'M1'});
title("Lick (first) Hit trial");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_hitFirstLick_Gng_m1.pdf'), '-dpdf', '-vector', '-bestfit')
close(hitFstLickFig)
fprintf("Saved hitFirstLick-aligned dff plot!\n")

%% first water (M1)
fstWaterDffM1 = cellWithNonEmptyColumns(rez.firstWater.m1);
fstWaterDffM1Ts = cellfun(@(a) linspace(-1, 1, length(a)), fstWaterDffM1(:, 2), 'UniformOutput', false);
[fstWaterDffM1Itp, fstWaterDffM1ItpTs] = temporalAlignInterp1(fstWaterDffM1 (:, 1), fstWaterDffM1Ts, 0.001);
[rezM1.meanFstWaterDff, ~, rezM1.semFstWaterDff] = meanstdsem(cell2mat(fstWaterDffM1Itp));

firstWaterFig = plotMeanSemColor(rezM1.meanFstWaterDff, rezM1.semFstWaterDff, fstWaterDffM1ItpTs, m1_color, {'M1'});
title("Water (first)");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); grid on;
print(fullfile(filePath, 'Figure', 'dff_firstWater_m1.pdf'), '-dpdf', '-vector', '-bestfit')
close(firstWaterFig)
fprintf("Saved firstWater-aligned dff plot!\n")

%% first air (M1)
if isfield(rez, 'firstAirpuff')
    fstAirDffM1 = cellWithNonEmptyColumns(rez.firstAirpuff.m1);
    fstAirDffM1Ts = cellfun(@(a) linspace(-1, 1, length(a)), fstAirDffM1(:, 2), 'UniformOutput', false);
    [fstAirDffM1Itp, fstAirDffM1ItpTs] = temporalAlignInterp1(fstAirDffM1 (:, 1), fstAirDffM1Ts);
    [rezM1.meanFstAirDff, ~, rezM1.semFstAirDff] = meanstdsem(cell2mat(fstAirDffM1Itp));

    fstAirFig = plotMeanSemColor(rezM1.meanFstAirDff, rezM1.semFstAirDff, fstAirDffM1ItpTs, m1_color, {'M1'});
    title("Air (first)");
    xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); ylim([-0.8 0.8])

    print(fullfile(filePath, 'Figure', 'dff_firstAirpuff_m1.pdf'), '-dpdf', '-vector', '-bestfit')
    close(fstAirFig)
    fprintf("Saved firstAir-aligned dff plot!\n")
end

%% FA first licks (M1)
faFstLickDffM1 = cellWithNonEmptyColumns(rez.faDffFirstLick.m1);
faFstLickDffM1Ts = cellfun(@(a) linspace(-1, 1, length(a)), faFstLickDffM1(:, 2), 'UniformOutput', false);
[faFstLickDffM1Itp, faFstLickDffM1ItpTs] = temporalAlignInterp1(faFstLickDffM1 (:, 1), faFstLickDffM1Ts, 0.001);
[rezM1.meanFaFstLickDff, ~, rezM1.semFaFstLickDff] = meanstdsem(cell2mat(faFstLickDffM1Itp));

faFstLickFig = plotMeanSem(rezM1.meanFaFstLickDff, rezM1.semFaFstLickDff, faFstLickDffM1ItpTs, {'M1'});
title("Lick (first) Fa trial");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); 
print(fullfile(filePath, 'Figure', 'dff_faFstLickFig_m1.pdf'), '-dpdf', '-vector', '-bestfit')
close(faFstLickFig)

%% Compare licks during Hit vs FA vs postStim (M1)
% Without normalization
hitFaPostStimFig = plotMeanSem([rezM1.meanHitFstLickDff; rezM1.meanFaFstLickDff; rezM1.meanPostStimLickDff], ...
    [rezM1.semHitFstLickDff; rezM1.semFaFstLickDff; rezM1.semPostStimLickDff], ...
    hitFstLickDffM1ItpTs, {'Hit', 'FA', 'postStim'}); 
title("Lick Hit vs FA vs postStim")
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); grid on
print(fullfile(filePath, 'Figure', 'dff_lick_Hit_Fa_postStim_m1.pdf'), '-dpdf', '-vector', '-bestfit')
close(hitFaPostStimFig)

% With normalization
rezM1.meanHitFstLickDffBaseSub = baseSubNorm(rezM1.meanHitFstLickDff, hitFstLickDffM1ItpTs, [-1 -0.5]);
rezM1.meanFaFstLickDffBaseSub = baseSubNorm(rezM1.meanFaFstLickDff, faFstLickDffM1ItpTs, [-1 -0.5]);
rezM1.meanPostStimLickDffBaseSub = baseSubNorm(rezM1.meanPostStimLickDff, postStimLickDffM1ItpTs, [-1 -0.5]);

hitFaPostStimFig_baseSub = plotMeanSem([rezM1.meanHitFstLickDffBaseSub; rezM1.meanFaFstLickDffBaseSub; rezM1.meanPostStimLickDffBaseSub], ...
    [rezM1.semHitFstLickDff; rezM1.semFaFstLickDff; rezM1.semPostStimLickDff], ...
    hitFstLickDffM1ItpTs, {'Hit', 'FA', 'postStim'}); 
title("Lick Hit vs FA vs ITI vs postStim (baseSub)")
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); grid on
print(fullfile(filePath, 'Figure', 'dff_lick_Hit_Fa_postStim_baseSub_m1.pdf'), '-dpdf', '-vector', '-bestfit')
close(hitFaPostStimFig_baseSub)

%% train svm to classify licks at different task epochs
[X_hit, X_hit_ts] = binAvg1msSpkCountMat(cell2mat(hitFstLickDffM1Itp), 50, 50); 
[X_fa, X_fa_ts] = binAvg1msSpkCountMat(cell2mat(faFstLickDffM1Itp), 50, 50); 

X_hit_bins = hitFstLickDffM1ItpTs(1:50:end); 
if length(X_hit_bins)-size(X_hit, 2)==1
    X_hit_bins = X_hit_bins(1:end-1)+0.025; 
end

% imagesc the hit and fa trials dff
X_hit_fig = imagescWithTimeInfo(smooth2a(X_hit, 0, 1), X_hit_bins); 
set(gca, 'XTick', -2:.5:4, 'TickDir', 'out');
caxis([-0.5 0.5])
print(fullfile(filePath, 'Figure', 'dff_hit_timeBin_by_trial_m1.pdf'), '-dpdf', '-vector', '-bestfit')
close(X_hit_fig)

X_fa_fig = imagescWithTimeInfo(smooth2a(X_fa, 0, 1), X_hit_bins); 
set(gca, 'XTick', -2:.5:4, 'TickDir', 'out');
caxis([-0.5 0.5])
print(fullfile(filePath, 'Figure', 'dff_fa_timeBin_by_trial_m1.pdf'), '-dpdf', '-vector', '-bestfit')
close(X_fa_fig)

Xs = [X_hit; X_fa];
y = [ones(size(X_hit, 1), 1)*1; ones(size(X_fa, 1), 1)*2];
rezM1.lickHitFaSvm = multiClass_svm_peth(Xs, y, 10);
rezM1.lickHitFaSvmTs = X_hit_bins; 
rezM1.lickHitFaNb = multiClass_naiveBayesClassifier_peth(Xs, y, 10);
rezM1.lickHitFaNbTs = X_hit_bins; 

%% Miss CueOn (M1)
% primary motor (stim onset) go no-go
missCueI = [tbytDat.rewardTrI] & ~waterI & ~lickI;  
[rezM1.meanStimOnMissDff, rezM1.semStimOnMissDff] = trialGroupMeanSem(stimOnDff_m1_itp, {missCueI});
stimOnMissFig = plotMeanSem(rezM1.meanStimOnMissDff, rezM1.semStimOnMissDff, timepts, {'M1 miss'});
title("Go Miss M1 (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_cueOn_goMiss_m1.pdf'), '-dpdf', '-vector', '-bestfit')
close(stimOnMissFig)
fprintf("Saved tone(stim)-aligned dff of Miss trials!\n")

%% CR CueOn (M1)
% primary motor (stim onset) go no-go
crCueI = [tbytDat.punishTrI] & ~airpuffI & ~lickI;  
[rezM1.meanStimOnCrDff, rezM1.semStimOnCrDff] = trialGroupMeanSem(stimOnDff_m1_itp, {crCueI});
stimOnCrFig = plotMeanSem(rezM1.meanStimOnCrDff, rezM1.semStimOnCrDff, timepts, {'M1 correct rejection'});
title("NoGo CR M1 (cue onset at time=0)");
set(gca, 'TickDir', 'out')
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:5, 'TickDir', 'out', 'YTick', -1:0.1:5); xlim([-1 5]); %ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_cueOn_nogoCR_m1.pdf'), '-dpdf', '-vector', '-bestfit')
close(stimOnCrFig)
fprintf("Saved tone(stim)-aligned dff of Correct Rejection trials!\n")

% %% plot stim onset PETHs together
% plotMeanSem([rez.mStimOnDffM1; rez.mStimOnDffM2; rez.mStimOnDffS1; rez.mStimOnDffRs; rez.mStimOnDffV], ...
%     [rez.semStimOnDffM1; rez.semStimOnDffM2; rez.semStimOnDffS1; rez.semStimOnDffRs; rez.semStimOnDffV], ...
%     timepts, ...
%     {'M1', 'M2', 'S1', 'RS', 'V'});
% xlabel('Time (s)')
% ylabel('DFF')
% set(gca, 'XTick', -2:1:4)
% xlim([-2 2])
% print(fullfile(filePath, 'Figure', 'dff_cueOn_Collect.pdf'), '-dpdf', '-vector', '-bestfit')

%% plot classification accuracy
% plotMeanSem(smooth2a([rez.logitGngCueOnM1; rez.logitGngCueOnM2; rez.logitGngCueOnS1; rez.logitGngCueOnRs; rez.logitGngCueOnV], 0, 2), ...
%     zeros(5, length(timepts)), timepts, {'M1', 'M2', 'S1', 'RS', 'V'});
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
% title('mean baseline M1')
% set(gca, 'XTickLabel', {})
% set(gca, 'YTickLabel', {})
% print(fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA003/DA003_083023/Figure', 'dff_base_M1_Img.pdf'), '-dpdf', '-vector', '-bestfit')
% 
% % reward dff image
% rez.meanRwdDffImgMotor = nanmean(cell2mat(reshape(rwdDffC_motor(:, 4), [1, 1, length(tbytDat)])), 3);
% imageFrameWithNaNsEdgeMap(rez.meanRwdDffImgMotor, [-.3 .3], dorsalMaps.edgeMap, [0, 0, 1], 0.3)
% print(fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA003/DA003_083023/Figure', 'dff_reward_M1_Img.pdf'), '-dpdf', '-vector', '-bestfit')
% title('mean reward M1')
% set(gca, 'XTickLabel', {})
% set(gca, 'YTickLabel', {})
% print(fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA003/DA003_083023/Figure', 'dff_reward_M1_Img.pdf'), '-dpdf', '-vector', '-bestfit')

%showAllenEdgesOnTransformedWF(rez.meanBaseDffImgMotor, dorsalMaps, 0, 1)

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


