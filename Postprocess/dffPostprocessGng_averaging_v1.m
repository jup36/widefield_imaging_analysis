function dffPostprocev1Gng_averaging_v1(filePath)
%filePath = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA008/DA008_101823';
% This function performs further analyses and plotting using data in
% 'rez'(outcome from 'dffPostprocev1Gng_averaging.m'). 
% 

%% load rez tbytDat
[~, header] = fileparts(filePath);
load(fullfile(filePath, 'Matfiles', strcat(header, '_dff_evtAligned_regionMask.mat')), 'rez'); 
load(fullfile(filePath, 'Matfiles', strcat(header, '_tbytDat_pav1eGng')), 'tbytDat')

% user variables 
v1_color = [255 0 255]./255; 

%% crov1 correlogram lick bouts and dff
[rez.meanXcorrDffLick.v1, ~, rez.semXcorrDffLick.v1] = meanstdsem(cell2mat(rez.xcorrDffLick.V1(:, 1))); 
plotMeanSemColor(rez.meanXcorrDffLick.v1, rez.semXcorrDffLick.v1, rez.xcorrDffLick.V1{1, 2}, v1_color, {'xcorr v1'}); 
xlabel('Time (s)'); ylabel('xcorr'); set(gca, 'TickDir', 'out', 'XTick', -2000:1000:2000); grid on

print(fullfile(filePath, 'Figure', 'xcorr_dff_lickbouts_v1.pdf'), '-dpdf', '-vector', '-bestfit')

% [rez.meanXcorrDffLickBox.v1, ~, rez.semXcorrDffLickBox.v1] = meanstdsem(cell2mat(rez.xcorrDffLickBox.v1(:, 1))); 
% plotMeanSemColor(rez.meanXcorrDffLickBox.v1, rez.semXcorrDffLickBox.v1, rez.xcorrDffLickBox.v1{1, 2}, {'xcorr v1'}); 
% xlabel('Time (s)'); ylabel('xcorr'); set(gca, 'TickDir', 'out', 'XTick', -2000:1000:2000); 
% print(fullfile(filePath, 'Figure', 'xcorr_dff_lickbouts_v1.pdf'), '-dpdf', '-vector', '-bestfit')
% 
% [rez.meanXcorrDsDffLickBox.v1, ~, rez.semXcorrDsDffLickBox.v1] = meanstdsem(cell2mat(rez.xcorrDsDffLickBox.v1(:, 1))); 
% plotMeanSemColor(rez.meanXcorrDsDffLickBox.v1, rez.semXcorrDsDffLickBox.v1, -2:0.05:2, {'xcorr v1'}); 
% xlabel('Time (s)'); ylabel('xcorr'); set(gca, 'TickDir', 'out', 'XTick', -2:1:2); 
% print(fullfile(filePath, 'Figure', 'xcorr_downsampled_dff_lickbouts_v1.pdf'), '-dpdf', '-vector', '-bestfit')

%% example lick and dff traces
exampleTrials = [1, 10, 11, 21, 26]; 
x1 = 1; 
figure; hold on; 
for tt = 1:length(exampleTrials) 
    
    t = exampleTrials(tt); 
    licks = x1+find(rez.lickOnTfBin{t, 1}); 
    xend = x1+length(rez.dffsOnTfItpV1{t, 1})-1; 

    plot(x1:x1+length(rez.dffsOnTfItpV1{t, 1})-1, rez.dffsOnTfItpV1{t, 1}, 'Color', v1_color, 'LineWidth', 1); 
    vertline(licks, [0 1], 'lineWidth', 1, 'alpha', 0.5); 
    vertline(xend-length(rez.dffsOnTfItpV1{t, 1})/2, ylim, 'lineWidth', 0.75, 'lineStyle', ':', 'alpha', 0.5); 
    
    x1 = x1 + length(rez.dffsOnTfItpV1{t, 1}) + 1; 
end
title("dff and lick traces");
set(gca, 'TickDir', 'out'); 
axis tight
print(fullfile(filePath, 'Figure', 'dff_lickbouts_v1_example5trials.pdf'), '-dpdf', '-vector', '-bestfit')

%% Stim onset aligned activity (V1)
% align dffs take the mean and sem
% primary motor (stim onset)
[stimOnDff_v1_itp, timepts] = temporalAlignInterp1(rez.stimOnDffC.v1(:, 1), rez.stimOnDffC.v1(:, 2));
[rez.meanStimOnDff.v1, ~, rez.semStimOnDff.v1] = meanstdsem(cell2mat(stimOnDff_v1_itp));

% primary motor (stim onset) go no-go
[rez.meanStimOnDffGng.v1, rez.semStimOnDffGng.v1] = trialGroupMeanSem(stimOnDff_v1_itp, {[tbytDat.rewardTrI], [tbytDat.punishTrI]});
plotMeanSem(rez.meanStimOnDffGng.v1, rez.semStimOnDffGng.v1, timepts, {'Go', 'NoGo'});
title("Go Nogo V1 (cue onset at time=0)");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:4, 'TickDir', 'out', 'YTick', -1:0.5:1); xlim([-2 2]); ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_cueOn_Gng_v1.pdf'), '-dpdf', '-vector', '-bestfit')

% primary motor (stim onset) go no-go (logistic regrev1ion)
rez.logitGngCueOnV1 = logireg(cell2mat(stimOnDff_v1_itp), [tbytDat.rewardTrI]); % logistic regrev1ion

%% iti licks (V1)
itiLickDffV1 = flatCell(rez.itiLickDff.v1); % unnest the cell arrays
itiLickDffV1Ts = cellfun(@(a) linspace(-1, 1, length(a)), itiLickDffV1(:, 2), 'UniformOutput', false); 
[itiLickDffV1Itp, itiLickDffV1ItpTs] = temporalAlignInterp1(itiLickDffV1 (:, 1), itiLickDffV1Ts);
[rez.meanItiLickDff.v1, ~, rez.semItiLickDff.v1] = meanstdsem(cell2mat(itiLickDffV1Itp));

plotMeanSemColor(rez.meanItiLickDff.v1, rez.semItiLickDff.v1, itiLickDffV1ItpTs, v1_color, {'V1'});
title("lick during ITI");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); ylim([-.3 .3])

%% postStim licks (V1)
postStimLickDffV1 = flatCell(rez.postStimLickDff.v1); % unnest the cell arrays 
postStimLickDffV1Ts = cellfun(@(a) linspace(-1, 1, length(a)), postStimLickDffV1(:, 2), 'UniformOutput', false); 
[postStimLickDffV1Itp, postStimLickDffV1ItpTs] = temporalAlignInterp1(postStimLickDffV1 (:, 1), postStimLickDffV1Ts);
[rez.meanPostStimLickDff.v1, ~, rez.semPostStimLickDff.v1] = meanstdsem(cell2mat(postStimLickDffV1Itp));

plotMeanSemColor(rez.meanPostStimLickDff.v1, rez.semPostStimLickDff.v1, postStimLickDffV1ItpTs, v1_color, {'V1'});
title("lick during postStim");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); ylim([-.3 .3])

%% hit first licks (V1)
hitFstLickDffV1 = cellWithNonEmptyColumns(rez.hitDffFirstLick.v1); 
hitFstLickDffV1Ts = cellfun(@(a) linspace(-1, 1, length(a)), hitFstLickDffV1(:, 2), 'UniformOutput', false); 
[hitFstLickDffV1Itp, hitFstLickDffV1ItpTs] = temporalAlignInterp1(hitFstLickDffV1 (:, 1), hitFstLickDffV1Ts);
[rez.meanHitFstLickDff.v1, ~, rez.semHitFstLickDff.v1] = meanstdsem(cell2mat(hitFstLickDffV1Itp));

plotMeanSemColor(rez.meanHitFstLickDff.v1, rez.semHitFstLickDff.v1, hitFstLickDffV1ItpTs, v1_color, {'V1'});
title("Lick (first) Hit trial");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); ylim([-1 1])

%% first water (V1)
fstWaterDffV1 = cellWithNonEmptyColumns(rez.firstWater.v1); 
fstWaterDffV1Ts = cellfun(@(a) linspace(-1, 1, length(a)), fstWaterDffV1(:, 2), 'UniformOutput', false); 
[fstWaterDffV1Itp, fstWaterDffV1ItpTs] = temporalAlignInterp1(fstWaterDffV1 (:, 1), fstWaterDffV1Ts);
[rez.meanFstWaterDff.v1, ~, rez.semFstWaterDff.v1] = meanstdsem(cell2mat(fstWaterDffV1Itp));

plotMeanSemColor(rez.meanFstWaterDff.v1, rez.semFstWaterDff.v1, fstWaterDffV1ItpTs, v1_color, {'V1'});
title("Water (first)");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); ylim([-0.8 0.8]); grid on; 

print(fullfile(filePath, 'Figure', 'dff_firstWater_v1.pdf'), '-dpdf', '-vector', '-bestfit')

%% first air (V1)
fstAirDffV1 = cellWithNonEmptyColumns(rez.firstAirpuff.v1); 
fstAirDffV1Ts = cellfun(@(a) linspace(-1, 1, length(a)), fstAirDffV1(:, 2), 'UniformOutput', false); 
[fstAirDffV1Itp, fstAirDffV1ItpTs] = temporalAlignInterp1(fstAirDffV1 (:, 1), fstAirDffV1Ts);
[rez.meanFstAirDff.v1, ~, rez.semFstAirDff.v1] = meanstdsem(cell2mat(fstAirDffV1Itp));

plotMeanSemColor(rez.meanFstAirDff.v1, rez.semFstAirDff.v1, fstAirDffV1ItpTs, v1_color, {'V1'});
title("Air (first)");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); ylim([-0.8 0.8])

print(fullfile(filePath, 'Figure', 'dff_firstAirpuff_v1.pdf'), '-dpdf', '-vector', '-bestfit')

%% FA first licks (V1)
faFstLickDffV1 = cellWithNonEmptyColumns(rez.faDffFirstLick.v1); 
faFstLickDffV1Ts = cellfun(@(a) linspace(-1, 1, length(a)), faFstLickDffV1(:, 2), 'UniformOutput', false); 
[faFstLickDffV1Itp, faFstLickDffV1ItpTs] = temporalAlignInterp1(faFstLickDffV1 (:, 1), faFstLickDffV1Ts);
[rez.meanFaFstLickDff.v1, ~, rez.semFaFstLickDff.v1] = meanstdsem(cell2mat(faFstLickDffV1Itp));

plotMeanSem(rez.meanFaFstLickDff.v1, rez.semFaFstLickDff.v1, faFstLickDffV1ItpTs, {'V1'});
title("Lick (first) Fa trial");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); ylim([-1 1])

%% Compare licks during Hit vs FA vs ITI vs postStim (V1)
% Without normalization
plotMeanSem([rez.meanHitFstLickDff.v1; rez.meanFaFstLickDff.v1; rez.meanItiLickDff.v1; rez.meanPostStimLickDff.v1], ...
    [rez.semHitFstLickDff.v1; rez.semFaFstLickDff.v1; rez.semItiLickDff.v1; rez.semPostStimLickDff.v1], ...
    hitFstLickDffV1ItpTs, {'Hit', 'FA', 'ITI', 'postStim'})
title("Lick Hit vs FA vs ITI vs postStim")
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); ylim([-.8 1]); grid on
print(fullfile(filePath, 'Figure', 'dff_lick_Hit_Fa_ITI_postStim_v1.pdf'), '-dpdf', '-vector', '-bestfit')

% With normalization
rez.meanHitFstLickDffBaseSub.v1 = baseSubNorm(rez.meanHitFstLickDff.v1, hitFstLickDffV1ItpTs, [-1 -0.5]); 
rez.meanFaFstLickDffBaseSub.v1 = baseSubNorm(rez.meanFaFstLickDff.v1, faFstLickDffV1ItpTs, [-1 -0.5]); 
rez.meanItiLickDffBaseSub.v1 = baseSubNorm(rez.meanItiLickDff.v1, itiLickDffV1ItpTs, [-1 -0.5]); 
rez.meanPostStimLickDffBaseSub.v1 = baseSubNorm(rez.meanPostStimLickDff.v1, postStimLickDffV1ItpTs, [-1 -0.5]); 

plotMeanSem([rez.meanHitFstLickDffBaseSub.v1; rez.meanFaFstLickDffBaseSub.v1; rez.meanItiLickDffBaseSub.v1; rez.meanPostStimLickDffBaseSub.v1], ...
    [rez.semHitFstLickDff.v1; rez.semFaFstLickDff.v1; rez.semItiLickDff.v1; rez.semPostStimLickDff.v1], ...
    hitFstLickDffV1ItpTs, {'Hit', 'FA', 'ITI', 'postStim'})
title("Lick Hit vs FA vs ITI vs postStim (baseSub)")
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); ylim([-.6 1.2]); grid on
print(fullfile(filePath, 'Figure', 'dff_lick_Hit_Fa_ITI_postStim_baseSub_v1.pdf'), '-dpdf', '-vector', '-bestfit')

%% train svm or naive bayes classifier to classify licks at different task epochs
X_hit = cell2mat(hitFstLickDffV1Itp); 
X_fa = cell2mat(faFstLickDffV1Itp); 
X_postStim = cell2mat(postStimLickDffV1Itp); 
X_iti = cell2mat(itiLickDffV1Itp); 

Xs = [X_hit; X_fa; X_postStim; X_iti];
y = [ones(size(X_hit, 1), 1)*1; ones(size(X_fa, 1), 1)*2; ones(size(X_postStim, 1), 1)*3; ones(size(X_iti, 1), 1)*3];
rez.lickHitFaRestSvm.v1 = multiClass_svm_peth(Xs, y, 10); 
rez.lickHitFaRestNb.v1 = multiClass_svm_peth(Xs, y, 10); 

Xs = [X_hit; X_fa];
y = [ones(size(X_hit, 1), 1)*1; ones(size(X_fa, 1), 1)*2];
rez.lickHitFaSvm.v1 = multiClass_svm_peth(Xs, y, 10); 
rez.lickHitFaNb.v1 = multiClass_svm_peth(Xs, y, 10); 

figure; hold on; 
plot(mean(rez.lickHitFaRestSvm.v1)); 
plot(mean(rez.lickHitFaRestNb.v1)); 

figure; hold on; 
plot(mean(rez.lickHitFaSvm.v1)); 
plot(mean(rez.lickHitFaNb.v1)); 

%% Miss CueOn (V1)
missCueDffV1 = cellWithNonEmptyColumns(rez.missCueDff.v1); 
missCueDffV1Ts = cellfun(@(a) linspace(-1, 1, length(a)), missCueDffV1(:, 2), 'UniformOutput', false); 
[missCueDffV1Itp, missCueDffV1ItpTs] = temporalAlignInterp1(missCueDffV1 (:, 1), missCueDffV1Ts);
[rez.meanMiv1CueDff.v1, ~, rez.semMiv1CueDff.v1] = meanstdsem(cell2mat(missCueDffV1Itp));

plotMeanSem(rez.meanMiv1CueDff.v1, rez.semMiv1CueDff.v1, missCueDffV1ItpTs, {'V1'});
title("Miv1 CueOn trial");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); ylim([-1 1])

print(fullfile(filePath, 'Figure', 'dff_miss_cueOn_v1.pdf'), '-dpdf', '-vector', '-bestfit')

%% CR CueOn (V1)
crCueDffV1 = cellWithNonEmptyColumns(rez.crCueDff.v1); 
crCueDffV1Ts = cellfun(@(a) linspace(-1, 1, length(a)), crCueDffV1(:, 2), 'UniformOutput', false); 
[crCueDffV1Itp, crCueDffV1ItpTs] = temporalAlignInterp1(crCueDffV1 (:, 1), crCueDffV1Ts);
[rez.meanCrCueDff.v1, ~, rez.semCrCueDff.v1] = meanstdsem(cell2mat(crCueDffV1Itp));

plotMeanSem(rez.meanCrCueDff.v1, rez.semCrCueDff.v1, crCueDffV1ItpTs, {'V1'});
title("Correct rejection CueOn trial");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); ylim([-1 1])

print(fullfile(filePath, 'Figure', 'dff_cr_cueOn_v1.pdf'), '-dpdf', '-vector', '-bestfit')

%% save rez
save(fullfile(filePath, 'Matfiles', strcat(header, '_dff_evtAligned_regionMask.mat')), 'rez'); 

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
        % 'tbytPsthC' contains trial-by-trial data in each cell (it is av1umed that data are temporally aligned acrov1 trials already).
        % 'groupC' contains logical for each grouping, logicals are av1umed to be of same lengths as the number of trials.
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

        % Determine the number of samples for each clav1 in the test set
        numClav10 = sum(Y == 0);
        numClav11 = sum(Y == 1);
        testSetSize = min(numClav10, numClav11) / 2; % Or any other criterion

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




end


