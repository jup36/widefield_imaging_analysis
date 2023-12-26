function dffPostprocessGng_averaging_m2(filePath)
%filePath = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA008/DA008_101823';
% This function performs further analyses and plotting using data in
% 'rez'(outcome from 'dffPostprocessGng_averaging.m'). 
% 

%% load rez tbytDat
[~, header] = fileparts(filePath);
load(fullfile(filePath, 'Matfiles', strcat(header, '_dff_evtAligned_regionMask.mat')), 'rez'); 
load(fullfile(filePath, 'Matfiles', strcat(header, '_tbytDat_pam2eGng')), 'tbytDat')

% user variables 
m2_color = [64 191 255]./255; 

%% cross correlogram lick bouts and dff
[rez.meanXcorrDffLick.m2, ~, rez.semXcorrDffLick.m2] = meanstdsem(cell2mat(rez.xcorrDffLick.m2(:, 1))); 
plotMeanSemColor(rez.meanXcorrDffLick.m2, rez.semXcorrDffLick.m2, rez.xcorrDffLick.m2{1, 2}, m2_color, {'xcorr m2'}); 
xlabel('Time (s)'); ylabel('xcorr'); set(gca, 'TickDir', 'out', 'XTick', -2000:1000:2000); grid on

print(fullfile(filePath, 'Figure', 'xcorr_dff_lickbouts_m2.pdf'), '-dpdf', '-vector', '-bestfit')

% [rez.meanXcorrDffLickBox.m2, ~, rez.semXcorrDffLickBox.m2] = meanstdsem(cell2mat(rez.xcorrDffLickBox.m2(:, 1))); 
% plotMeanSemColor(rez.meanXcorrDffLickBox.m2, rez.semXcorrDffLickBox.m2, rez.xcorrDffLickBox.m2{1, 2}, {'xcorr m2'}); 
% xlabel('Time (s)'); ylabel('xcorr'); set(gca, 'TickDir', 'out', 'XTick', -2000:1000:2000); 
% print(fullfile(filePath, 'Figure', 'xcorr_dff_lickbouts_m2.pdf'), '-dpdf', '-vector', '-bestfit')
% 
% [rez.meanXcorrDsDffLickBox.m2, ~, rez.semXcorrDsDffLickBox.m2] = meanstdsem(cell2mat(rez.xcorrDsDffLickBox.m2(:, 1))); 
% plotMeanSemColor(rez.meanXcorrDsDffLickBox.m2, rez.semXcorrDsDffLickBox.m2, -2:0.05:2, {'xcorr m2'}); 
% xlabel('Time (s)'); ylabel('xcorr'); set(gca, 'TickDir', 'out', 'XTick', -2:1:2); 
% print(fullfile(filePath, 'Figure', 'xcorr_downsampled_dff_lickbouts_m2.pdf'), '-dpdf', '-vector', '-bestfit')

%% example lick and dff traces
exampleTrials = [1, 10, 11, 21, 26]; 
x1 = 1; 
figure; hold on; 
for tt = 1:length(exampleTrials) 
    
    t = exampleTrials(tt); 
    licks = x1+find(rez.lickOnTfBin{t, 1}); 
    xend = x1+length(rez.dffsOnTfItpM2{t, 1})-1; 

    plot(x1:x1+length(rez.dffsOnTfItpM2{t, 1})-1, rez.dffsOnTfItpM2{t, 1}, 'Color', m2_color, 'LineWidth', 1); 
    vertline(licks, [0 1], 'lineWidth', 1, 'alpha', 0.5); 
    vertline(xend-length(rez.dffsOnTfItpM2{t, 1})/2, ylim, 'lineWidth', 0.75, 'lineStyle', ':', 'alpha', 0.5); 
    
    x1 = x1 + length(rez.dffsOnTfItpM2{t, 1}) + 1; 
end
title("dff and lick traces");
set(gca, 'TickDir', 'out'); 
axis tight
print(fullfile(filePath, 'Figure', 'dff_lickbouts_m2_example5trials.pdf'), '-dpdf', '-vector', '-bestfit')

%% Stim onset aligned activity (M2)
% align dffs take the mean and sem
% primary motor (stim onset)
[stimOnDff_m2_itp, timepts] = temporalAlignInterp1(rez.stimOnDffC.m2(:, 1), rez.stimOnDffC.m2(:, 2));
[rez.meanStimOnDff.m2, ~, rez.semStimOnDff.m2] = meanstdsem(cell2mat(stimOnDff_m2_itp));

% primary motor (stim onset) go no-go
[rez.meanStimOnDffGng.m2, rez.semStimOnDffGng.m2] = trialGroupMeanSem(stimOnDff_m2_itp, {[tbytDat.rewardTrI], [tbytDat.punishTrI]});
plotMeanSem(rez.meanStimOnDffGng.m2, rez.semStimOnDffGng.m2, timepts, {'Go', 'NoGo'});
title("Go Nogo M2 (cue onset at time=0)");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:4, 'TickDir', 'out', 'YTick', -1:0.5:1); xlim([-2 2]); ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_cueOn_Gng_m2.pdf'), '-dpdf', '-vector', '-bestfit')

% primary motor (stim onset) go no-go (logistic regression)
rez.logitGngCueOnM2 = logireg(cell2mat(stimOnDff_m2_itp), [tbytDat.rewardTrI]); % logistic regression

%% iti licks (M2)
itiLickDffM2 = flatCell(rez.itiLickDff.m2); % unnest the cell arrays
itiLickDffM2Ts = cellfun(@(a) linspace(-1, 1, length(a)), itiLickDffM2(:, 2), 'UniformOutput', false); 
[itiLickDffM2Itp, itiLickDffM2ItpTs] = temporalAlignInterp1(itiLickDffM2 (:, 1), itiLickDffM2Ts);
[rez.meanItiLickDff.m2, ~, rez.semItiLickDff.m2] = meanstdsem(cell2mat(itiLickDffM2Itp));

plotMeanSemColor(rez.meanItiLickDff.m2, rez.semItiLickDff.m2, itiLickDffM2ItpTs, m2_color, {'M2'});
title("lick during ITI");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); ylim([-.3 .3])

%% postStim licks (M2)
postStimLickDffM2 = flatCell(rez.postStimLickDff.m2); % unnest the cell arrays 
postStimLickDffM2Ts = cellfun(@(a) linspace(-1, 1, length(a)), postStimLickDffM2(:, 2), 'UniformOutput', false); 
[postStimLickDffM2Itp, postStimLickDffM2ItpTs] = temporalAlignInterp1(postStimLickDffM2 (:, 1), postStimLickDffM2Ts);
[rez.meanPostStimLickDff.m2, ~, rez.semPostStimLickDff.m2] = meanstdsem(cell2mat(postStimLickDffM2Itp));

plotMeanSemColor(rez.meanPostStimLickDff.m2, rez.semPostStimLickDff.m2, postStimLickDffM2ItpTs, m2_color, {'M2'});
title("lick during ITI");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); ylim([-.3 .3])

%% hit fim2t licks (M2)
hitFstLickDffM2 = cellWithNonEmptyColumns(rez.hitDffFim2tLick.m2); 
hitFstLickDffM2Ts = cellfun(@(a) linspace(-1, 1, length(a)), hitFstLickDffM2(:, 2), 'UniformOutput', false); 
[hitFstLickDffM2Itp, hitFstLickDffM2ItpTs] = temporalAlignInterp1(hitFstLickDffM2 (:, 1), hitFstLickDffM2Ts);
[rez.meanHitFstLickDff.m2, ~, rez.semHitFstLickDff.m2] = meanstdsem(cell2mat(hitFstLickDffM2Itp));

plotMeanSemColor(rez.meanHitFstLickDff.m2, rez.semHitFstLickDff.m2, hitFstLickDffM2ItpTs, m2_color, {'M2'});
title("Lick (fim2t) Hit trial");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); ylim([-1 1])

%% fim2t water (M2)
fstWaterDffM2 = cellWithNonEmptyColumns(rez.fim2tWater.m2); 
fstWaterDffM2Ts = cellfun(@(a) linspace(-1, 1, length(a)), fstWaterDffM2(:, 2), 'UniformOutput', false); 
[fstWaterDffM2Itp, fstWaterDffM2ItpTs] = temporalAlignInterp1(fstWaterDffM2 (:, 1), fstWaterDffM2Ts);
[rez.meanFstWaterDff.m2, ~, rez.semFstWaterDff.m2] = meanstdsem(cell2mat(fstWaterDffM2Itp));

plotMeanSemColor(rez.meanFstWaterDff.m2, rez.semFstWaterDff.m2, fstWaterDffM2ItpTs, m2_color, {'M2'});
title("Water (fim2t)");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); ylim([-0.8 0.8]); grid on; 

print(fullfile(filePath, 'Figure', 'dff_fim2tWater_m2.pdf'), '-dpdf', '-vector', '-bestfit')

%% fim2t air (M2)
fstAirDffM2 = cellWithNonEmptyColumns(rez.fim2tAirpuff.m2); 
fstAirDffM2Ts = cellfun(@(a) linspace(-1, 1, length(a)), fstAirDffM2(:, 2), 'UniformOutput', false); 
[fstAirDffM2Itp, fstAirDffM2ItpTs] = temporalAlignInterp1(fstAirDffM2 (:, 1), fstAirDffM2Ts);
[rez.meanFstAirDff.m2, ~, rez.semFstAirDff.m2] = meanstdsem(cell2mat(fstAirDffM2Itp));

plotMeanSemColor(rez.meanFstAirDff.m2, rez.semFstAirDff.m2, fstAirDffM2ItpTs, m2_color, {'M2'});
title("Air (fim2t)");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); ylim([-0.8 0.8])

print(fullfile(filePath, 'Figure', 'dff_fim2tAirpuff_m2.pdf'), '-dpdf', '-vector', '-bestfit')

%% FA fim2t licks (M2)
faFstLickDffM2 = cellWithNonEmptyColumns(rez.faDffFim2tLick.m2); 
faFstLickDffM2Ts = cellfun(@(a) linspace(-1, 1, length(a)), faFstLickDffM2(:, 2), 'UniformOutput', false); 
[faFstLickDffM2Itp, faFstLickDffM2ItpTs] = temporalAlignInterp1(faFstLickDffM2 (:, 1), faFstLickDffM2Ts);
[rez.meanFaFstLickDff.m2, ~, rez.semFaFstLickDff.m2] = meanstdsem(cell2mat(faFstLickDffM2Itp));

plotMeanSem(rez.meanFaFstLickDff.m2, rez.semFaFstLickDff.m2, faFstLickDffM2ItpTs, {'M2'});
title("Lick (fim2t) Fa trial");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); ylim([-1 1])

%% Compare licks during Hit vs FA vs ITI vs postStim (M2)
% Without normalization
plotMeanSem([rez.meanHitFstLickDff.m2; rez.meanFaFstLickDff.m2; rez.meanItiLickDff.m2; rez.meanPostStimLickDff.m2], ...
    [rez.semHitFstLickDff.m2; rez.semFaFstLickDff.m2; rez.semItiLickDff.m2; rez.semPostStimLickDff.m2], ...
    hitFstLickDffM2ItpTs, {'Hit', 'FA', 'ITI', 'postStim'})
title("Lick Hit vs FA vs ITI vs postStim")
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); ylim([-.8 1]); grid on
print(fullfile(filePath, 'Figure', 'dff_lick_Hit_Fa_ITI_postStim_m2.pdf'), '-dpdf', '-vector', '-bestfit')

% With normalization
rez.meanHitFstLickDffBaseSub.m2 = baseSubNorm(rez.meanHitFstLickDff.m2, hitFstLickDffM2ItpTs, [-1 -0.5]); 
rez.meanFaFstLickDffBaseSub.m2 = baseSubNorm(rez.meanFaFstLickDff.m2, faFstLickDffM2ItpTs, [-1 -0.5]); 
rez.meanItiLickDffBaseSub.m2 = baseSubNorm(rez.meanItiLickDff.m2, itiLickDffM2ItpTs, [-1 -0.5]); 
rez.meanPostStimLickDffBaseSub.m2 = baseSubNorm(rez.meanPostStimLickDff.m2, postStimLickDffM2ItpTs, [-1 -0.5]); 

plotMeanSem([rez.meanHitFstLickDffBaseSub.m2; rez.meanFaFstLickDffBaseSub.m2; rez.meanItiLickDffBaseSub.m2; rez.meanPostStimLickDffBaseSub.m2], ...
    [rez.semHitFstLickDff.m2; rez.semFaFstLickDff.m2; rez.semItiLickDff.m2; rez.semPostStimLickDff.m2], ...
    hitFstLickDffM2ItpTs, {'Hit', 'FA', 'ITI', 'postStim'})
title("Lick Hit vs FA vs ITI vs postStim (baseSub)")
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); ylim([-.6 1.2]); grid on
print(fullfile(filePath, 'Figure', 'dff_lick_Hit_Fa_ITI_postStim_baseSub_m2.pdf'), '-dpdf', '-vector', '-bestfit')

%% train svm or naive bayes classifier to classify licks at different task epochs
X_hit = cell2mat(hitFstLickDffM2Itp); 
X_fa = cell2mat(faFstLickDffM2Itp); 
X_postStim = cell2mat(postStimLickDffM2Itp); 
X_iti = cell2mat(itiLickDffM2Itp); 

Xs = [X_hit; X_fa; X_postStim; X_iti];
y = [ones(size(X_hit, 1), 1)*1; ones(size(X_fa, 1), 1)*2; ones(size(X_postStim, 1), 1)*3; ones(size(X_iti, 1), 1)*3];
rez.lickHitFaRestSvm.m2 = multiClass_svm_peth(Xs, y, 10); 
rez.lickHitFaRestNb.m2 = multiClass_svm_peth(Xs, y, 10); 

Xs = [X_hit; X_fa];
y = [ones(size(X_hit, 1), 1)*1; ones(size(X_fa, 1), 1)*2];
rez.lickHitFaSvm.m2 = multiClass_svm_peth(Xs, y, 10); 
rez.lickHitFaNb.m2 = multiClass_svm_peth(Xs, y, 10); 

figure; hold on; 
plot(mean(rez.lickHitFaRestSvm.m2)); 
plot(mean(rez.lickHitFaRestNb.m2)); 

figure; hold on; 
plot(mean(rez.lickHitFaSvm.m2)); 
plot(mean(rez.lickHitFaNb.m2)); 

%% Miss CueOn (M2)
missCueDffM2 = cellWithNonEmptyColumns(rez.missCueDff.m2); 
missCueDffM2Ts = cellfun(@(a) linspace(-1, 1, length(a)), missCueDffM2(:, 2), 'UniformOutput', false); 
[missCueDffM2Itp, missCueDffM2ItpTs] = temporalAlignInterp1(missCueDffM2 (:, 1), missCueDffM2Ts);
[rez.meanMissCueDff.m2, ~, rez.semMissCueDff.m2] = meanstdsem(cell2mat(missCueDffM2Itp));

plotMeanSem(rez.meanMissCueDff.m2, rez.semMissCueDff.m2, missCueDffM2ItpTs, {'M2'});
title("Miss CueOn trial");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); ylim([-1 1])

print(fullfile(filePath, 'Figure', 'dff_miss_cueOn_m2.pdf'), '-dpdf', '-vector', '-bestfit')

%% CR CueOn (M2)
crCueDffM2 = cellWithNonEmptyColumns(rez.crCueDff.m2); 
crCueDffM2Ts = cellfun(@(a) linspace(-1, 1, length(a)), crCueDffM2(:, 2), 'UniformOutput', false); 
[crCueDffM2Itp, crCueDffM2ItpTs] = temporalAlignInterp1(crCueDffM2 (:, 1), crCueDffM2Ts);
[rez.meanCrCueDff.m2, ~, rez.semCrCueDff.m2] = meanstdsem(cell2mat(crCueDffM2Itp));

plotMeanSem(rez.meanCrCueDff.m2, rez.semCrCueDff.m2, crCueDffM2ItpTs, {'M2'});
title("Correct rejection CueOn trial");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); ylim([-1 1])

print(fullfile(filePath, 'Figure', 'dff_cr_cueOn_m2.pdf'), '-dpdf', '-vector', '-bestfit')

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




end


