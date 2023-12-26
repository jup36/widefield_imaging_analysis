function dffPostprocessGng_averaging_ss(filePath)
%filePath = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA008/DA008_101823';
% This function performs further analyses and plotting using data in
% 'rez'(outcome from 'dffPostprocessGng_averaging.m'). 
% 

%% load rez tbytDat
[~, header] = fileparts(filePath);
load(fullfile(filePath, 'Matfiles', strcat(header, '_dff_evtAligned_regionMask.mat')), 'rez'); 
load(fullfile(filePath, 'Matfiles', strcat(header, '_tbytDat_parseGng')), 'tbytDat')

% user variables 
ss_color = [128 127 255]./255; 

%% cross correlogram lick bouts and dff
[rez.meanXcorrDffLick.ss, ~, rez.semXcorrDffLick.ss] = meanstdsem(cell2mat(rez.xcorrDffLick.ss(:, 1))); 
plotMeanSemColor(rez.meanXcorrDffLick.ss, rez.semXcorrDffLick.ss, rez.xcorrDffLick.ss{1, 2}, ss_color, {'xcorr ss'}); 
xlabel('Time (s)'); ylabel('xcorr'); set(gca, 'TickDir', 'out', 'XTick', -2000:1000:2000); grid on

print(fullfile(filePath, 'Figure', 'xcorr_dff_lickbouts_ss.pdf'), '-dpdf', '-vector', '-bestfit')

% [rez.meanXcorrDffLickBox.ss, ~, rez.semXcorrDffLickBox.ss] = meanstdsem(cell2mat(rez.xcorrDffLickBox.ss(:, 1))); 
% plotMeanSemColor(rez.meanXcorrDffLickBox.ss, rez.semXcorrDffLickBox.ss, rez.xcorrDffLickBox.ss{1, 2}, {'xcorr ss'}); 
% xlabel('Time (s)'); ylabel('xcorr'); set(gca, 'TickDir', 'out', 'XTick', -2000:1000:2000); 
% print(fullfile(filePath, 'Figure', 'xcorr_dff_lickbouts_ss.pdf'), '-dpdf', '-vector', '-bestfit')
% 
% [rez.meanXcorrDsDffLickBox.ss, ~, rez.semXcorrDsDffLickBox.ss] = meanstdsem(cell2mat(rez.xcorrDsDffLickBox.ss(:, 1))); 
% plotMeanSemColor(rez.meanXcorrDsDffLickBox.ss, rez.semXcorrDsDffLickBox.ss, -2:0.05:2, {'xcorr ss'}); 
% xlabel('Time (s)'); ylabel('xcorr'); set(gca, 'TickDir', 'out', 'XTick', -2:1:2); 
% print(fullfile(filePath, 'Figure', 'xcorr_downsampled_dff_lickbouts_ss.pdf'), '-dpdf', '-vector', '-bestfit')

%% example lick and dff traces
exampleTrials = [1, 10, 11, 21, 26]; 
x1 = 1; 
figure; hold on; 
for tt = 1:length(exampleTrials) 
    
    t = exampleTrials(tt); 
    licks = x1+find(rez.lickOnTfBin{t, 1}); 
    xend = x1+length(rez.dffsOnTfItpSs{t, 1})-1; 

    plot(x1:x1+length(rez.dffsOnTfItpSs{t, 1})-1, rez.dffsOnTfItpSs{t, 1}, 'Color', ss_color, 'LineWidth', 1); 
    vertline(licks, [0 1], 'lineWidth', 1, 'alpha', 0.5); 
    vertline(xend-length(rez.dffsOnTfItpSs{t, 1})/2, ylim, 'lineWidth', 0.75, 'lineStyle', ':', 'alpha', 0.5); 
    
    x1 = x1 + length(rez.dffsOnTfItpSs{t, 1}) + 1; 
end
title("dff and lick traces");
set(gca, 'TickDir', 'out'); 
axis tight
print(fullfile(filePath, 'Figure', 'dff_lickbouts_ss_example5trials.pdf'), '-dpdf', '-vector', '-bestfit')

%% Stim onset aligned activity (Ss)
% align dffs take the mean and sem
% primary motor (stim onset)
[stimOnDff_ss_itp, timepts] = temporalAlignInterp1(rez.stimOnDffC.ss(:, 1), rez.stimOnDffC.ss(:, 2));
[rez.meanStimOnDff.ss, ~, rez.semStimOnDff.ss] = meanstdsem(cell2mat(stimOnDff_ss_itp));

% primary motor (stim onset) go no-go
[rez.meanStimOnDffGng.ss, rez.semStimOnDffGng.ss] = trialGroupMeanSem(stimOnDff_ss_itp, {[tbytDat.rewardTrI], [tbytDat.punishTrI]});
plotMeanSem(rez.meanStimOnDffGng.ss, rez.semStimOnDffGng.ss, timepts, {'Go', 'NoGo'});
title("Go Nogo Ss (cue onset at time=0)");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:4, 'TickDir', 'out', 'YTick', -1:0.5:1); xlim([-2 2]); ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_cueOn_Gng_ss.pdf'), '-dpdf', '-vector', '-bestfit')

% primary motor (stim onset) go no-go (logistic regression)
rez.logitGngCueOnSs = logireg(cell2mat(stimOnDff_ss_itp), [tbytDat.rewardTrI]); % logistic regression

%% iti licks (Ss)
itiLickDffSs = flatCell(rez.itiLickDff.ss); % unnest the cell arrays
itiLickDffSsTs = cellfun(@(a) linspace(-1, 1, length(a)), itiLickDffSs(:, 2), 'UniformOutput', false); 
[itiLickDffSsItp, itiLickDffSsItpTs] = temporalAlignInterp1(itiLickDffSs (:, 1), itiLickDffSsTs);
[rez.meanItiLickDff.ss, ~, rez.semItiLickDff.ss] = meanstdsem(cell2mat(itiLickDffSsItp));

plotMeanSemColor(rez.meanItiLickDff.ss, rez.semItiLickDff.ss, itiLickDffSsItpTs, ss_color, {'Ss'});
title("lick during ITI");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); ylim([-.3 .3])

%% postStim licks (Ss)
postStimLickDffSs = flatCell(rez.postStimLickDff.ss); % unnest the cell arrays 
postStimLickDffSsTs = cellfun(@(a) linspace(-1, 1, length(a)), postStimLickDffSs(:, 2), 'UniformOutput', false); 
[postStimLickDffSsItp, postStimLickDffSsItpTs] = temporalAlignInterp1(postStimLickDffSs (:, 1), postStimLickDffSsTs);
[rez.meanPostStimLickDff.ss, ~, rez.semPostStimLickDff.ss] = meanstdsem(cell2mat(postStimLickDffSsItp));

plotMeanSemColor(rez.meanPostStimLickDff.ss, rez.semPostStimLickDff.ss, postStimLickDffSsItpTs, ss_color, {'Ss'});
title("lick during ITI");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); ylim([-.3 .3])

%% hit first licks (Ss)
hitFstLickDffSs = cellWithNonEmptyColumns(rez.hitDffFirstLick.ss); 
hitFstLickDffSsTs = cellfun(@(a) linspace(-1, 1, length(a)), hitFstLickDffSs(:, 2), 'UniformOutput', false); 
[hitFstLickDffSsItp, hitFstLickDffSsItpTs] = temporalAlignInterp1(hitFstLickDffSs (:, 1), hitFstLickDffSsTs);
[rez.meanHitFstLickDff.ss, ~, rez.semHitFstLickDff.ss] = meanstdsem(cell2mat(hitFstLickDffSsItp));

plotMeanSemColor(rez.meanHitFstLickDff.ss, rez.semHitFstLickDff.ss, hitFstLickDffSsItpTs, ss_color, {'Ss'});
title("Lick (first) Hit trial");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); ylim([-1 1])

%% first water (Ss)
fstWaterDffSs = cellWithNonEmptyColumns(rez.firstWater.ss); 
fstWaterDffSsTs = cellfun(@(a) linspace(-1, 1, length(a)), fstWaterDffSs(:, 2), 'UniformOutput', false); 
[fstWaterDffSsItp, fstWaterDffSsItpTs] = temporalAlignInterp1(fstWaterDffSs (:, 1), fstWaterDffSsTs);
[rez.meanFstWaterDff.ss, ~, rez.semFstWaterDff.ss] = meanstdsem(cell2mat(fstWaterDffSsItp));

plotMeanSemColor(rez.meanFstWaterDff.ss, rez.semFstWaterDff.ss, fstWaterDffSsItpTs, ss_color, {'Ss'});
title("Water (first)");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); ylim([-0.8 0.8]); grid on; 

print(fullfile(filePath, 'Figure', 'dff_firstWater_ss.pdf'), '-dpdf', '-vector', '-bestfit')

%% first air (Ss)
fstAirDffSs = cellWithNonEmptyColumns(rez.firstAirpuff.ss); 
fstAirDffSsTs = cellfun(@(a) linspace(-1, 1, length(a)), fstAirDffSs(:, 2), 'UniformOutput', false); 
[fstAirDffSsItp, fstAirDffSsItpTs] = temporalAlignInterp1(fstAirDffSs (:, 1), fstAirDffSsTs);
[rez.meanFstAirDff.ss, ~, rez.semFstAirDff.ss] = meanstdsem(cell2mat(fstAirDffSsItp));

plotMeanSemColor(rez.meanFstAirDff.ss, rez.semFstAirDff.ss, fstAirDffSsItpTs, ss_color, {'Ss'});
title("Air (first)");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); ylim([-0.8 0.8])

print(fullfile(filePath, 'Figure', 'dff_firstAirpuff_ss.pdf'), '-dpdf', '-vector', '-bestfit')

%% FA first licks (Ss)
faFstLickDffSs = cellWithNonEmptyColumns(rez.faDffFirstLick.ss); 
faFstLickDffSsTs = cellfun(@(a) linspace(-1, 1, length(a)), faFstLickDffSs(:, 2), 'UniformOutput', false); 
[faFstLickDffSsItp, faFstLickDffSsItpTs] = temporalAlignInterp1(faFstLickDffSs (:, 1), faFstLickDffSsTs);
[rez.meanFaFstLickDff.ss, ~, rez.semFaFstLickDff.ss] = meanstdsem(cell2mat(faFstLickDffSsItp));

plotMeanSem(rez.meanFaFstLickDff.ss, rez.semFaFstLickDff.ss, faFstLickDffSsItpTs, {'Ss'});
title("Lick (first) Fa trial");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); ylim([-1 1])

%% Compare licks during Hit vs FA vs ITI vs postStim (Ss)
% Without normalization
plotMeanSem([rez.meanHitFstLickDff.ss; rez.meanFaFstLickDff.ss; rez.meanItiLickDff.ss; rez.meanPostStimLickDff.ss], ...
    [rez.semHitFstLickDff.ss; rez.semFaFstLickDff.ss; rez.semItiLickDff.ss; rez.semPostStimLickDff.ss], ...
    hitFstLickDffSsItpTs, {'Hit', 'FA', 'ITI', 'postStim'})
title("Lick Hit vs FA vs ITI vs postStim")
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); ylim([-.8 1]); grid on
print(fullfile(filePath, 'Figure', 'dff_lick_Hit_Fa_ITI_postStim_ss.pdf'), '-dpdf', '-vector', '-bestfit')

% With normalization
rez.meanHitFstLickDffBaseSub.ss = baseSubNorm(rez.meanHitFstLickDff.ss, hitFstLickDffSsItpTs, [-1 -0.5]); 
rez.meanFaFstLickDffBaseSub.ss = baseSubNorm(rez.meanFaFstLickDff.ss, faFstLickDffSsItpTs, [-1 -0.5]); 
rez.meanItiLickDffBaseSub.ss = baseSubNorm(rez.meanItiLickDff.ss, itiLickDffSsItpTs, [-1 -0.5]); 
rez.meanPostStimLickDffBaseSub.ss = baseSubNorm(rez.meanPostStimLickDff.ss, postStimLickDffSsItpTs, [-1 -0.5]); 

plotMeanSem([rez.meanHitFstLickDffBaseSub.ss; rez.meanFaFstLickDffBaseSub.ss; rez.meanItiLickDffBaseSub.ss; rez.meanPostStimLickDffBaseSub.ss], ...
    [rez.semHitFstLickDff.ss; rez.semFaFstLickDff.ss; rez.semItiLickDff.ss; rez.semPostStimLickDff.ss], ...
    hitFstLickDffSsItpTs, {'Hit', 'FA', 'ITI', 'postStim'})
title("Lick Hit vs FA vs ITI vs postStim (baseSub)")
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); ylim([-.6 1.2]); grid on
print(fullfile(filePath, 'Figure', 'dff_lick_Hit_Fa_ITI_postStim_baseSub_ss.pdf'), '-dpdf', '-vector', '-bestfit')

%% train svm or naive bayes classifier to classify licks at different task epochs
X_hit = cell2mat(hitFstLickDffSsItp); 
X_fa = cell2mat(faFstLickDffSsItp); 
X_postStim = cell2mat(postStimLickDffSsItp); 
X_iti = cell2mat(itiLickDffSsItp); 

Xs = [X_hit; X_fa; X_postStim; X_iti];
y = [ones(size(X_hit, 1), 1)*1; ones(size(X_fa, 1), 1)*2; ones(size(X_postStim, 1), 1)*3; ones(size(X_iti, 1), 1)*3];
rez.lickHitFaRestSvm.ss = multiClass_svm_peth(Xs, y, 10); 
rez.lickHitFaRestNb.ss = multiClass_svm_peth(Xs, y, 10); 

Xs = [X_hit; X_fa];
y = [ones(size(X_hit, 1), 1)*1; ones(size(X_fa, 1), 1)*2];
rez.lickHitFaSvm.ss = multiClass_svm_peth(Xs, y, 10); 
rez.lickHitFaNb.ss = multiClass_svm_peth(Xs, y, 10); 

figure; hold on; 
plot(mean(rez.lickHitFaRestSvm.ss)); 
plot(mean(rez.lickHitFaRestNb.ss)); 

figure; hold on; 
plot(mean(rez.lickHitFaSvm.ss)); 
plot(mean(rez.lickHitFaNb.ss)); 

%% Miss CueOn (Ss)
missCueDffSs = cellWithNonEmptyColumns(rez.missCueDff.ss); 
missCueDffSsTs = cellfun(@(a) linspace(-1, 1, length(a)), missCueDffSs(:, 2), 'UniformOutput', false); 
[missCueDffSsItp, missCueDffSsItpTs] = temporalAlignInterp1(missCueDffSs (:, 1), missCueDffSsTs);
[rez.meanMissCueDff.ss, ~, rez.semMissCueDff.ss] = meanstdsem(cell2mat(missCueDffSsItp));

plotMeanSem(rez.meanMissCueDff.ss, rez.semMissCueDff.ss, missCueDffSsItpTs, {'Ss'});
title("Miss CueOn trial");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); ylim([-1 1])

print(fullfile(filePath, 'Figure', 'dff_miss_cueOn_ss.pdf'), '-dpdf', '-vector', '-bestfit')

%% CR CueOn (Ss)
crCueDffSs = cellWithNonEmptyColumns(rez.crCueDff.ss); 
crCueDffSsTs = cellfun(@(a) linspace(-1, 1, length(a)), crCueDffSs(:, 2), 'UniformOutput', false); 
[crCueDffSsItp, crCueDffSsItpTs] = temporalAlignInterp1(crCueDffSs (:, 1), crCueDffSsTs);
[rez.meanCrCueDff.ss, ~, rez.semCrCueDff.ss] = meanstdsem(cell2mat(crCueDffSsItp));

plotMeanSem(rez.meanCrCueDff.ss, rez.semCrCueDff.ss, crCueDffSsItpTs, {'Ss'});
title("Correct rejection CueOn trial");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); ylim([-1 1])

print(fullfile(filePath, 'Figure', 'dff_cr_cueOn_ss.pdf'), '-dpdf', '-vector', '-bestfit')

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

        % tbytPsthC = stimOnDff_ss_itp;
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
        % %   e.g., cell2mat(stimOnDff_ss_itp);
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


