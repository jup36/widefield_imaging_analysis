function dffPostprocessGng_averaging_rs(filePath)
%filePath = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA008/DA008_101823';
% This function performs further analyses and plotting using data in
% 'rez'(outcome from 'dffPostprocersGng_averaging.m'). 
% 

%% load rez tbytDat
[~, header] = fileparts(filePath);
load(fullfile(filePath, 'Matfiles', strcat(header, '_dff_evtAligned_regionMask.mat')), 'rez'); 
load(fullfile(filePath, 'Matfiles', strcat(header, '_tbytDat_parseGng')), 'tbytDat')

% user variables 
rs_color = [191 64 255]./255; 

%% crors correlogram lick bouts and dff
[rez.meanXcorrDffLick.rs, ~, rez.semXcorrDffLick.rs] = meanstdsem(cell2mat(rez.xcorrDffLick.rs(:, 1))); 
plotMeanSemColor(rez.meanXcorrDffLick.rs, rez.semXcorrDffLick.rs, rez.xcorrDffLick.rs{1, 2}, rs_color, {'xcorr rs'}); 
xlabel('Time (s)'); ylabel('xcorr'); set(gca, 'TickDir', 'out', 'XTick', -2000:1000:2000); grid on

print(fullfile(filePath, 'Figure', 'xcorr_dff_lickbouts_rs.pdf'), '-dpdf', '-vector', '-bestfit')

% [rez.meanXcorrDffLickBox.rs, ~, rez.semXcorrDffLickBox.rs] = meanstdsem(cell2mat(rez.xcorrDffLickBox.rs(:, 1))); 
% plotMeanSemColor(rez.meanXcorrDffLickBox.rs, rez.semXcorrDffLickBox.rs, rez.xcorrDffLickBox.rs{1, 2}, {'xcorr rs'}); 
% xlabel('Time (s)'); ylabel('xcorr'); set(gca, 'TickDir', 'out', 'XTick', -2000:1000:2000); 
% print(fullfile(filePath, 'Figure', 'xcorr_dff_lickbouts_rs.pdf'), '-dpdf', '-vector', '-bestfit')
% 
% [rez.meanXcorrDsDffLickBox.rs, ~, rez.semXcorrDsDffLickBox.rs] = meanstdsem(cell2mat(rez.xcorrDsDffLickBox.rs(:, 1))); 
% plotMeanSemColor(rez.meanXcorrDsDffLickBox.rs, rez.semXcorrDsDffLickBox.rs, -2:0.05:2, {'xcorr rs'}); 
% xlabel('Time (s)'); ylabel('xcorr'); set(gca, 'TickDir', 'out', 'XTick', -2:1:2); 
% print(fullfile(filePath, 'Figure', 'xcorr_downsampled_dff_lickbouts_rs.pdf'), '-dpdf', '-vector', '-bestfit')

%% example lick and dff traces
exampleTrials = [1, 10, 11, 21, 26]; 
x1 = 1; 
figure; hold on; 
for tt = 1:length(exampleTrials) 
    
    t = exampleTrials(tt); 
    licks = x1+find(rez.lickOnTfBin{t, 1}); 
    xend = x1+length(rez.dffsOnTfItpRs{t, 1})-1; 

    plot(x1:x1+length(rez.dffsOnTfItpRs{t, 1})-1, rez.dffsOnTfItpRs{t, 1}, 'Color', rs_color, 'LineWidth', 1); 
    vertline(licks, [0 1], 'lineWidth', 1, 'alpha', 0.5); 
    vertline(xend-length(rez.dffsOnTfItpRs{t, 1})/2, ylim, 'lineWidth', 0.75, 'lineStyle', ':', 'alpha', 0.5); 
    
    x1 = x1 + length(rez.dffsOnTfItpRs{t, 1}) + 1; 
end
title("dff and lick traces");
set(gca, 'TickDir', 'out'); 
axis tight
print(fullfile(filePath, 'Figure', 'dff_lickbouts_rs_example5trials.pdf'), '-dpdf', '-vector', '-bestfit')

%% Stim onset aligned activity (Rs)
% align dffs take the mean and sem
% primary motor (stim onset)
[stimOnDff_rs_itp, timepts] = temporalAlignInterp1(rez.stimOnDffC.rs(:, 1), rez.stimOnDffC.rs(:, 2));
[rez.meanStimOnDff.rs, ~, rez.semStimOnDff.rs] = meanstdsem(cell2mat(stimOnDff_rs_itp));

% primary motor (stim onset) go no-go
[rez.meanStimOnDffGng.rs, rez.semStimOnDffGng.rs] = trialGroupMeanSem(stimOnDff_rs_itp, {[tbytDat.rewardTrI], [tbytDat.punishTrI]});
plotMeanSem(rez.meanStimOnDffGng.rs, rez.semStimOnDffGng.rs, timepts, {'Go', 'NoGo'});
title("Go Nogo Rs (cue onset at time=0)");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:4, 'TickDir', 'out', 'YTick', -1:0.5:1); xlim([-2 2]); ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_cueOn_Gng_rs.pdf'), '-dpdf', '-vector', '-bestfit')

% primary motor (stim onset) go no-go (logistic regrersion)
rez.logitGngCueOnRs = logireg(cell2mat(stimOnDff_rs_itp), [tbytDat.rewardTrI]); % logistic regrersion

%% iti licks (Rs)
itiLickDffRs = flatCell(rez.itiLickDff.rs); % unnest the cell arrays
itiLickDffRsTs = cellfun(@(a) linspace(-1, 1, length(a)), itiLickDffRs(:, 2), 'UniformOutput', false); 
[itiLickDffRsItp, itiLickDffRsItpTs] = temporalAlignInterp1(itiLickDffRs (:, 1), itiLickDffRsTs);
[rez.meanItiLickDff.rs, ~, rez.semItiLickDff.rs] = meanstdsem(cell2mat(itiLickDffRsItp));

plotMeanSemColor(rez.meanItiLickDff.rs, rez.semItiLickDff.rs, itiLickDffRsItpTs, rs_color, {'Rs'});
title("lick during ITI");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); ylim([-.3 .3])

%% postStim licks (Rs)
postStimLickDffRs = flatCell(rez.postStimLickDff.rs); % unnest the cell arrays 
postStimLickDffRsTs = cellfun(@(a) linspace(-1, 1, length(a)), postStimLickDffRs(:, 2), 'UniformOutput', false); 
[postStimLickDffRsItp, postStimLickDffRsItpTs] = temporalAlignInterp1(postStimLickDffRs (:, 1), postStimLickDffRsTs);
[rez.meanPostStimLickDff.rs, ~, rez.semPostStimLickDff.rs] = meanstdsem(cell2mat(postStimLickDffRsItp));

plotMeanSemColor(rez.meanPostStimLickDff.rs, rez.semPostStimLickDff.rs, postStimLickDffRsItpTs, rs_color, {'Rs'});
title("lick during postStim");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.1:3); xlim([-1 1]); ylim([-.3 .3])

%% hit first licks (Rs)
hitFstLickDffRs = cellWithNonEmptyColumns(rez.hitDffFirstLick.rs); 
hitFstLickDffRsTs = cellfun(@(a) linspace(-1, 1, length(a)), hitFstLickDffRs(:, 2), 'UniformOutput', false); 
[hitFstLickDffRsItp, hitFstLickDffRsItpTs] = temporalAlignInterp1(hitFstLickDffRs (:, 1), hitFstLickDffRsTs);
[rez.meanHitFstLickDff.rs, ~, rez.semHitFstLickDff.rs] = meanstdsem(cell2mat(hitFstLickDffRsItp));

plotMeanSemColor(rez.meanHitFstLickDff.rs, rez.semHitFstLickDff.rs, hitFstLickDffRsItpTs, rs_color, {'Rs'});
title("Lick (first) Hit trial");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); ylim([-1 1])

%% first water (Rs)
fstWaterDffRs = cellWithNonEmptyColumns(rez.firstWater.rs); 
fstWaterDffRsTs = cellfun(@(a) linspace(-1, 1, length(a)), fstWaterDffRs(:, 2), 'UniformOutput', false); 
[fstWaterDffRsItp, fstWaterDffRsItpTs] = temporalAlignInterp1(fstWaterDffRs (:, 1), fstWaterDffRsTs);
[rez.meanFstWaterDff.rs, ~, rez.semFstWaterDff.rs] = meanstdsem(cell2mat(fstWaterDffRsItp));

plotMeanSemColor(rez.meanFstWaterDff.rs, rez.semFstWaterDff.rs, fstWaterDffRsItpTs, rs_color, {'Rs'});
title("Water (first)");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); ylim([-0.8 0.8]); grid on; 

print(fullfile(filePath, 'Figure', 'dff_firstWater_rs.pdf'), '-dpdf', '-vector', '-bestfit')

%% first air (Rs)
fstAirDffRs = cellWithNonEmptyColumns(rez.firstAirpuff.rs); 
fstAirDffRsTs = cellfun(@(a) linspace(-1, 1, length(a)), fstAirDffRs(:, 2), 'UniformOutput', false); 
[fstAirDffRsItp, fstAirDffRsItpTs] = temporalAlignInterp1(fstAirDffRs (:, 1), fstAirDffRsTs);
[rez.meanFstAirDff.rs, ~, rez.semFstAirDff.rs] = meanstdsem(cell2mat(fstAirDffRsItp));

plotMeanSemColor(rez.meanFstAirDff.rs, rez.semFstAirDff.rs, fstAirDffRsItpTs, rs_color, {'Rs'});
title("Air (first)");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); ylim([-0.8 0.8])

print(fullfile(filePath, 'Figure', 'dff_firstAirpuff_rs.pdf'), '-dpdf', '-vector', '-bestfit')

%% FA first licks (Rs)
faFstLickDffRs = cellWithNonEmptyColumns(rez.faDffFirstLick.rs); 
faFstLickDffRsTs = cellfun(@(a) linspace(-1, 1, length(a)), faFstLickDffRs(:, 2), 'UniformOutput', false); 
[faFstLickDffRsItp, faFstLickDffRsItpTs] = temporalAlignInterp1(faFstLickDffRs (:, 1), faFstLickDffRsTs);
[rez.meanFaFstLickDff.rs, ~, rez.semFaFstLickDff.rs] = meanstdsem(cell2mat(faFstLickDffRsItp));

plotMeanSem(rez.meanFaFstLickDff.rs, rez.semFaFstLickDff.rs, faFstLickDffRsItpTs, {'Rs'});
title("Lick (first) Fa trial");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); ylim([-1 1])

%% Compare licks during Hit vs FA vs ITI vs postStim (Rs)
% Without normalization
plotMeanSem([rez.meanHitFstLickDff.rs; rez.meanFaFstLickDff.rs; rez.meanItiLickDff.rs; rez.meanPostStimLickDff.rs], ...
    [rez.semHitFstLickDff.rs; rez.semFaFstLickDff.rs; rez.semItiLickDff.rs; rez.semPostStimLickDff.rs], ...
    hitFstLickDffRsItpTs, {'Hit', 'FA', 'ITI', 'postStim'})
title("Lick Hit vs FA vs ITI vs postStim")
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); ylim([-.8 1]); grid on
print(fullfile(filePath, 'Figure', 'dff_lick_Hit_Fa_ITI_postStim_rs.pdf'), '-dpdf', '-vector', '-bestfit')

% With normalization
rez.meanHitFstLickDffBaseSub.rs = baseSubNorm(rez.meanHitFstLickDff.rs, hitFstLickDffRsItpTs, [-1 -0.5]); 
rez.meanFaFstLickDffBaseSub.rs = baseSubNorm(rez.meanFaFstLickDff.rs, faFstLickDffRsItpTs, [-1 -0.5]); 
rez.meanItiLickDffBaseSub.rs = baseSubNorm(rez.meanItiLickDff.rs, itiLickDffRsItpTs, [-1 -0.5]); 
rez.meanPostStimLickDffBaseSub.rs = baseSubNorm(rez.meanPostStimLickDff.rs, postStimLickDffRsItpTs, [-1 -0.5]); 

plotMeanSem([rez.meanHitFstLickDffBaseSub.rs; rez.meanFaFstLickDffBaseSub.rs; rez.meanItiLickDffBaseSub.rs; rez.meanPostStimLickDffBaseSub.rs], ...
    [rez.semHitFstLickDff.rs; rez.semFaFstLickDff.rs; rez.semItiLickDff.rs; rez.semPostStimLickDff.rs], ...
    hitFstLickDffRsItpTs, {'Hit', 'FA', 'ITI', 'postStim'})
title("Lick Hit vs FA vs ITI vs postStim (baseSub)")
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); ylim([-.6 1.2]); grid on
print(fullfile(filePath, 'Figure', 'dff_lick_Hit_Fa_ITI_postStim_baseSub_rs.pdf'), '-dpdf', '-vector', '-bestfit')

%% train svm or naive bayes classifier to classify licks at different task epochs
X_hit = cell2mat(hitFstLickDffRsItp); 
X_fa = cell2mat(faFstLickDffRsItp); 
X_postStim = cell2mat(postStimLickDffRsItp); 
X_iti = cell2mat(itiLickDffRsItp); 

Xs = [X_hit; X_fa; X_postStim; X_iti];
y = [ones(size(X_hit, 1), 1)*1; ones(size(X_fa, 1), 1)*2; ones(size(X_postStim, 1), 1)*3; ones(size(X_iti, 1), 1)*3];
rez.lickHitFaRestSvm.rs = multiClass_svm_peth(Xs, y, 10); 
rez.lickHitFaRestNb.rs = multiClass_svm_peth(Xs, y, 10); 

Xs = [X_hit; X_fa];
y = [ones(size(X_hit, 1), 1)*1; ones(size(X_fa, 1), 1)*2];
rez.lickHitFaSvm.rs = multiClass_svm_peth(Xs, y, 10); 
rez.lickHitFaNb.rs = multiClass_svm_peth(Xs, y, 10); 

figure; hold on; 
plot(mean(rez.lickHitFaRestSvm.rs)); 
plot(mean(rez.lickHitFaRestNb.rs)); 

figure; hold on; 
plot(mean(rez.lickHitFaSvm.rs)); 
plot(mean(rez.lickHitFaNb.rs)); 

%% Miss CueOn (Rs)
missCueDffRs = cellWithNonEmptyColumns(rez.missCueDff.rs); 
missCueDffRsTs = cellfun(@(a) linspace(-1, 1, length(a)), missCueDffRs(:, 2), 'UniformOutput', false); 
[missCueDffRsItp, missCueDffRsItpTs] = temporalAlignInterp1(missCueDffRs (:, 1), missCueDffRsTs);
[rez.meanMirsCueDff.rs, ~, rez.semMirsCueDff.rs] = meanstdsem(cell2mat(missCueDffRsItp));

plotMeanSem(rez.meanMirsCueDff.rs, rez.semMirsCueDff.rs, missCueDffRsItpTs, {'Rs'});
title("Mirs CueOn trial");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); ylim([-1 1])

print(fullfile(filePath, 'Figure', 'dff_miss_cueOn_rs.pdf'), '-dpdf', '-vector', '-bestfit')

%% CR CueOn (Rs)
crCueDffRs = cellWithNonEmptyColumns(rez.crCueDff.rs); 
crCueDffRsTs = cellfun(@(a) linspace(-1, 1, length(a)), crCueDffRs(:, 2), 'UniformOutput', false); 
[crCueDffRsItp, crCueDffRsItpTs] = temporalAlignInterp1(crCueDffRs (:, 1), crCueDffRsTs);
[rez.meanCrCueDff.rs, ~, rez.semCrCueDff.rs] = meanstdsem(cell2mat(crCueDffRsItp));

plotMeanSem(rez.meanCrCueDff.rs, rez.semCrCueDff.rs, crCueDffRsItpTs, {'Rs'});
title("Correct rejection CueOn trial");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:.5:4, 'TickDir', 'out', 'YTick', -3:0.2:3); xlim([-1 1]); ylim([-1 1])

print(fullfile(filePath, 'Figure', 'dff_cr_cueOn_rs.pdf'), '-dpdf', '-vector', '-bestfit')

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
        % 'tbytPsthC' contains trial-by-trial data in each cell (it is arsumed that data are temporally aligned acrors trials already).
        % 'groupC' contains logical for each grouping, logicals are arsumed to be of same lengths as the number of trials.
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

        % Determine the number of samples for each clars in the test set
        numClars0 = sum(Y == 0);
        numClars1 = sum(Y == 1);
        testSetSize = min(numClars0, numClars1) / 2; % Or any other criterion

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


