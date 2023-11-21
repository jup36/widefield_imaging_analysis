function dffPostprocessGng_averaging(filePath)
%filePath = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA008/DA008_101823';
%tbytDat
%   evtType:
%       1 or 2: visual (common, uncommon)
%           evtOn: photoDiodeOn
%           evtOff: photoDiodeOff
%           periEvtWin: evtOn-1:evtOff+1 (e.g. 4s)
%       3: reward
%           evtOn: water delivery
%           evtOff: 4-s after water delivery
%           periEvtWin: evtOn-1:evtOff (e.g. 5s)
%           Note that 5-s peri-reward window was used -1 to 4s relative to reward

% load trial-by-trial behavior & task data tbytDat.mat
%parentDir = fileparts(filePath);
[~, header] = fileparts(filePath);
fileBeh = GrabFiles_sort_trials('tbytDat_dff', 0, {fullfile(filePath, 'Matfiles')});
if isempty(fileBeh)
    fileBeh = GrabFiles_sort_trials('tbytDat_dff', 1, {filePath});
end
load(fullfile(fileBeh{1}), 'tbytDat')

% load the preprocessed dffs collect
load(fullfile(filePath, 'Matfiles', strcat(header, '_dffsmCollect.mat')), 'dffsmCell');

% load Allen dorsalMap
load('allenDorsalMap', 'dorsalMaps');

% load Allen transformation object
folder_list_imgTrial = GrabFiles_sort_trials(header, 0, {fullfile(filePath, strcat(header, '_img'))}); % use GrabFiles_sort_trials to sort both files and folders
load(fullfile(folder_list_imgTrial{1}, 'transParamsAllen.mat'), 'transParams')

% file directory for trials
filePathTrials = fullfile(filePath, strcat(header, '_trials'));
if exist(filePathTrials, 'dir') == 0
    mkdir(filePathTrials);
end

% apply region mask
% M1
motorId = cell2mat(cellfun(@(a) strcmpi(a, 'Primary motor area, Layer 1'), dorsalMaps.labelTable.name, 'UniformOutput', false));
motorI = dorsalMaps.labelTable.id(motorId);
motorMask = dorsalMaps.dorsalMap==motorI;

% M2
smotorId = cell2mat(cellfun(@(a) strcmpi(a, 'Secondary motor area, Layer 1'), dorsalMaps.labelTable.name, 'UniformOutput', false));
smotorI = dorsalMaps.labelTable.id(smotorId);
smotorMask = dorsalMaps.dorsalMap==smotorI;

% S1
ssMask = getSSmask(dorsalMaps);

% Visual area
vMask = getVmask(dorsalMaps);

% Retrosplenial area
rsMask = getRSmask(dorsalMaps);

rwdDffC_motor = cell(length(tbytDat), 4);
rwdDffC_smotor = cell(length(tbytDat), 4);
rwdDffC_ss = cell(length(tbytDat), 4);
rwdDffC_v = cell(length(tbytDat), 4);
rwdDffC_rs = cell(length(tbytDat), 4);


% trial index
% masks
% time index (baseline and event times)



for t = 1:length(tbytDat)

    % align image stack to the AllenCCF
    wfa = alignStackToAllenKabsch(dffsmCell{t}, dorsalMaps.dorsalMap, transParams.tformObj);

    timeStamps = tbytDat(t).frameTrel; % store timestamps

    % Align to stim Onset
    [stimOnDffC.m1{t, 1}, stimOnDffC.m1{t, 2}] = alignToEvent(wfa, 0, tbytDat(t).frameTrel, motorMask, [-3 2.5]);
    [stimOnDffC.m2{t, 1}, stimOnDffC.m2{t, 2}] = alignToEvent(wfa, 0, tbytDat(t).frameTrel, smotorMask, [-3 2.5]);
    [stimOnDffC.s1{t, 1}, stimOnDffC.s1{t, 2}] = alignToEvent(wfa, 0, tbytDat(t).frameTrel, ssMask, [-3 2.5]);
    [stimOnDffC.v{t, 1}, stimOnDffC.v{t, 2}] = alignToEvent(wfa, 0, tbytDat(t).frameTrel, vMask, [-3 2.5]);
    [stimOnDffC.rs{t, 1}, stimOnDffC.rs{t, 2}] = alignToEvent(wfa, 0, tbytDat(t).frameTrel, rsMask, [-3 2.5]);

    fprintf("processed trial #%d\n", t)
end


%% Stim onset aligned activity (M1)
% align dffs take the mean and sem
% primary motor (stim onset)
[stimOnDff_m1_itp, timepts] = temporalAlignInterp1(stimOnDffC.m1(:, 1), stimOnDffC.m1(:, 2));
[rez.mStimOnDffM1, ~, rez.semStimOnDffM1] = meanstdsem(cell2mat(stimOnDff_m1_itp));

% primary motor (stim onset) go no-go
[rez.mStimOnDffGngM1, rez.sStimOnDffGngM1] = trialGroupMeanSem(stimOnDff_m1_itp, {[tbytDat.rewardTrI], [tbytDat.punishTrI]});
plotMeanSem(rez.mStimOnDffGngM1, rez.sStimOnDffGngM1, timepts, {'Go', 'NoGo'});
title("Go Nogo M1 (cue onset at time=0)");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:4); xlim([-2 2]); ylim([-0.6 1])
print(fullfile(filePath, 'Figure', 'dff_cueOn_Gng_m1.pdf'), '-dpdf', '-vector', '-bestfit')

% primary motor (stim onset) go no-go (logistic regression)
rez.logitGngCueOnM1 = logireg(cell2mat(stimOnDff_m1_itp), [tbytDat.rewardTrI]); % logistic regression

%% Stim onset aligned activity (M2)
% supplementary motor (stim onset)
stimOnDff_m2_itp = temporalAlignInterp1(stimOnDffC.m2(:, 1), stimOnDffC.m2(:, 2));
[rez.mStimOnDffM2, ~, rez.semStimOnDffM2] = meanstdsem(cell2mat(stimOnDff_m2_itp));

% supplementary motor (stim onset) go no-go
[rez.mStimOnDffGngM2, rez.sStimOnDffGngM2] = trialGroupMeanSem(stimOnDff_m2_itp, {[tbytDat.rewardTrI], [tbytDat.punishTrI]});
plotMeanSem(rez.mStimOnDffGngM2, rez.sStimOnDffGngM2, timepts, {'Go', 'NoGo'});
title("Go Nogo M2 (cue onset at time=0)");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:4); xlim([-2 2]); ylim([-0.6 1]);
print(fullfile(filePath, 'Figure', 'dff_cueOn_Gng_m2.pdf'), '-dpdf', '-vector', '-bestfit')

% supplementary motor (stim onset) go no-go (logistic regression)
rez.logitGngCueOnM2 = logireg(cell2mat(stimOnDff_m2_itp), [tbytDat.rewardTrI]); % logistic regression

%% Stim onset aligned activity (S1)
% somatosensory (stim onset)
stimOnDff_s1_itp = temporalAlignInterp1(stimOnDffC.s1(:, 1), stimOnDffC.s1(:, 2));
[rez.mStimOnDffS1, ~, rez.semStimOnDffS1] = meanstdsem(cell2mat(stimOnDff_s1_itp));

% somato sensory (stim onset) go no-go
[rez.mStimOnDffGngS1, rez.sStimOnDffGngS1] = trialGroupMeanSem(stimOnDff_s1_itp, {[tbytDat.rewardTrI], [tbytDat.punishTrI]});
plotMeanSem(rez.mStimOnDffGngS1, rez.sStimOnDffGngS1, timepts, {'Go', 'NoGo'});
title("Go Nogo S1 (cue onset at time=0)");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:4); xlim([-2 2]); ylim([-0.6 1]);
print(fullfile(filePath, 'Figure', 'dff_cueOn_Gng_s1.pdf'), '-dpdf', '-vector', '-bestfit')

% somato sensory (stim onset) go no-go (logistic regression)
rez.logitGngCueOnS1 = logireg(cell2mat(stimOnDff_s1_itp), [tbytDat.rewardTrI]); % logistic regression

%% Stim onset aligned activity (V)
% visual (stim onset)
stimOnDff_v_itp = temporalAlignInterp1(stimOnDffC.v(:, 1), stimOnDffC.v(:, 2));
[rez.mStimOnDffV, ~, rez.semStimOnDffV] = meanstdsem(cell2mat(stimOnDff_v_itp));

% visual (stim onset) go no-go
[rez.mStimOnDffGngV, rez.sStimOnDffGngV] = trialGroupMeanSem(stimOnDff_v_itp, {[tbytDat.rewardTrI], [tbytDat.punishTrI]});
plotMeanSem(rez.mStimOnDffGngV, rez.sStimOnDffGngV, timepts, {'Go', 'NoGo'});
title("Go Nogo V (cue onset at time=0)");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:4); xlim([-2 2]); ylim([-0.6 1]);
print(fullfile(filePath, 'Figure', 'dff_cueOn_Gng_v.pdf'), '-dpdf', '-vector', '-bestfit')

% visual (stim onset) go no-go (logistic regression)
rez.logitGngCueOnV = logireg(cell2mat(stimOnDff_v_itp), [tbytDat.rewardTrI]); % logistic regression

%% retrosplenial cortex (stim onset)
stimOnDff_rs_itp = temporalAlignInterp1(stimOnDffC.rs(:, 1), stimOnDffC.rs(:, 2));
[rez.mStimOnDffRs, ~, rez.semStimOnDffRs] = meanstdsem(cell2mat(stimOnDff_rs_itp));

% retrosplenial (stim onset) go no-go
[rez.mStimOnDffGngRs, rez.sStimOnDffGngRs] = trialGroupMeanSem(stimOnDff_rs_itp, {[tbytDat.rewardTrI], [tbytDat.punishTrI]});
plotMeanSem(rez.mStimOnDffGngRs, rez.sStimOnDffGngRs, timepts, {'Go', 'NoGo'});
title("Go Nogo Rs (cue onset at time=0)");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:4); xlim([-2 2]); ylim([-0.6 1]);
print(fullfile(filePath, 'Figure', 'dff_cueOn_Gng_rs.pdf'), '-dpdf', '-vector', '-bestfit')

% retrosplenial (stim onset) go no-go (logistic regression)
rez.logitGngCueOnRs = logireg(cell2mat(stimOnDff_rs_itp), [tbytDat.rewardTrI]); % logistic regression

%% plot stim onset PETHs together
plotMeanSem([rez.mStimOnDffM1; rez.mStimOnDffM2; rez.mStimOnDffS1; rez.mStimOnDffRs; rez.mStimOnDffV], ...
    [rez.semStimOnDffM1; rez.semStimOnDffM2; rez.semStimOnDffS1; rez.semStimOnDffRs; rez.semStimOnDffV], ...
    timepts, ...
    {'M1', 'M2', 'S1', 'RS', 'V'});
xlabel('Time (s)')
ylabel('DFF')
set(gca, 'XTick', -2:1:4)
xlim([-2 2])
print(fullfile(filePath, 'Figure', 'dff_cueOn_Collect.pdf'), '-dpdf', '-vector', '-bestfit')

%% plot classification accuracy
plotMeanSem(smooth2a([rez.logitGngCueOnM1; rez.logitGngCueOnM2; rez.logitGngCueOnS1; rez.logitGngCueOnRs; rez.logitGngCueOnV], 0, 1), ...
    zeros(5, length(timepts)), timepts, {'M1', 'M2', 'S1', 'RS', 'V'});
title("Cross-validated classification accuracy")
xlabel('Time (s)')
ylabel('Classification Accuracy')
set(gca, 'XTick', -2:1:4)
xlim([-2 2])
print(fullfile(filePath, 'Figure', 'logisticRegression_Gng_cueOn_Collect.pdf'), '-dpdf', '-vector', '-bestfit')






% baseline dff image
rez.meanBaseDffImgMotor = nanmean(cell2mat(reshape(rwdDffC_motor(:, 3), [1, 1, length(tbytDat)])), 3);
imageFrameWithNaNsEdgeMap(rez.meanBaseDffImgMotor, [-.3 .3], dorsalMaps.edgeMap, [0, 0, 1], 0.3)
title('mean baseline M1')
set(gca, 'XTickLabel', {})
set(gca, 'YTickLabel', {})
print(fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA003/DA003_083023/Figure', 'dff_base_M1_Img.pdf'), '-dpdf', '-vector', '-bestfit')

% reward dff image
rez.meanRwdDffImgMotor = nanmean(cell2mat(reshape(rwdDffC_motor(:, 4), [1, 1, length(tbytDat)])), 3);
imageFrameWithNaNsEdgeMap(rez.meanRwdDffImgMotor, [-.3 .3], dorsalMaps.edgeMap, [0, 0, 1], 0.3)
print(fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA003/DA003_083023/Figure', 'dff_reward_M1_Img.pdf'), '-dpdf', '-vector', '-bestfit')
title('mean reward M1')
set(gca, 'XTickLabel', {})
set(gca, 'YTickLabel', {})
print(fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA003/DA003_083023/Figure', 'dff_reward_M1_Img.pdf'), '-dpdf', '-vector', '-bestfit')



showAllenEdgesOnTransformedWF(rez.meanBaseDffImgMotor, dorsalMaps, 0, 1)

imageFrameWithNaNs(rez.meanRwdDffImgMotor, [-.3 .3])



% save tbytDat without dffs
tbytDat = rmfield(tbytDat, {'dff', 'dffsm'});
save(fullfile(parentDir, Matfiles, [header, '_tbytDat_dff']), 'tbytDat')



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [evtAlignedWfa, evtAlignedTs] = alignToEvent(wfa, eventTimeToAlign, frameT, mask, timeWin)
        wfaMask = apply2DMaskTo3DStack(wfa, mask);

        if eventTimeToAlign >= min(timeWin) && eventTimeToAlign <= max(timeWin)
            frameI = frameT >= min(timeWin) & frameT <= max(timeWin);
            evtAlignedWfa = wfaMask(frameI, 1);
            evtAlignedTs = frameT(frameI);
        else
            error("The event(s) to align is out of the time window!")
        end
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
        % X: observations-by-feature
        %   e.g., cell2mat(stimOnDff_m1_itp);
        % Y: label
        %   e.g., [tbytDat.rewardTrI];
        Y = Y(:);

        % Assuming X is your feature matrix and Y is your labels
        [trainInd, testInd] = crossvalind('HoldOut', Y, 0.3); % 70-30 split

        accuracy = nan(1, size(X, 2));
        % Training
        for jj = 1:size(X, 2) % iterate features
            model = fitglm(X(trainInd, jj), Y(trainInd), 'Distribution', 'binomial');
            % Testing
            Y_pred = predict(model, X(testInd,jj)) > 0.5; % Thresholding at 0.5
            accuracy(1, jj) = mean(Y_pred == Y(testInd));
        end

    end



end


