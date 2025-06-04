function dffPostprocessAuditoryGng_averaging_m1_dual_6OHDA(filePath, channel)
%filePath = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA008/DA008_101823';
% This function performs further analyses and plotting using data in
% 'rez'(outcome from 'dffPostprocessGng_averaging.m').
%

%% load rez tbytDat
match_header = regexp(filePath, 'm\d{1,4}_\d{6}', 'match');
header = match_header{1};
load(fullfile(filePath, 'Matfiles', strcat(header, '_', channel, '_dff_evtAligned_regionMask.mat')), 'rez', 'voluntaryHitLickI');

% user variables
m1_color = [0 240 240]./255;

if exist(fullfile(filePath, 'Figure'), 'dir')~=7
    mkdir(fullfile(filePath, 'Figure'))
end

%% Stim onset aligned activity (m1)
% align dffs take the mean and sem
% m1 stim onset all trials
[stimOn_m1L_itp, ts_ext] = temporalAlignInterp1(rez.stimOnDffC.m1L(:, 1), rez.stimOnDffC.m1L(:, 2), 0.01);
[rezM1.mStimOnL, ~, rezM1.eStimOnL] = meanstdsem(cell2mat(stimOn_m1L_itp));

[stimOn_m1R_itp, ~] = temporalAlignInterp1(rez.stimOnDffC.m1R(:, 1), rez.stimOnDffC.m1R(:, 2), 0.01);
[rezM1.mStimOnDffR, ~, rezM1.eStimOnDffR] = meanstdsem(cell2mat(stimOn_m1R_itp));

% m1 stim onset voluntary lick trials
[stimOnVolLick_m1L_itp, ~] = temporalAlignInterp1(rez.stimOnDffC.m1L(voluntaryHitLickI, 1), rez.stimOnDffC.m1L(voluntaryHitLickI, 2), 0.01);
[rezM1.mStimOnVolLickDffL, ~, rezM1.eStimOnVolLickDffL] = meanstdsem(cell2mat(stimOnVolLick_m1L_itp));

[stimOnVolLick_m1R_itp, ~] = temporalAlignInterp1(rez.stimOnDffC.m1R(voluntaryHitLickI, 1), rez.stimOnDffC.m1R(voluntaryHitLickI, 2), 0.01);
[rezM1.mStimOnVolLickDffR, ~, rezM1.eStimOnVolLickDffR] = meanstdsem(cell2mat(stimOnVolLick_m1R_itp));

%% Hit aligned activity (m1)
% m1 hit all trials 
[hitLickFst_m1L_itp, ts] = temporalAlignInterp1(rez.hitLickFst.m1L(:, 1), rez.hitLickFst.m1L(:, 2), 0.01);
[rezM1.meanHitLickFstL, ~, rezM1.semHitLickFstL] = meanstdsem(cell2mat(hitLickFst_m1L_itp));

[hitLickFst_m1R_itp, ~] = temporalAlignInterp1(rez.hitLickFst.m1R(:, 1), rez.hitLickFst.m1R(:, 2), 0.01);
[rezM1.meanHitLickFstR, ~, rezM1.semHitLickFstR] = meanstdsem(cell2mat(hitLickFst_m1R_itp));

% m1 hit voluntary lick trials
[hitVolLickFst_m1L_itp, ~] = temporalAlignInterp1(rez.hitLickFst.m1L(voluntaryHitLickI, 1), rez.hitLickFst.m1L(voluntaryHitLickI, 2), 0.01);
[rezM1.meanHitVolLickFstL, ~, rezM1.semHitVolLickFstL] = meanstdsem(cell2mat(hitVolLickFst_m1L_itp));

[hitVolLickFst_m1R_itp, ~] = temporalAlignInterp1(rez.hitLickFst.m1R(voluntaryHitLickI, 1), rez.hitLickFst.m1R(voluntaryHitLickI, 2), 0.01);
[rezM1.meanHitVolLickFstR, ~, rezM1.semHitVolLickFstR] = meanstdsem(cell2mat(hitVolLickFst_m1R_itp));

%% Water aligned activity (m1)
% m1 water all trials
[water_m1L_itp, ~] = temporalAlignInterp1(rez.firstWater.m1L(:, 1), rez.firstWater.m1L(:, 2), 0.01);
[rezM1.mWaterL, ~, rezM1.eWaterL] = meanstdsem(cell2mat(water_m1L_itp));

[water_m1R_itp, ~] = temporalAlignInterp1(rez.firstWater.m1R(:, 1), rez.firstWater.m1R(:, 2), 0.01);
[rezM1.mWaterR, ~, rezM1.eWaterR] = meanstdsem(cell2mat(water_m1R_itp));

% m1 water voluntary lick trials 
[waterVolLick_m1L_itp, ~] = temporalAlignInterp1(rez.firstWater.m1L(voluntaryHitLickI, 1), rez.firstWater.m1L(voluntaryHitLickI, 2), 0.01);
[rezM1.mWaterVolLickL, ~, rezM1.eWaterVolLickL] = meanstdsem(cell2mat(waterVolLick_m1L_itp));

[waterVolLick_m1R_itp, ~] = temporalAlignInterp1(rez.firstWater.m1R(voluntaryHitLickI, 1), rez.firstWater.m1R(voluntaryHitLickI, 2), 0.01);
[rezM1.mWaterVolLickR, ~, rezM1.eWaterVolLickR] = meanstdsem(cell2mat(waterVolLick_m1R_itp));

%% save rez
save(fullfile(filePath, 'Matfiles', strcat(header, '_dff_evtAligned_regionMask_m1_', channel, '.mat')), 'rezM1');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
end


