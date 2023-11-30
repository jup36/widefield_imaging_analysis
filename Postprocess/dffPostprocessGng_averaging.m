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

%% load trial-by-trial behavior & task data tbytDat.mat
%parentDir = fileparts(filePath);
[~, header] = fileparts(filePath);
fileBeh = GrabFiles_sort_trials('tbytDat_dff', 0, {fullfile(filePath, 'Matfiles')});
if isempty(fileBeh)
    fileBeh = GrabFiles_sort_trials('tbytDat_dff', 1, {filePath});
end
load(fullfile(fileBeh{1}), 'tbytDat')

tbytDat = parseGngTrials(tbytDat);

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

for t = 1:length(tbytDat)

    % align image stack to the AllenCCF
    dff = alignStackToAllenKabsch(dffsmCell{t}, dorsalMaps.dorsalMap, transParams.tformObj);

    frT = tbytDat(t).frameTrel; % store timestamps

    % apply region masks
    dffM1 = apply2DMaskTo3DStack(dff, motorMask);
    dffM2 = apply2DMaskTo3DStack(dff, smotorMask);
    dffSs = apply2DMaskTo3DStack(dff, ssMask);
    dffV1 = apply2DMaskTo3DStack(dff, vMask);
    dffRs = apply2DMaskTo3DStack(dff, rsMask);

    %% Align to stim Onset
    [rez.stimOnDffC.m1{t, 1}, rez.stimOnDffC.m1{t, 2}] = alignToEvent(dffM1, 0, frT, [-3 2.5]);
    [rez.stimOnDffC.m2{t, 1}, rez.stimOnDffC.m2{t, 2}] = alignToEvent(dffM2, 0, frT, [-3 2.5]);
    [rez.stimOnDffC.Ss{t, 1}, rez.stimOnDffC.s1{t, 2}] = alignToEvent(dffSs, 0, frT, [-3 2.5]);
    [rez.stimOnDffC.v1{t, 1}, rez.stimOnDffC.v1{t, 2}] = alignToEvent(dffV1, 0, frT, [-3 2.5]);
    [rez.stimOnDffC.rs{t, 1}, rez.stimOnDffC.rs{t, 2}] = alignToEvent(dffRs, 0, frT, [-3 2.5]);

    %% Align to licks during ITI
    for ii = 1:length(tbytDat(t).itiLickChunk) % there can be multiple bouts
        [rez.itiLickDff.m1{t, 1}{ii, 1}, rez.itiLickDff.m1{t, 1}{ii, 2}] = alignToEvent(dffM1, tbytDat(t).itiLickChunk{ii}(1), frT, [-1 1]); % M1 align to the ITI lick bout use the 1st lick of each bout
        [rez.itiLickDff.m2{t, 1}{ii, 1}, rez.itiLickDff.m2{t, 1}{ii, 2}]  = alignToEvent(dffM2, tbytDat(t).itiLickChunk{ii}(1), frT, [-1 1]); % M2 align to the ITI lick bout use the 1st lick of each bout
        [rez.itiLickDff.ss{t, 1}{ii, 1}, rez.itiLickDff.ss{t, 1}{ii, 2}] = alignToEvent(dffSs, tbytDat(t).itiLickChunk{ii}(1), frT, [-1 1]); % Ss align to the ITI lick bout use the 1st lick of each bout
        [rez.itiLickDff.v1{t, 1}{ii, 1}, rez.itiLickDff.v1{t, 1}{ii, 2}] = alignToEvent(dffV1, tbytDat(t).itiLickChunk{ii}(1), frT, [-1 1]); % V1 align to the ITI lick bout use the 1st lick of each bout
        [rez.itiLickDff.rs{t, 1}{ii, 1}, rez.itiLickDff.rs{t, 1}{ii, 2}] = alignToEvent(dffRs, tbytDat(t).itiLickChunk{ii}(1), frT, [-1 1]); % Rs align to the ITI lick bout use the 1st lick of each bout
    end

    
    %% Align to water (1st water only)
    if ~isempty(tbytDat(t).water)
        firstWater = tbytDat(t).water(1) - tbytDat(t).stimOn; 
        % m1
        [rez.firstWater.m1{t, 1}, rez.firstWater.m1{t, 2}] = alignToEvent(dffM1, firstWater, frT, [-1 1]); % consume lick 1st in the bout
        % m2
        [rez.firstWater.m2{t, 1}, rez.firstWater.m2{t, 2}] = alignToEvent(dffM2, firstWater, frT, [-1 1]); % consume lick 1st in the bout
        % ss
        [rez.firstWater.ss{t, 1}, rez.firstWater.ss{t, 2}] = alignToEvent(dffSs, firstWater, frT, [-1 1]); % consume lick 1st in the bout
        % v1
        [rez.firstWater.v1{t, 1}, rez.firstWater.v1{t, 2}] = alignToEvent(dffV1, firstWater, frT, [-1 1]); % consume lick 1st in the bout
        % rs
        [rez.firstWater.rs{t, 1}, rez.firstWater.rs{t, 2}]  = alignToEvent(dffRs, firstWater, frT, [-1 1]); % consume lick 1st in the bout
    end

    %% Align to airpuff (1st airpuff only)
    if ~isempty(tbytDat(t).airpuff)
        firstAirpuff = tbytDat(t).airpuff(1) - tbytDat(t).stimOn; 
        % m1
        [rez.firstAirpuff.m1{t, 1}, rez.firstAirpuff.m1{t, 2}] = alignToEvent(dffM1, firstAirpuff, frT, [-1 1]); % consume lick 1st in the bout
        % m2
        [rez.firstAirpuff.m2{t, 1}, rez.firstAirpuff.m2{t, 2}] = alignToEvent(dffM2, firstAirpuff, frT, [-1 1]); % consume lick 1st in the bout
        % ss
        [rez.firstAirpuff.ss{t, 1}, rez.firstAirpuff.ss{t, 2}] = alignToEvent(dffSs, firstAirpuff, frT, [-1 1]); % consume lick 1st in the bout
        % v1
        [rez.firstAirpuff.v1{t, 1}, rez.firstAirpuff.v1{t, 2}] = alignToEvent(dffV1, firstAirpuff, frT, [-1 1]); % consume lick 1st in the bout
        % rs
        [rez.firstAirpuff.rs{t, 1}, rez.firstAirpuff.rs{t, 2}]  = alignToEvent(dffRs, firstAirpuff, frT, [-1 1]); % consume lick 1st in the bout
    end

    %% Align to licks during Cue
    % Hit
    % ToDo: Align to the first or last lick? 
    if ~isempty(tbytDat(t).hitLicks)
        % m1
        [rez.hitDffFirstLick.m1{t, 1}, rez.hitDffFirstLick.m1{t, 2}] = alignToEvent(dffM1, tbytDat(t).hitLicks(1), frT, [-1 1]); % hit lick 1st in the bout
        [rez.hitDffLastLick.m1{t, 1}, rez.hitDffLastLick.m1{t, 2}] = alignToEvent(dffM1, tbytDat(t).hitLicks(end), frT, [-1 1]); % hit lick last in the bout
        [rez.hitDffCueOn.m1{t, 1}, rez.hitDffCueOn.m1{t, 2}] = alignToEvent(dffM1, 0, frT, [-1 1]); % hit lick cueOn
        % m2
        [rez.hitDffFirstLick.m2{t, 1}, rez.hitDffFirstLick.m2{t, 2}] = alignToEvent(dffM2, tbytDat(t).hitLicks(1), frT, [-1 1]);
        [rez.hitDffLastLick.m2{t, 1}, rez.hitDffLastLick.m2{t, 2}] = alignToEvent(dffM2, tbytDat(t).hitLicks(end), frT, [-1 1]);
        [rez.hitDffCueOn.m2{t, 1}, rez.hitDffCueOn.m2{t, 2}] = alignToEvent(dffM2, 0, frT, [-1 1]);
        % ss
        [rez.hitDffFirstLick.ss{t, 1}, rez.hitDffFirstLick.ss{t, 2}] = alignToEvent(dffSs, tbytDat(t).hitLicks(1), frT, [-1 1]);
        [rez.hitDffLastLick.ss{t, 1}, rez.hitDffLastLick.ss{t, 2}] = alignToEvent(dffSs, tbytDat(t).hitLicks(end), frT, [-1 1]);
        [rez.hitDffCueOn.ss{t, 1}, rez.hitDffCueOn.ss{t, 2}] = alignToEvent(dffSs, 0, frT, [-1 1]);
        % v1
        [rez.hitDffFirstLick.v1{t, 1}, rez.hitDffFirstLick.v1{t, 2}] = alignToEvent(dffV1, tbytDat(t).hitLicks(1), frT, [-1 1]);
        [rez.hitDffLastLick.v1{t, 1}, rez.hitDffLastLick.v1{t, 2}] = alignToEvent(dffV1, tbytDat(t).hitLicks(end), frT, [-1 1]);
        [rez.hitDffCueOn.v1{t, 1}, rez.hitDffCueOn.v1{t, 2}] = alignToEvent(dffV1, 0, frT, [-1 1]);
        % rs
        [rez.hitDffFirstLick.rs{t, 1}, rez.hitDffFirstLick.rs{t, 2}] = alignToEvent(dffRs, tbytDat(t).hitLicks(1), frT, [-1 1]);
        [rez.hitDffLastLick.rs{t, 1}, rez.hitDffFirstLick.rs{t, 2}] = alignToEvent(dffRs, tbytDat(t).hitLicks(end), frT, [-1 1]);
        [rez.hitDffCueOn.rs{t, 1}, rez.hitDffCueOn.rs{t, 2}] = alignToEvent(dffRs, 0, frT, [-1 1]);
                
        if ~isempty(tbytDat(t).consumeLicks)
            % m1
            [rez.hitDffConsumeLick.m1{t, 1}, rez.hitDffConsumeLick.m1{t, 2}] = alignToEvent(dffM1, tbytDat(t).consumeLicks(1), frT, [-1 1]); % consume lick 1st in the bout
            % m2
            [rez.hitDffConsumeLick.m2{t, 1}, rez.hitDffConsumeLick.m2{t, 2}] = alignToEvent(dffM2, tbytDat(t).consumeLicks(1), frT, [-1 1]); % consume lick 1st in the bout
            % ss
            [rez.hitDffConsumeLick.ss{t, 1}, rez.hitDffConsumeLick.ss{t, 2}] = alignToEvent(dffSs, tbytDat(t).consumeLicks(1), frT, [-1 1]); % consume lick 1st in the bout
            % v1
            [rez.hitDffConsumeLick.v1{t, 1}, rez.hitDffConsumeLick.v1{t, 2}] = alignToEvent(dffV1, tbytDat(t).consumeLicks(1), frT, [-1 1]); % consume lick 1st in the bout
            % rs
            [rez.hitDffConsumeLick.rs{t, 1}, rez.hitDffConsumeLick.rs{t, 2}]  = alignToEvent(dffRs, tbytDat(t).consumeLicks(1), frT, [-1 1]); % consume lick 1st in the bout
        end

    % Miss
    elseif tbytDat(t).rewardTrI && isempty(tbytDat(t).water)
        % m1
        [rez.missCueDff.m1{t, 1}, rez.missCueDff.m1{t, 2}] = alignToEvent(dffM1, 0, frT, [-1 1]); % align to cue onset
        % m2
        [rez.missCueDff.m2{t, 1}, rez.missCueDff.m2{t, 2}] = alignToEvent(dffM2, 0, frT, [-1 1]); % align to cue onset
        % ss
        [rez.missCueDff.ss{t, 1}, rez.missCueDff.ss{t, 2}] = alignToEvent(dffSs, 0, frT, [-1 1]); % align to cue onset
        % v1
        [rez.missCueDff.v1{t, 1}, rez.missCueDff.v1{t, 2}] = alignToEvent(dffV1, 0, frT, [-1 1]); % align to cue onset
        % rs
        [rez.missCueDff.rs{t, 1}, rez.missCueDff.rs{t, 2}] = alignToEvent(dffRs, 0, frT, [-1 1]); % align to cue onset
   
    % False Alarm
    elseif ~isempty(tbytDat(t).faLicks)
        % m1
        [rez.faDffFirstLick.m1{t, 1}, rez.faDffFirstLick.m1{t, 2}] = alignToEvent(dffM1, tbytDat(t).faLicks(1), frT, [-1 1]); % fa lick 1st in the bout
        [rez.faDffLastLick.m1{t, 1}, rez.faDffLastLick.m1{t, 2}] = alignToEvent(dffM1, tbytDat(t).faLicks(end), frT, [-1 1]); % fa lick last in the bout
        [rez.faDffCueOn.m1{t, 1}, rez.faDffCueOn.m1{t, 2}] = alignToEvent(dffM1, 0, frT, [-1 1]); % fa lick cueOn
        % m2
        [rez.faDffFirstLick.m2{t, 1}, rez.faDffFirstLick.m2{t, 2}] = alignToEvent(dffM2, tbytDat(t).faLicks(1), frT, [-1 1]);
        [rez.faDffLastLick.m2{t, 1}, rez.faDffLastLick.m2{t, 2}] = alignToEvent(dffM2, tbytDat(t).faLicks(end), frT, [-1 1]);
        [rez.faDffCueOn.m2{t, 1}, rez.faDffCueOn.m2{t, 2}] = alignToEvent(dffM2, 0, frT, [-1 1]);
        % ss
        [rez.faDffFirstLick.ss{t, 1}, rez.faDffFirstLick.ss{t, 2}] = alignToEvent(dffSs, tbytDat(t).faLicks(1), frT, [-1 1]);
        [rez.faDffLastLick.ss{t, 1}, rez.faDffLastLick.ss{t, 2}] = alignToEvent(dffSs, tbytDat(t).faLicks(end), frT, [-1 1]);
        [rez.faDffCueOn.ss{t, 1}, rez.faDffCueOn.ss{t, 2}] = alignToEvent(dffSs, 0, frT, [-1 1]);
        % v1
        [rez.faDffFirstLick.v1{t, 1}, rez.faDffFirstLick.v1{t, 2}] = alignToEvent(dffV1, tbytDat(t).faLicks(1), frT, [-1 1]);
        [rez.faDffLastLick.v1{t, 1}, rez.faDffLastLick.v1{t, 2}] = alignToEvent(dffV1, tbytDat(t).faLicks(end), frT, [-1 1]);
        [rez.faDffCueOn.v1{t, 1}, rez.faDffCueOn.v1{t, 2}] = alignToEvent(dffV1, 0, frT, [-1 1]);
        % rs
        [rez.faDffFirstLick.rs{t, 1}, rez.faDffFirstLick.rs{t, 2}] = alignToEvent(dffRs, tbytDat(t).faLicks(1), frT, [-1 1]);
        [rez.faDffLastLick.rs{t, 1}, rez.faDffFirstLick.rs{t, 2}] = alignToEvent(dffRs, tbytDat(t).faLicks(end), frT, [-1 1]);
        [rez.faDffCueOn.rs{t, 1}, rez.faDffCueOn.rs{t, 2}] = alignToEvent(dffRs, 0, frT, [-1 1]);
                
        if ~isempty(tbytDat(t).postAirLicks)
            % m1
            [rez.faDffPostAirLick.m1{t, 1}, rez.faDffPostAirLick.m1{t, 2}] = alignToEvent(dffM1, tbytDat(t).postAirLicks(1), frT, [-1 1]); % postAir lick 1st in the bout
            % m2
            [rez.faDffPostAirLick.m2{t, 1}, rez.faDffPostAirLick.m2{t, 2}] = alignToEvent(dffM2, tbytDat(t).postAirLicks(1), frT, [-1 1]); % postAir lick 1st in the bout
            % ss
            [rez.faDffPostAirLick.ss{t, 1}, rez.faDffPostAirLick.ss{t, 2}] = alignToEvent(dffSs, tbytDat(t).postAirLicks(1), frT, [-1 1]); % postAir lick 1st in the bout
            % v1
            [rez.faDffPostAirLick.v1{t, 1}, rez.faDffPostAirLick.v1{t, 2}] = alignToEvent(dffV1, tbytDat(t).postAirLicks(1), frT, [-1 1]); % postAir lick 1st in the bout
            % rs
            [rez.faDffPostAirLick.rs{t, 1}, rez.faDffPostAirLick.rs{t, 2}]  = alignToEvent(dffRs, tbytDat(t).postAirLicks(1), frT, [-1 1]); % postAir lick 1st in the bout
        end

    % Correct Rejection
    elseif tbytDat(t).punishTrI && isempty(tbytDat(t).airpuff)   
        % m1
        [rez.crCueDff.m1{t, 1}, rez.crCueDff.m1{t, 2}] = alignToEvent(dffM1, 0, frT, [-1 1]); % align to cue onset
        % m2
        [rez.crCueDff.m2{t, 1}, rez.crCueDff.m2{t, 2}]  = alignToEvent(dffM2, 0, frT, [-1 1]); % align to cue onset
        % ss
        [rez.crCueDff.ss{t, 1}, rez.crCueDff.ss{t, 2}] = alignToEvent(dffSs, 0, frT, [-1 1]); % align to cue onset
        % v1
        [rez.crCueDff.v1{t, 1}, rez.crCueDff.v1{t, 2}] = alignToEvent(dffV1, 0, frT, [-1 1]); % align to cue onset
        % rs
        [rez.crCueDff.rs{t, 1}, rez.crCueDff.rs{t, 2}] = alignToEvent(dffRs, 0, frT, [-1 1]); % align to cue onset
    end
    fprintf("processed trial #%d\n", t)
end

% save rez
save(fullfile(filePath, 'Matfiles', strcat(header, '_dff_evtAligned_regionMask.mat')), 'rez'); 

%% Stim onset aligned activity (M1)
% align dffs take the mean and sem
% primary motor (stim onset)
[stimOnDff_m1_itp, timepts] = temporalAlignInterp1(rez.stimOnDffC.m1(:, 1), rez.stimOnDffC.m1(:, 2));
[rez.meanStimOnDff.m1, ~, rez.semStimOnDff.m1] = meanstdsem(cell2mat(stimOnDff_m1_itp));

% primary motor (stim onset) go no-go
[rez.meanStimOnDffGng.m1, rez.semStimOnDffGng.m1] = trialGroupMeanSem(stimOnDff_m1_itp, {[tbytDat.rewardTrI], [tbytDat.punishTrI]});
plotMeanSem(rez.meanStimOnDffGng.m1, rez.semStimOnDffGng.m1, timepts, {'Go', 'NoGo'});
title("Go Nogo M1 (cue onset at time=0)");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:4, 'TickDir', 'out', 'YTick', -1:0.5:1); xlim([-2 2]); ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_cueOn_Gng_m1.pdf'), '-dpdf', '-vector', '-bestfit')

% primary motor (stim onset) go no-go (logistic regression)
rez.logitGngCueOnM1 = logireg(cell2mat(stimOnDff_m1_itp), [tbytDat.rewardTrI]); % logistic regression

%% iti licks (M1)
itiLickDffM1 = flatCell(rez.itiLickDff.m1); 
itiLickDffM1Ts = cellfun(@(a) linspace(-1, 1, length(a)), itiLickDffM1(:, 2), 'UniformOutput', false); 
[itiLickDffM1Itp, itiLickDffM1ItpTs] = temporalAlignInterp1(itiLickDffM1 (:, 1), itiLickDffM1Ts);
[rez.meanItiLickDff.m1, ~, rez.semItiLickDff.m1] = meanstdsem(cell2mat(itiLickDffM1Itp));



cellfun(@length, itiLickDffM1, 'UniformOutput', false)


rez.itiLickDff.m1



itiLickTsM1 = cellfun(@(a) a(:, 2), rez.itiLickDff.m1(itiLickI), 'UniformOutput', false); 







%% Stim onset aligned activity (M2)
% supplementary motor (stim onset)
stimOnDff_m2_itp = temporalAlignInterp1(stimOnDffC.m2(:, 1), stimOnDffC.m2(:, 2));
[rez.mStimOnDffM2, ~, rez.semStimOnDffM2] = meanstdsem(cell2mat(stimOnDff_m2_itp));

% supplementary motor (stim onset) go no-go
[rez.mStimOnDffGngM2, rez.sStimOnDffGngM2] = trialGroupMeanSem(stimOnDff_m2_itp, {[tbytDat.rewardTrI], [tbytDat.punishTrI]});
plotMeanSem(rez.mStimOnDffGngM2, rez.sStimOnDffGngM2, timepts, {'Go', 'NoGo'});
title("Go Nogo M2 (cue onset at time=0)");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:4, 'TickDir', 'out', 'YTick', -1:0.5:1); xlim([-2 2]); ylim([-1 1])
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
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:4, 'TickDir', 'out', 'YTick', -1:0.5:1); xlim([-2 2]); ylim([-1 1])
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
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:4, 'TickDir', 'out', 'YTick', -1:0.5:1); xlim([-2 2]); ylim([-1 1])
print(fullfile(filePath, 'Figure', 'dff_cueOn_Gng_v.pdf'), '-dpdf', '-vector', '-bestfit')

% visual (stim onset) go no-go (logistic regression)
rez.logitGngCueOnV = logireg(cell2mat(stimOnDff_v_itp), [tbytDat.rewardTrI]); % logistic regression

%% retrosplenial cortex (RS)
stimOnDff_rs_itp = temporalAlignInterp1(stimOnDffC.rs(:, 1), stimOnDffC.rs(:, 2));
[rez.mStimOnDffRs, ~, rez.semStimOnDffRs] = meanstdsem(cell2mat(stimOnDff_rs_itp));

% retrosplenial (stim onset) go no-go
[rez.mStimOnDffGngRs, rez.sStimOnDffGngRs] = trialGroupMeanSem(stimOnDff_rs_itp, {[tbytDat.rewardTrI], [tbytDat.punishTrI]});
plotMeanSem(rez.mStimOnDffGngRs, rez.sStimOnDffGngRs, timepts, {'Go', 'NoGo'});
title("Go Nogo Rs (cue onset at time=0)");
xlabel('Time (s)'); ylabel('DFF'); set(gca, 'XTick', -2:1:4, 'TickDir', 'out', 'YTick', -1:0.5:1); xlim([-2 2]); ylim([-1 1])
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
plotMeanSem(smooth2a([rez.logitGngCueOnM1; rez.logitGngCueOnM2; rez.logitGngCueOnS1; rez.logitGngCueOnRs; rez.logitGngCueOnV], 0, 2), ...
    zeros(5, length(timepts)), timepts, {'M1', 'M2', 'S1', 'RS', 'V'});
title("Cross-validated classification accuracy")
xlabel('Time (s)')
ylabel('Classification Accuracy')
set(gca, 'XTick', -2:1:4, 'YTick', -1:0.1:1, 'TickDir', 'out')
xlim([-2 2])
ylim([0.45 0.85])
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




end


