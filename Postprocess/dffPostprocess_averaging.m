function dffPostprocess_averaging(filePath)
%filePath = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA003/DA003_083023/DA003_083023_img';
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
parentDir = fileparts(filePath);
[~, header] = fileparts(parentDir);
fileBeh = GrabFiles_sort_trials('tbytDat', 1, {parentDir});
load(fullfile(fileBeh{1}), 'tbytDat')



%parentDir = fileparts(filePath);
[~, header] = fileparts(filePath);
fileBeh = GrabFiles_sort_trials('tbytDat_dff', 0, {fullfile(filePath, 'Matfiles')});
if isempty(fileBeh)
    fileBeh = GrabFiles_sort_trials('tbytDat_dff', 1, {filePath});
end
load(fullfile(fileBeh{1}), 'tbytDat')



% load the preprocessed dffs
[file_list_dff, folder_list_dff] = GrabFiles_sort_trials('dff_combined.mat', 1, {filePath}); % use GrabFiles_sort_trials to sort both files and folders
dffC = cell(1, length(file_list_dff));
for ff = 1:length(file_list_dff)
    load(file_list_dff{ff}, 'dff')
    dffC{1, ff} = dff;
    clearvars dff
    fprintf("finished loading dff file #%d\n", ff)
end

% load Allen dorsalMap
load('allenDorsalMap', 'dorsalMaps');

% load Allen transformation object
load(fullfile(folder_list_dff{1}, 'transParamsAllen.mat'), 'transParams')

% file directory for trials
filePathTrials = fullfile(parentDir, strcat(header, '_trials'));
if exist(filePathTrials, 'dir') == 0
    mkdir(filePathTrials);
end

%refCmosFrameIdx = 1:2700; % there must be 2700 cmos exposure pulses / frames recorded

% take corresponding frames with for each trial with 2D gaussian filtering
for tt = 1:length(tbytDat)
    if ~isempty(tbytDat(tt).cmosExp)
        dff = dffC{1, tbytDat(tt).cmosExpTrainI};
        tbytDat(tt).frameT = tbytDat(tt).cmosExp(1:2:end); % frame time (needs to alternate due to interleaved violet frames for hemodynamic correction)
        temp1stFrameI = floor(tbytDat(tt).cmosExpPulsesOfTrain{1}./2); % 1st frame to take in dff (indices must be divided by 2 since dff already taken excluding violet frames)
        tbytDat(tt).frameI = temp1stFrameI:temp1stFrameI+length(tbytDat(tt).frameT)-1;
        tbytDat(tt).dff = dff(:,:,tbytDat(tt).frameI); % aligned dff
        tbytDat(tt).dffsm = applyImgaussfilt(tbytDat(tt).dff);

        % map temporal events to cmosExp pulses
        tbytDat(tt).frameLickI = check_timestamp_overlap(tbytDat(tt).frameT, tbytDat(tt).licks);
        tbytDat(tt).frameWaterI = check_timestamp_overlap(tbytDat(tt).frameT, tbytDat(tt).water);

        % timestamp each frame relative to the stim onset
        if tbytDat(tt).evtType < 3  % 1 or 2: visual (common, uncommon)
            tbytDat(tt).frameStimOnI = check_timestamp_overlap(tbytDat(tt).frameT, tbytDat(tt).evtOn);
            tbytDat(tt).frameStimOffI = check_timestamp_overlap(tbytDat(tt).frameT, tbytDat(tt).evtOff);
            tbytDat(tt).frameStimI = zeros(length(tbytDat(tt).frameT), 1);
            tbytDat(tt).frameStimI(find(tbytDat(tt).frameStimOnI, 1):find(tbytDat(tt).frameStimOffI, 1), 1) = 1;
        end
        tbytDat(tt).frameTrel = tbytDat(tt).frameT-tbytDat(tt).evtOn;
        tbytDat(tt).dir = folder_list_dff{tbytDat(tt).cmosExpTrainI};
    end
    fprintf('processed dff of trial#%d\n', tt)
end

% Trial indices
rwdI = cell2mat(cellfun(@(a) a==3, {tbytDat.evtType}, 'UniformOutput', false)); % reward trials
dffI = cell2mat(cellfun(@(a) ~isempty(a), {tbytDat.dffsm}, 'UniformOutput', false)); % trials with dff
rwdDffI = rwdI & dffI; 

% just temporally align to water delivery 
%tbytDat(21).dffsm
%tbytDat(21).dffsm(isnan(tbytDat(21).dffsm))


% transform frames 
%wfa = alignStackToAllenKabsch(tbytDat(21).dffsm, dorsalMaps.dorsalMap, transParams.tformObj); 

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


rwdDffC_motor = cell(length(tbytDat), 4); 
rwdDffC_smotor = cell(length(tbytDat), 4); 
rwdDffC_ss = cell(length(tbytDat), 4); 
rwdDffC_v = cell(length(tbytDat), 4); 

for t = 1:length(tbytDat)
    if rwdDffI(t) && tbytDat(t).cmosExpTrainI~=2 % for some reason, the block 2 was problematic, so just exclude now for DA002_083023
        wfa = alignStackToAllenKabsch(tbytDat(t).dffsm, dorsalMaps.dorsalMap, transParams.tformObj); 
       
        timeStamps = tbytDat(t).frameTrel; % store timestamps
        baseI = timeStamps<0;
        rwdI = timeStamps>=0 & timeStamps<1; 
        
        % Motor (M1) 
        rwdDffC_motor{t, 1} = apply2DMaskTo3DStack(wfa, motorMask); % get the mean masked data across time
        rwdDffC_motor{t, 2} = timeStamps; 
        rwdDffC_motor{t, 3} = nanmean(wfa(:, :, baseI), 3); % base dff
        rwdDffC_motor{t, 4} = nanmean(wfa(:, :, rwdI), 3); % rwd dff

        % Supplementary Motor (M2) 
        rwdDffC_smotor{t, 1} = apply2DMaskTo3DStack(wfa, smotorMask); % get the mean masked data across time
        rwdDffC_smotor{t, 2} = timeStamps; 

        % Somatosensory (S1) 
        rwdDffC_ss{t, 1} = apply2DMaskTo3DStack(wfa, ssMask); % get the mean masked data across time
        rwdDffC_ss{t, 2} = timeStamps; 

        % Visual (V1) 
        rwdDffC_v{t, 1} = apply2DMaskTo3DStack(wfa, vMask); % get the mean masked data across time
        rwdDffC_v{t, 2} = timeStamps; 
    end
    fprintf("processed trial #%d\n", t)
end 
save(fullfile(parentDir, 'Matfiles', 'rwdDffC_motor.mat'), 'rwdDffC_motor', '-v7.3')
save(fullfile(parentDir, 'Matfiles', 'rwdDffC_smotor.mat'), 'rwdDffC_smotor', '-v7.3')
save(fullfile(parentDir, 'Matfiles', 'rwdDffC_ss.mat'), 'rwdDffC_ss', '-v7.3')
save(fullfile(parentDir, 'Matfiles', 'rwdDffC_v.mat'), 'rwdDffC_v', '-v7.3')

% align dffs take the mean and sem 
% primary motor
[rwdDffC_motor_itp, timepts] = temporalAlignInterp1(rwdDffC_motor(:, 1), rwdDffC_motor(:, 2)); 
[rez.meanRwdDffMotor, ~, rez.semRwdDffMotor] = meanstdsem(cell2mat(rwdDffC_motor_itp)); 
plotMeanSEM(rez.meanRwdDffMotor, rez.semRwdDffMotor, timepts);
title('dff aligned to reward (M1)')
xlabel('Time (s)')
ylabel('DFF')
set(gca, 'XTick', -1:1:4)
print(fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA003/DA003_083023/Figure', 'dff_reward_M1.pdf'), '-dpdf', '-vector', '-bestfit')

% supplementary motor
[rwdDffC_smotor_itp, timepts] = temporalAlignInterp1(rwdDffC_smotor(:, 1), rwdDffC_smotor(:, 2)); 
[rez.meanRwdDffSMotor, ~, rez.semRwdDffSMotor] = meanstdsem(cell2mat(rwdDffC_smotor_itp)); 
plotMeanSEM(rez.meanRwdDffSMotor, rez.semRwdDffSMotor, timepts);
title('dff aligned to reward (M2)')
xlabel('Time (s)')
ylabel('DFF')
set(gca, 'XTick', -1:1:4)
print(fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA003/DA003_083023/Figure', 'dff_reward_M2.pdf'), '-dpdf', '-vector', '-bestfit')

% somatosensory
[rwdDffC_ss_itp, timepts] = temporalAlignInterp1(rwdDffC_ss(:, 1), rwdDffC_ss(:, 2)); 
[rez.meanRwdDffSensory, ~, rez.semRwdDffSensory] = meanstdsem(cell2mat(rwdDffC_ss_itp)); 
plotMeanSEM(rez.meanRwdDffSensory, rez.semRwdDffSensory, timepts);
title('dff aligned to reward (Somatosensory)')
xlabel('Time (s)')
ylabel('DFF')
set(gca, 'XTick', -1:1:4)
print(fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA003/DA003_083023/Figure', 'dff_reward_SS.pdf'), '-dpdf', '-vector', '-bestfit')

% visual
[rwdDffC_v_itp, timepts] = temporalAlignInterp1(rwdDffC_v(:, 1), rwdDffC_v(:, 2)); 
[rez.meanRwdDffVisual, ~, rez.semRwdDffVisual] = meanstdsem(cell2mat(rwdDffC_v_itp)); 
plotMeanSEM(rez.meanRwdDffVisual, rez.semRwdDffVisual, timepts);
title('dff aligned to reward (Visual)')
xlabel('Time (s)')
ylabel('DFF')
set(gca, 'XTick', -1:1:4)
print(fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA003/DA003_083023/Figure', 'dff_reward_V.pdf'), '-dpdf', '-vector', '-bestfit')

plotMeanSEM([rez.meanRwdDffMotor; rez.meanRwdDffSMotor; rez.meanRwdDffSensory; rez.meanRwdDffVisual], ... 
            [rez.semRwdDffMotor; rez.semRwdDffSMotor; rez.semRwdDffSensory; rez.semRwdDffVisual], timepts);
print(fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA003/DA003_083023/Figure', 'dff_reward_Collect.pdf'), '-dpdf', '-vector', '-bestfit')



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

