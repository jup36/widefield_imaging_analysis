function [tbytEvt, params, options] = dffCombinedProcessingDual(filePath, dffC, tbytEvt, options, params)
%This function takes the preprocessed raw image 'stack' (ome_stack) as a cell and process them
% using the parameters specified in the options and other user defined inputs.
% INPUT:
%   dffC: a cell array containing preprocessed/dff image stacks, each stack
%       corresponds to a trial or a block of trials depending on the task and image acquisition scheme.
%   cmosExp: timestamps for cmos exposure pulses (1st col: timestamps, 2nd col: chunk (trial or block of trial) id).
%   options: prepro_log.mat
%   params: parameters for stack processing

% whereabouts
match_header = regexp(filePath, 'm\d{1,4}_\d{6}', 'match');
header = match_header{1};

if params.saveTbytDff
    % file directory for trials
    filePathTrials = fullfile(filePath, strcat(header, '_trials'));
    if exist(filePathTrials, 'dir') == 0
        mkdir(filePathTrials);
    end
end

%% align frames per each trial and get trial-by-trial dff
for tt = 1:length(tbytEvt) % increment trials
    switch params.channelOfInterest
        case 'green'
            if ~isempty(tbytEvt(tt).blueLED)
                dff = dffC{1, tbytEvt(tt).blueLEDTrainI};
                tempBlueLEDIdx = tbytEvt(tt).blueLEDPulsesOfTrain{1}:...
                    min(tbytEvt(tt).blueLEDPulsesOfTrain{end}, size(dff,3));
                tbytEvt(tt).dff = dff(:,:,tempBlueLEDIdx);
                tbytEvt(tt).dffsm = applyImgaussfilt(tbytEvt(tt).dff);
                tbytEvt(tt).frameT = tbytEvt(tt).blueLED(1:length(tempBlueLEDIdx)); % frame time
                tbytEvt(tt).frameTrel = tbytEvt(tt).frameT-tbytEvt(tt).evtOn; % frame timestamps relative to the event onset
                if params.saveTbytDff
                    trSubDir = fullfile(filePathTrials, sprintf('dffG_block_%d_trial_%d', tbytEvt(tt).blueLEDTrainI, tt));
                    if exist(trSubDir, 'dir') ~= 7
                        mkdir(trSubDir);
                    end
                    tbytDff = tbytEvt(tt).dff;
                    tbytDffsm = tbytEvt(tt).dffsm;
                    save(fullfile(trSubDir, 'tbytDff.mat'), 'tbytDff', 'tbytDffsm')
                    clearvars tbytDff tbytDffsm
                end
            end

        case 'red'
            if ~isempty(tbytEvt(tt).limeLED)
                dff = dffC{1, tbytEvt(tt).limeLEDTrainI};
                tempLimeLEDIdx = tbytEvt(tt).limeLEDPulsesOfTrain{1}:...
                    min(tbytEvt(tt).limeLEDPulsesOfTrain{end}, size(dff,3));
                tbytEvt(tt).dff = dff(:,:,tempLimeLEDIdx);
                tbytEvt(tt).dffsm = applyImgaussfilt(tbytEvt(tt).dff);
                tbytEvt(tt).frameT = tbytEvt(tt).limeLED(1:length(tempLimeLEDIdx)); % frame time
                tbytEvt(tt).frameTrel = tbytEvt(tt).frameT-tbytEvt(tt).evtOn; % frame timestamps relative to the event onset
                if params.saveTbytDff
                    trSubDir = fullfile(filePathTrials, sprintf('dffR_block_%d_trial_%d', tbytEvt(tt).limeLEDTrainI, tt));
                    if exist(trSubDir, 'dir') ~= 7
                        mkdir(trSubDir);
                    end
                    tbytDff = tbytEvt(tt).dff;
                    tbytDffsm = tbytEvt(tt).dffsm;
                    save(fullfile(trSubDir, 'tbytDff.mat'), 'tbytDff', 'tbytDffsm')
                    clearvars tbytDff tbytDffsm
                end
            end
    end
    fprintf('completed trial %d/%d\n', tt, length(tbytEvt))
end
end