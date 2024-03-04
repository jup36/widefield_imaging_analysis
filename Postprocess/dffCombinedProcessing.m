function [tbytEvt, params, options] = dffCombinedProcessing(filePath, dffC, tbytEvt, options, cmosExp, params)
%This function takes the preprocessed raw image 'stack' (ome_stack) as a cell and process them
% using the parameters specified in the options and other user defined inputs.
% INPUT:
%   dffC: a cell array containing preprocessed/dff image stacks, each stack
%       corresponds to a trial or a block of trials depending on the task and image acquisition scheme.
%   cmosExp: timestamps for cmos exposure pulses (1st col: timestamps, 2nd col: chunk (trial or block of trial) id).
%   options: prepro_log.mat
%   params: parameters for stack processing


% whereabouts
[~, header] = fileparts(filePath);
if params.saveTbytDff
    % file directory for trials
    filePathTrials = fullfile(filePath, strcat(header, '_trials'));
    if exist(filePathTrials, 'dir') == 0
        mkdir(filePathTrials);
    end
end

%% align frames per each trial and get trial-by-trial dff
for tt = 1:length(tbytEvt) % increment trials
    if ~isempty(tbytEvt(tt).cmosExp)
        cmosI = tbytEvt(tt).cmosExpTrainI; % cmos pulse train index
        dff = dffC{1, tbytEvt(tt).cmosExpTrainI};
        tbytEvt(tt).frameT = tbytEvt(tt).cmosExp(1:2:end); % frame time (needs to alternate due to interleaved violet frames for hemodynamic correction)
        temp1stFrameI = floor(tbytEvt(tt).cmosExpPulsesOfTrain{1}./2); % 1st frame to take in dff (indices must be divided by 2 since dff already taken excluding violet frames)
        tbytEvt(tt).frameI = temp1stFrameI:temp1stFrameI+length(tbytEvt(tt).frameT)-1;
        tbytEvt(tt).dff = dff(:,:,tbytEvt(tt).frameI); % aligned dff
        tbytEvt(tt).dffsm = applyImgaussfilt(tbytEvt(tt).dff);
        tbytEvt(tt).frameTrel = tbytEvt(tt).frameT-tbytEvt(tt).evtOn; % frame timestamps relative to the event onset

        % save trial-by-trial dff
        if params.saveTbytDff
            trSubDir = fullfile(filePathTrials, sprintf('block_%d_trial_%d', tbytEvt(tt).cmosExpTrainI, tt));
            if exist(trSubDir, 'dir') ~= 7
                mkdir(trSubDir);
            end
            tbytDff = tbytEvt(tt).dff;
            tbytDffsm = tbytEvt(tt).dffsm;
            save(fullfile(trSubDir, 'tbytDff.mat'), 'tbytDff', 'tbytDffsm')
            clearvars tbytDff tbytDffsm
        end
    end
    fprintf('completed trial %d/%d\n', tt, length(tbytEvt))
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [stack_b, stack_v, bFrameI, vFrameI] = getBVframes(stack)

        numbFrames = size(stack, 3);

        stack_b = stack(:, :, 1:2:end);
        bFrameI = 1:2:numbFrames;

        stack_v = stack(:, :, 2:2:end);
        vFrameI = 2:2:numbFrames;

        %Trim until both the same length (may be one frame discrepancy)
        min_length = min(cell2mat(cellfun(@(x) size(x,3), {stack_b,stack_v},'UniformOutput',0)));
        stack_b = stack_b(:,:,1:min_length);
        bFrameI = bFrameI(1:min_length);

        stack_v = stack_v(:,:,1:min_length);
        vFrameI = vFrameI(1:min_length);

        % to ensure correct assignment of blue and violet frames check their mean intensities
        if nanmean(stack_b(:)) < nanmean(stack_v(:))
            stack_b_copy = stack_b;
            bFrameI_copy = bFrameI;
            stack_b = stack_v;
            bFrameI = vFrameI;
            stack_v = stack_b_copy;
            vFrameI = bFrameI_copy;
            clearvars stack_b_copy
        end

        %remove the masked pixels by setting to NaN
        masked_pxls = (stack_b==0);
        stack_b(masked_pxls)=NaN;
        stack_v(masked_pxls)=NaN;

    end



end