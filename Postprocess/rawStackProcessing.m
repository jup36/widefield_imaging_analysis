function [tbytEvt, params, options] = rawStackProcessing(filePath, stackC, tbytEvt, options, cmosExp, params)
%This function takes the preprocessed raw image 'stack' (ome_stack) as a cell and process them
% using the parameters specified in the options and other user defined inputs.
% INPUT:
%   stackC: a cell array containing preprocessed image stacks, each stack
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

% sanity check on input parameters regarding dff calculation
if ismember(params.dffMethod, {'movingavg', 'mean', 'median', 'mode'})
    % override (for any potential discrepancy between user-defined vs. prepro_log options) the params for dff in options with the user-defined
    % parameters
    options.method = params.dffMethod;
    options.method_window = params.movingWinF;
    fprintf(strcat('dff is calculated across trials using', ' movingavg', '!\n'))
elseif ismember(params.dffMethod, {'tbytBase'})
    fprintf(strcat('dff is calculated trial by trial!\n'))
else
    error('dffMethod must be set as tbytBase or movingavg!')
end

%% get signal (blue) and correction (violet) frames separately, as an option perform dff if using 'movingavg'
for i = 1:length(stackC) % increment stacks
    [stackS(i).stack_b, stackS(i).stack_v, stackS(i).bFrameI, stackS(i).vFrameI] = getBVframes(stackC{i});
    %there are large transients in both signals for the first 10 seconds of the recording, likely due to LED warm up (should add 'burn in' frames in the future);
    %this can effects downstream results so replace the first 10 seconds of both with the average of the first minute; You should remove these at the end of analysis (after timelocking to behavior/ephys)
    stackS(i).stack_b(:,:,1:5*options.fps) = repmat(nanmean(stackS(i).stack_b(:,:,1:30*options.fps),3),1,1,5*options.fps);
    stackS(i).stack_v(:,:,1:5*options.fps) = repmat(nanmean(stackS(i).stack_v(:,:,1:30*options.fps),3),1,1,5*options.fps);

    % get fitted (scaled) violet stack
    [stackS(i).stack_b, stackS(i).stack_v_fitted] = getPixelwiseScaledV(stackS(i).stack_b, stackS(i).stack_v, options);

    % in case calculating dff using movingavg rolling over each stack
    if ismember(params.dffMethod, {'movingavg', 'mean', 'median', 'mode'})
        stackB_block_dff = makeDFF(stackS(i).stack_b, options);
        if params.bvCorrectLogic
            stackV_block_dff = makeDFF(stackS(i).stack_v_fitted, options);
            stackS(i).dff = stackB_block_dff-stackV_block_dff;
        else
            stackS(i).dff = stackB_block_dff;
        end
    end
    fprintf('completed stackS %d/%d\n', i, length(stackC))
end

%% align frames per each trial and get trial-by-trial dff
for tt = 1:length(tbytEvt) % increment trials
    if ~isempty(tbytEvt(tt).cmosExp)
        cmosI = tbytEvt(tt).cmosExpTrainI; % cmos pulse train index
        cmosExpBlock = cmosExp(cmosExp(:, 2)==cmosI, 1); % cmos timestamps of the cmosI
        pulseEdges = tbytEvt(tt).cmosExpPulsesOfTrain{1}:tbytEvt(tt).cmosExpPulsesOfTrain{end}; % edges defining cmos pulses relavant to this trial

        bFrameI_trial = ismember(stackS(cmosI).bFrameI, pulseEdges); % blue frames of the cmosI that correspond to this trial
        vFrameI_trial = ismember(stackS(cmosI).vFrameI, pulseEdges); % violet frames of the cmosI that correspond to this trial

        % match frame numbers (crop)
        if sum(bFrameI_trial)>sum(vFrameI_trial)
            bFrameI_trial(find(bFrameI_trial, 1, 'last'))=false;
        elseif sum(bFrameI_trial)<sum(vFrameI_trial)
            vFrameI_trial(find(vFrameI_trial, 1, 'last'))=false;
        end

        % frame ID
        tbytEvt(tt).bFramesInTrain = stackS(cmosI).bFrameI(bFrameI_trial); % frame ID in the cmos pulse train
        tbytEvt(tt).vFramesInTrain = stackS(cmosI).bFrameI(bFrameI_trial); % frame ID in the cmos pulse train

        % frame timestamps
        bFrameTs_trial = cmosExpBlock(stackS(cmosI).bFrameI(bFrameI_trial)); % blue frame timestamps
        vFrameTs_trial = cmosExpBlock(stackS(cmosI).vFrameI(vFrameI_trial)); % violet frame timestamps

        tbytEvt(tt).frameT = bFrameTs_trial; % frame timestamps (blue frames)
        tbytEvt(tt).frameTrel = tbytEvt(tt).frameT-tbytEvt(tt).evtOn; % frame timestamps relative to the event onset

        % get dff
        % sanity check on input parameters
        if strcmpi(params.dffMethod, 'tbytBase') % dff using the trial-by-trial baseline
            stackB_trial = stackS(cmosI).stack_b(:, :, bFrameI_trial);
            stackV_trial = stackS(cmosI).stack_v_fitted(:, :, vFrameI_trial); % ensure to use stack_v_fitted for proper hemodynamic correction

            % make trial-by-trial dff
            if strcmpi(params.tbytBaseWinF, 'all')
                stack_base_logic = tbytEvt(tt).frameTrel<0; % use all negatively timestamped frames, i.e., all frames before the event onset
            elseif isnumeric(params.tbytBaseWinF)
                stack_base_logic = tbytEvt(tt).frameTrel < min(tbytEvt(tt).frameTrel(1)+params.tbytBaseWinF, 0); % use negatively timestamped left-most frames within the specified window width
            end

            stackB_trial_dff = makeDFF_trial(stackB_trial, stack_base_logic, 'dffMethod', 'mean', 'type', 'dff', 'detrend', false);

            if params.bvCorrectLogic
                stackV_trial_dff = makeDFF_trial(stackV_trial, stack_base_logic, 'dffMethod', 'mean', 'type', 'dff', 'detrend', false);
                tbytEvt(tt).dff = stackB_trial_dff-stackV_trial_dff;
            else
                tbytEvt(tt).dff = stackB_trial_dff;
            end
        elseif ismember(params.dffMethod, {'movingavg', 'mean', 'median', 'mode'}) % dff using movingavg or others that use across-trial (e.g., blockwise) F
            tbytEvt(tt).dff = stackS(cmosI).dff(:, :, bFrameI_trial);
        else
            error('dffMethod must be set as tbytBase or movingavg!')
        end

        % pixel-wise smoothing
        tbytEvt(tt).dffsm = applyImgaussfilt(tbytEvt(tt).dff);

        if isempty(tbytEvt(tt).dff)
            tbytEvt(tt).dff = NaN; 
        end

        if isempty(tbytEvt(tt).dffsm)
            tbytEvt(tt).dffsm = NaN; 
        end

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