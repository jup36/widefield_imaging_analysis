%% Refit Chunked Motif Data and Combine Results
% ------------------------------------------------------------------------------
% This script processes dual-imaging data across multiple sessions and animals,
% specifically for calcium motifs. It performs the following key steps:
%
% 1. Iterates over animal folders containing dual-imaging data.
% 2. Identifies session folders and verifies that train/test motif chunk files 
%    are organized in matched pairs (odd = train, even = test).
% 3. Loads refitted CNMF motifs (W, H) and associated metadata (bad pixels, stats).
% 4. Extracts and conditions raw dF/F traces (removes NaNs from bad pixels).
% 5. Shifts each motif’s temporal weight (H) according to its spatial peak lag in W,
%    so that H reflects the peak expression time.
% 6. Timestamps H using LED pulse train metadata (via tbytDat).
% 7. Saves the compiled chunk-wise data into a single .mat file per session:
%    `*_refitChunks_red_dff_combined.mat`
%
% Prerequisite:
%   This script assumes that all motif chunks have been independently refit 
%   using the script `batchProcess_refitBasisMotifsScotty_redCalcium.m`.
%
% Output per session:
%   - hC:      Cell array of motif temporal weights [original, shifted, timestamps]
%   - wC:      Cell array of motif spatial components
%   - nanpxsC: Cell array of bad-pixel masks
%   - dffC:    Conditioned dF/F data with NaNs removed
%   - stats_refitC: CNMF statistics for each chunk
%   - tbytDat: Trial-aligned behavioral data
%   - trI:     Trial-type index matrix
%
% Author: Junchol Park, 7/31/25
% ------------------------------------------------------------------------------
%% whereabouts
filePathBase = compatiblepath('/Volumes/buschman/Rodent Data/dualImaging_parkj');
filePreprocessed = compatiblepath('/Volumes/buschman/Rodent Data/Wide Field Microscopy/ExampleData/Preprocessed');
fileKeyword = '_red_dff_combined.mat';

% Search for all animal folders containing dual-imaging data
mListC = GrabFiles_sort_trials('m*jRGECO', 0, {filePathBase});

%% Main (Batch)
for j = 1:numel(mListC) % Iterate over each animal

    % Extract mouse ID (e.g., m1045) from folder name
    mId = cell2mat(regexp(mListC{j}, 'm\d{4}', 'match'));

    % Find session folders within the animal folder
    filePathSessions = find_keyword_containing_folder(mListC{j}, mId, 'recursive', false);

    for jj = 1:numel(filePathSessions)
        filePath = cell2mat(find_keyword_containing_folder(filePathSessions{jj}, 'task', 'recursive', false));
        filePath_img = cell2mat(find_keyword_containing_folder(filePath, '_img', 'recursive', false));
        [~, header] = fileparts(filePath_img);
        filePathPreprop = cell2mat(find_keyword_containing_folder(filePreprocessed, header, 'recursive', false));

        if ~isempty(filePathPreprop) % proceed only if the motifs exist
            filePath_mat = fullfile(filePath, 'Matfiles');
            filePath_fnC = cell2mat(GrabFiles_sort_trials(['_img_list' fileKeyword], 0, {filePath_img}));
            traintestLogic = verifyTrainTestSets(filePath_fnC); % check the train (odd numbers) and test (even numbers) set organization
            assert(traintestLogic)
            fileList_tbytDat = cell2mat(GrabFiles_sort_trials('_tbytDat_parseGng', 0, {filePath_mat}));
            filePathPre = cell2mat(GrabFiles_sort_trials(fileKeyword, 0, {filePathPreprop}));
            preS = load(filePathPre, 'data_train', 'data_test');
            filePath_figure = fullfile(filePath, 'Figure');

            %% load tbytDat
            load(fileList_tbytDat, 'tbytDat')
            trI = trialTypeInfoAuditoryGng(filePath);

            %% load refit data
            fileList = GrabFiles_sort_trials('_red_dff_combined_refitchunk_', 0, {filePath_mat});
            chunkN = ceil(numel(fileList)/2);

            hC = cell(1, numel(fileList));
            wC = cell(1, numel(fileList));
            nanpxsC = cell(1, numel(fileList));
            stats_refitC = cell(1, numel(fileList));
            dffC = cell(1, numel(fileList));

            %% Main loop
            for i = 1:chunkN
                temp_train_file = cell2mat(GrabFiles_sort_trials(sprintf('_red_dff_combined_refitchunk_%dtrain.mat', i), 0, {filePath_mat}));
                temp_test_file = cell2mat(GrabFiles_sort_trials(sprintf('_red_dff_combined_refitchunk_%dtest.mat', i), 0, {filePath_mat}));

                if ~isempty(temp_train_file)
                    temp_train = load(temp_train_file, 'H', 'w', 'bad_pxl', 'stats_refit');
                    hC{1, 2*i-1} = temp_train.H;
                    wC{1, 2*i-1} = temp_train.w;
                    nanpxsC{1, 2*i-1} = temp_train.bad_pxl;
                    stats_refitC{1, 2*i-1} = temp_train.stats_refit;
                    dffC{1, 2*i-1} = conditionDffMat(preS.data_train(:,:,i)', find(nanpxsC{i}));
                    clearvars temp_train
                end

                if ~isempty(temp_test_file)
                    temp_test  = load(temp_test_file, 'H', 'w', 'bad_pxl', 'stats_refit');
                    hC{1, 2*i} = temp_test.H;
                    wC{1, 2*i} = temp_test.w;
                    nanpxsC{1, 2*i} = temp_test.bad_pxl;
                    stats_refitC{1, 2*i} = temp_test.stats_refit;
                    dffC{1, 2*i} = conditionDffMat(preS.data_test(:,:,i)', find(nanpxsC{i}));
                    clearvars temp_test
                end
            end

            %% temporal shifting of H by the peak expression lags in W
            hC_shift = cellfun(@(a, b) motifExpressionPeak(a, b), wC(1,:), hC(1,:), 'UniformOutput', false);
            hC = [hC; hC_shift];

            %% timestamping pulses using LED pulses
            for tr = 1:length(hC)
                hC{3, tr} = NaN(1, size(hC{1, tr}, 2));

                % iterate through pulse trains
                trialI = cellfun(@(a) a==tr, {tbytDat.limeLEDTrainI});
                ledPulsesOfTrain = {tbytDat(trialI).limeLEDPulsesOfTrain};
                ledPulsesTimeOfTrain = {tbytDat(trialI).limeLED};

                for tt = 1:sum(trialI)
                    pulseEdge = cell2mat(ledPulsesOfTrain{tt});
                    pulseTime = ledPulsesTimeOfTrain{tt};
                    hC{3, tr}(pulseEdge(1):pulseEdge(end)) = pulseTime;
                end
                hC{3, tr} = hC{3, tr}(1, 1:size(hC{1, tr}, 2));
                hC{3, tr} = fillNaNtimestamps(hC{3, tr});
            end

            %% save
            saveName = strcat(header, "_refitChunks", fileKeyword); % e.g., "m1044_121924_task_day4-5_img_refitChunks_red_dff_combined.mat"
            save(fullfile(filePath_mat, saveName), "hC", "wC", "nanpxsC", "dffC", "stats_refitC", "tbytDat", "trI", "-v7.3") 

        end
        fprintf("Completed saving the refitted motif chunk data for session#%d of %s\n", jj, mId)
    end
end


%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%
function outLogic = verifyTrainTestSets(filePath_fnC)
load(filePath_fnC, 'fnC')
fileNames = cellfun(@getFileNameOnly, fnC, 'UniformOutput', false);
pattern = '_([0-9]{1,2})_';
numberStr = cellfun(@(f) cell2mat(regexp(f, pattern, 'tokens', 'once')), fileNames, 'UniformOutput', false);
numberVal = cellfun(@(x) str2double(x), numberStr);

% To assert first column is odd second column equals to 'first column + 1'
outLogic = all(mod(numberVal(:,1), 2) == 1) & isequal(numberVal(:,1) + 1, numberVal(:,2));

    function fname = getFileNameOnly(fpath)
        [~, name, ext] = fileparts(fpath);
        fname = [name, ext];
    end
end

