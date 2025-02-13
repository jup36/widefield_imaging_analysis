function [accuracy, betaValues, finalBetaValues] = multiClass_svm_with_cv(Xs, y, numResample, numFolds, validPixels)
% multiClass_svm_with_cv: Runs resampling and 5-fold cross-validation SVM classification.
%
% Inputs:
% - Xs: 1×numTimeBins cell array containing feature matrices for each time bin.
% - y: Labels for each class.
% - numResample: Number of resampling iterations for balancing class sizes.
% - numFolds: Number of cross-validation folds (e.g., 5).
% - validPixels: Indices of valid pixels (non-NaN).
%
% Outputs:
% - accuracy: numResample × numFolds × numTimeBins matrix of classification accuracies.
% - betaValues: SVM feature weights mapped back to 64×64.
% - finalBetaValues: betaValues averaged across folds and resamples.

%% **Step 1: Check if Resampling is Needed**
if isscalar(unique(cellfun(@length, y)))
    numResample = 1; % If equal number of samples per class, do not resample
    doResample = false;
else
    doResample = true;
end

accuracy = NaN(numResample, numFolds, length(Xs)); % Store accuracy across resampling & CV folds
betaValues = cell(numResample, numFolds, length(Xs)); % Store beta values

for rs = 1:numResample  % Outer loop: Resampling to balance class sizes
    for t = 1:length(Xs)  % Loop over time bins
        X1_bin = cell2mat(Xs{t}{1}); % Flatten trials for class 1
        X2_bin = cell2mat(Xs{t}{2}); % Flatten trials for class 2

        % Ensure proper reshaping: (numTrials × numPixels)
        X1_bin = reshape(X1_bin, 64*64, []).'; % Convert to (numTrials × numPixels)
        X2_bin = reshape(X2_bin, 64*64, []).'; % Convert to (numTrials × numPixels)

        % Apply valid pixel selection
        X1_bin = X1_bin(:, validPixels); % numTrials × numValidPixels
        X2_bin = X2_bin(:, validPixels); % numTrials × numValidPixels

        % Combine feature matrix and labels
        X = [X1_bin; X2_bin];
        y_combined = [y{1}; y{2}];

        % **Step 2: Perform Resampling (Only if Needed)**
        if doResample
            minSize = min(cellfun(@length, y)); % Smallest class size
            X_balanced = [];
            y_balanced = [];

            uniqueClasses = unique(y_combined);
            for i = 1:length(uniqueClasses)
                classIndices = find(y_combined == uniqueClasses(i));
                randIndices = randsample(classIndices, minSize);
                X_balanced = [X_balanced; X(randIndices, :)];
                y_balanced = [y_balanced; y_combined(randIndices)];
            end
            X = X_balanced;
            y_combined = y_balanced;
        end

        % **Step 3: Perform 5-Fold Cross-Validation**
        cv = cvpartition(length(y_combined), 'KFold', numFolds);

        for fold = 1:numFolds
            % Define train/test indices
            trainIdx = training(cv, fold);
            testIdx = test(cv, fold);

            % Split data into training and testing sets
            XTrain = X(trainIdx, :);
            YTrain = y_combined(trainIdx);
            XTest  = X(testIdx, :);
            YTest  = y_combined(testIdx);

            % Train SVM
            SVMModel = fitcecoc(XTrain, YTrain); % By default, fitcecoc uses linear SVMs if a kernel is not specified

            % **Step 4: Extract Beta Values and Map Back to 64x64 Grid**
            % Extract beta values from the trained SVM
            if isprop(SVMModel, 'BinaryLearners')
                betaTemp = cellfun(@(m) m.Beta, SVMModel.BinaryLearners, 'UniformOutput', false);
                betaMatrix = mean(cell2mat(betaTemp), 2); % Average across binary classifiers if needed

                % Initialize full 64×64 matrix with NaNs
                betaFull = NaN(64*64, 1); % Vector of all pixels, set to NaN by default
                betaFull(validPixels) = betaMatrix; % Insert valid beta values at correct indices

                % Reshape into 64×64 grid
                betaMap = reshape(betaFull, 64, 64);

                % Store the remapped beta values
                betaValues{rs, fold, t} = betaMap;
            else
                betaValues{rs, fold, t} = NaN(64, 64); % Handle cases where beta extraction fails
            end


            % **Step 5: Evaluate Model Performance**
            YPred = predict(SVMModel, XTest);
            accuracy(rs, fold, t) = sum(YPred == YTest) / length(YTest);
        end
        fprintf('Completed Resample %d/%d, Time Bin %d/%d.\n', rs, numResample, t, length(Xs));
    end
end

% Initialize final averaged betaValues
finalBetaValues = cell(1, length(Xs)); % 1xnumTimeBins cell array

for t = 1:length(Xs) % Loop over time bins
    betaStack = []; % Placeholder for stacking beta values

    for rs = 1:numResample
        for fold = 1:numFolds
            betaCurrent = abs(betaValues{rs, fold, t}); % Extract current beta map
            betaStack = cat(3, betaStack, betaCurrent); % Stack across 3rd dimension
        end
    end

    % Compute the mean across resamples and folds, ignoring NaNs
    if ~isempty(betaStack)
        finalBetaValues{t} = nanmean(betaStack, 3);
    else
        finalBetaValues{t} = NaN(64, 64); % In case no valid beta maps exist
    end
end

disp('✅ Resampling & 5-Fold Cross-Validation SVM training completed.');
end