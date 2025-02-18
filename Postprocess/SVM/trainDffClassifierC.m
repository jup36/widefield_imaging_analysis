function [accuracy, betaValues, finalBetaValues] = trainDffClassifierC(dffC1, dffC2, resample, numFolds)
% trainDffClassifierC: Prepares mesoscopic df/f imaging data for SVM classification.
% It organizes the data into a structure suitable for classification.
%
% Inputs:
% - dffC1, dffC2: Cell arrays (1 × numTrials), each entry is a 64×64×numTimeBins matrix.
% - cvFold: Number of cross-validation folds.
%
% Xs: A cell array (1 x # of time bins)
%     - Each time-bin entry has 1 x # of class (e.g., Go vs NoGo) cell array
%     -- Each class entry has 1 x # of trials cell array
%     --- Each trial entry has 1 x # of pixels (e.g., 64x64) array
% y: A cell array (1 x # of class) cell array
%    - Each class entry has # of trial x 1 class labels (e.g., 1: Go, 2:NoGo)

%% Step 1: Validate and Extract Dimensions
dims = cell2mat(cellfun(@size, [dffC1, dffC2], 'UniformOutput', false));

% Ensure consistency across all trials
assert(isscalar(unique(dims(:,1))), 'Mismatch in number of rows!');
assert(isscalar(unique(dims(:,2))), 'Mismatch in number of columns!');
assert(isscalar(unique(dims(:,3))), 'Mismatch in number of time bins!');

% Store consistent dimensions
numRows = unique(dims(:,1));
numCols = unique(dims(:,2));
numTimeBins = unique(dims(:,3));
numPixels = numRows * numCols; % (64 x 64 = 4096)

%% Step 2: Flatten (Vectorize) Data
dffC1_vec = cellfun(@(x) reshape(x, numPixels, []), dffC1, 'UniformOutput', false);
dffC2_vec = cellfun(@(x) reshape(x, numPixels, []), dffC2, 'UniformOutput', false);

numTrialsC1 = length(dffC1_vec);
numTrialsC2 = length(dffC2_vec);

assert(all(cellfun(@(x) size(x,2), dffC1_vec) == numTimeBins), 'Mismatch in time bins for dffC1');
assert(all(cellfun(@(x) size(x,2), dffC2_vec) == numTimeBins), 'Mismatch in time bins for dffC2');

%% Step 3: Organize Data per Time Bin
Xs = cell(1, numTimeBins);

for t = 1:numTimeBins
    % Extract and vectorize trials for this time bin
    X1_bin = cellfun(@(x) x(:,t)', dffC1_vec, 'UniformOutput', false); % Each entry: 1 × numPixels
    X2_bin = cellfun(@(x) x(:,t)', dffC2_vec, 'UniformOutput', false); % Each entry: 1 × numPixels

    % Store separately in a nested cell array
    Xs{t} = {X1_bin, X2_bin}; % Xs{t}{1} -> dffC1 trials, Xs{t}{2} -> dffC2 trials
end

% Construct y (labels) as a separate cell array
y = {ones(numTrialsC1,1), ones(numTrialsC2,1) * 2}; % y{1} = dffC1 labels, y{2} = dffC2 labels

%% Step 4: Identify Valid Pixels within flattened pixels
validMask = ~isnan(sum(cell2mat(dffC1_vec), 2)) | ~isnan(sum(cell2mat(dffC2_vec), 2));
validPixels = find(validMask); % Indices of valid (non-NaN) pixels

disp('Data successfully prepared for SVM classification.');

%% Step 5: Train/Test SVM classifier
[accuracy, betaValues, finalBetaValues] = multiClass_svm_with_cv(Xs, y, resample, numFolds, validPixels);

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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





