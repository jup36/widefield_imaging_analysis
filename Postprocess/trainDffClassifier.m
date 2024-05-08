function [svmRez, svmRezTs, nbRez, nbRezTs, X1, X2] = trainDffClassifier(dffC1, dffC2, timepts, binSize, stepSize, cvFold)

[X1, ~] = binAvg1msSpkCountMat(cell2mat(dffC1), binSize, stepSize);
[X2, ~] = binAvg1msSpkCountMat(cell2mat(dffC2), binSize, stepSize);

X_bins = timepts(1:stepSize:end);

if length(X_bins)-size(X1, 2)==1
    X_bins = X_bins(1:end-1)+(stepSize/1000)/2;
end

assert(size(X1, 2)==size(X2, 2))

% train and test svm go vs nogo
Xs = [X1; X2];
y = [ones(size(X1, 1), 1)*1; ones(size(X2, 1), 1)*2];
svmRez = multiClass_svm_peth(Xs, y, cvFold);
svmRezTs = X_bins;
nbRez = multiClass_naiveBayesClassifier_peth(Xs, y, cvFold);
nbRezTs = X_bins;

end
