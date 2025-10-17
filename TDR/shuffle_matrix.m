function Xsh = shuffle_matrix(X, mode, trialIdx, blockSize)
% Apply the same shuffling policy column-wise.
Xsh = X;
for j = 1:size(X,2)
    Xsh(:,j) = shuffle_column(X(:,j), mode, trialIdx, blockSize);
end
end