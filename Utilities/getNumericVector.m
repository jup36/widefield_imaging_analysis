function v = getNumericVector(x)
% getNumericVector
%
% Converts numeric or cell-wrapped numeric content to a row vector.
if iscell(x)
    x = cell2mat(x);
end
if isempty(x)
    v = [];
else
    v = double(x(:)');
end
end