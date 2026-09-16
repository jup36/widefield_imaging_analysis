function x = getNumericScalar(x)
% getNumericScalar
%
% Converts numeric or cell-wrapped numeric content to a scalar.
if iscell(x)
    x = cell2mat(x);
end
if isempty(x)
    x = NaN;
else
    x = double(x(1));
end
end