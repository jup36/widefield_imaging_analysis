function vals = getNumericFieldColumn(s, fieldName)
% getNumericFieldColumn
%
% Extracts a scalar numeric field from a struct array as a column vector.
vals = NaN(numel(s), 1);
for i = 1:numel(s)
    if ~isfield(s, fieldName)
        continue;
    end
    v = s(i).(fieldName);
    if iscell(v)
        v = cell2mat(v);
    end
    if isempty(v)
        continue;
    end
    vals(i) = double(v(1));
end
end
