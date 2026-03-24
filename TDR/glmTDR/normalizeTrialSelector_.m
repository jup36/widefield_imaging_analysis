function idx = normalizeTrialSelector_(x, N, label)
if isempty(x)
    error('%s is empty.', label);
end
if islogical(x)
    if numel(x) ~= N
        error('%s is logical length=%d; expected N=%d.', label, numel(x), N);
    end
    idx = find(x);
else
    idx = x(:);
    if ~isnumeric(idx)
        error('%s must be numeric indices or logical mask.', label);
    end
    if any(~isfinite(idx)) || any(idx~=round(idx))
        error('%s contains non-integer or non-finite entries.', label);
    end
end
idx = idx(:);
if any(idx<1|idx>N)
    error('%s out of bounds for N=%d. min=%d max=%d', label, N, min(idx), max(idx));
end
end