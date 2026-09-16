function s = semLocal(v)
n = sum(isfinite(v));
if n < 2
    s = NaN;
else
    s = std(v, 'omitnan') / sqrt(n);
end
end