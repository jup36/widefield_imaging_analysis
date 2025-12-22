function v = normalizeVec(v)
v = v(:);
s = sum(v);
if s <= 0 || ~isfinite(s)
    v = ones(size(v)) / numel(v);
else
    v = v / s;
end
end
