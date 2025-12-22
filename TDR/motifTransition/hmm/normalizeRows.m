function A = normalizeRows(A)
rs = sum(A,2);
rs(~isfinite(rs) | rs<=0) = 1;
A = A ./ rs;
end