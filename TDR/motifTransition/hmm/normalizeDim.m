function Y = normalizeDim(X, dim)
Y = X ./ (sum(X, dim) + eps);
end