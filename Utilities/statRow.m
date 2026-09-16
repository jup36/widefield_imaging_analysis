function row = statRow(statsT, contrastLabel)
row = statsT(strcmp(statsT.contrast, contrastLabel), :);
assert(height(row) == 1, 'Expected exactly one row for contrast "%s".', contrastLabel);
end