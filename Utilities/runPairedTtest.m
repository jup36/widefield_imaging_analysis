function T = runPairedTtest(x, y, contrastLabel)
% Paired two-tailed t-test, one pair per animal. ttest(x,y) with two
% arguments is already the paired form. meanDiff = mean(x - y), so a
% positive value means the first-named condition is larger.
% Callers must pre-mask: a pair with either member missing is unusable.
x = x(:); y = y(:);
assert(numel(x) == numel(y), ...
    'Contrast "%s": paired test needs equal-length vectors -- check the mask.', contrastLabel);

d = x - y;
n = numel(d);

if n < 2
    warning('Contrast "%s": n=%d, too few pairs for a t-test.', contrastLabel, n);
    T = table(string(contrastLabel), n, NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
        'VariableNames', {'contrast','n','df','tstat','p','meanDiff','ci_lo','ci_hi','dz'});
    return;
end

[~, p, ci, st] = ttest(x, y);
dz = mean(d) / std(d);            % Cohen's dz, the paired-design effect size

T = table(string(contrastLabel), n, st.df, st.tstat, p, mean(d), ci(1), ci(2), dz, ...
    'VariableNames', {'contrast','n','df','tstat','p','meanDiff','ci_lo','ci_hi','dz'});

fprintf('  %-26s t(%d) = %7.3f, p = %8.4g, diff = %+.4f [%+.4f %+.4f], dz = %+.2f (n=%d)\n', ...
    contrastLabel, st.df, st.tstat, p, mean(d), ci(1), ci(2), dz, n);
end