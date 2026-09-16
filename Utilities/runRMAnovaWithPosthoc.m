function [rmRow, phRows] = runRMAnovaWithPosthoc(V, levelNames, conditionLabel)
% One-way repeated-measures ANOVA over the columns of V (animal x level),
% plus Bonferroni-corrected pairwise post-hocs.
%
% V must be complete -- one finite value per animal per level. The caller
% masks for this; fitrm would otherwise drop rows silently and the
% post-hocs would rest on a different sample than the omnibus test.
%
% Returns:
%   rmRow  : one-row table with the omnibus F, uncorrected p, and the
%            Greenhouse-Geisser corrected p and df, plus Mauchly's test.
%   phRows : one row per unordered pair of levels, with the Bonferroni
%            p-value and the mean difference (level1 - level2).
%
% GG is reported because the sphericity assumption is rarely testable with
% real power at these sample sizes. Quote p_GG unless Mauchly clearly
% passes.
[n, nLev] = size(V);
assert(nLev >= 2, 'RM-ANOVA needs at least two levels.');
assert(all(isfinite(V(:))), ...
    'Condition "%s": V contains non-finite values -- mask before calling.', conditionLabel);

varNames = matlab.lang.makeValidName(levelNames);

if n < 3
    warning('Condition "%s": n=%d animals, too few for an RM-ANOVA.', conditionLabel, n);
    rmRow = table(string(conditionLabel), n, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
        'VariableNames', {'condition','n','df1','df2','F','p','df1_GG','df2_GG','p_GG','p_Mauchly'});
    phRows = table();
    return;
end

T  = array2table(V, 'VariableNames', varNames);
WD = table(categorical(levelNames(:), levelNames(:)), 'VariableNames', {'partnerClass'});

rm = fitrm(T, sprintf('%s-%s ~ 1', varNames{1}, varNames{end}), 'WithinDesign', WD);
av = ranova(rm, 'WithinModel', 'partnerClass');

% ranova rows: '(Intercept)', '(Intercept):partnerClass', 'Error(partnerClass)'.
% Match by content rather than a literal, since the exact row name format
% varies by toolbox version.
rowNames = string(av.Properties.RowNames);
iEffect = find(contains(rowNames, 'partnerClass') & ~contains(rowNames, 'Error'), 1);
iError  = find(contains(rowNames, 'Error'), 1);
assert(~isempty(iEffect) && ~isempty(iError), ...
    'Condition "%s": could not locate the partnerClass effect/error rows in ranova output.', conditionLabel);

df1 = av.DF(iEffect);
df2 = av.DF(iError);
F   = av.F(iEffect);
p   = av.pValue(iEffect);
pGG = av.pValueGG(iEffect);

% Greenhouse-Geisser epsilon, for the corrected df actually being used
eps_gg = epsilon(rm);
if istable(eps_gg)
    ggVal = eps_gg.GreenhouseGeisser(1);
else
    ggVal = eps_gg(1);
end

% Mauchly's sphericity test (undefined for two levels -- sphericity is
% automatic there, so report NaN rather than a spurious value)
if nLev > 2
    mw = mauchly(rm);
    pMauchly = mw.pValue(1);
else
    pMauchly = NaN;
end

rmRow = table(string(conditionLabel), n, df1, df2, F, p, ...
    df1*ggVal, df2*ggVal, pGG, pMauchly, ...
    'VariableNames', {'condition','n','df1','df2','F','p','df1_GG','df2_GG','p_GG','p_Mauchly'});

fprintf('\n-- Partner class, %s (n=%d) --\n', conditionLabel, n);
fprintf('  RM-ANOVA: F(%d,%d) = %.3f, p = %.4g | GG-corrected F(%.2f,%.2f), p_GG = %.4g (eps = %.3f)\n', ...
    df1, df2, F, p, df1*ggVal, df2*ggVal, pGG, ggVal);
if ~isnan(pMauchly)
    fprintf('  Mauchly sphericity: p = %.4g%s\n', pMauchly, ...
        ternary_local(pMauchly < 0.05, '  <- sphericity violated, use p_GG', ''));
end

% -------- Bonferroni-corrected pairwise post-hocs --------
mc = multcompare(rm, 'partnerClass', 'ComparisonType', 'bonferroni');

% multcompare returns each pair twice (A-B and B-A). Keep one direction,
% ordered by the level order given by the caller.
lvl1 = string(mc.partnerClass_1);
lvl2 = string(mc.partnerClass_2);
lvlOrder = string(levelNames(:))';
[~, i1] = ismember(lvl1, lvlOrder);
[~, i2] = ismember(lvl2, lvlOrder);
keep = i1 < i2;

phRows = table(repmat(string(conditionLabel), sum(keep), 1), lvl1(keep), lvl2(keep), ...
    mc.Difference(keep), mc.StdErr(keep), mc.pValue(keep), mc.Lower(keep), mc.Upper(keep), ...
    'VariableNames', {'condition','level1','level2','meanDiff','stdErr','p_bonferroni','ci_lo','ci_hi'});

fprintf('  Post-hoc (Bonferroni):\n');
for k = 1:height(phRows)
    fprintf('    %-6s vs %-6s : diff = %+.4f, p = %8.4g, 95%% CI [%+.4f %+.4f]\n', ...
        phRows.level1(k), phRows.level2(k), phRows.meanDiff(k), phRows.p_bonferroni(k), ...
        phRows.ci_lo(k), phRows.ci_hi(k));
end
end