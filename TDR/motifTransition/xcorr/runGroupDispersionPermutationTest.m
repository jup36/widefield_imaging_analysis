function dispOut = runGroupDispersionPermutationTest(Y, sessInfo, groupDefs, varargin)
%RUNGROUPDISPERSIONPERMUTATIONTEST
%   Tests whether one group's animals are more DISPERSED (in MDS space)
%   than another's -- e.g. "slow learners' animals occupy more separate,
%   scattered regions than fast learners' do." Unlike the convergence/
%   ellipsoid analyses elsewhere in this project, this metric does NOT
%   build any reference/centroid from one group and test that SAME group
%   against it -- so there is no circularity concern here at all; a
%   simple group-label permutation test is valid and sufficient.
%
%   Dispersion is decomposed additively (same logic as a one-way ANOVA):
%     TOTAL dispersion (within a group) = BETWEEN-animal + WITHIN-animal
%   where, pooling ALL sessions of ALL animals in a group into one cloud:
%     - BETWEEN-animal: how far each animal's OWN mean (across its
%       sessions) sits from the group's grand mean -- "do the animals
%       themselves occupy different regions?"
%     - WITHIN-animal: how far each animal's individual SESSIONS scatter
%       around that animal's OWN mean -- "is any one animal's trajectory
%       itself erratic?"
%   These are reported SEPARATELY because they can point in different
%   directions -- e.g. animals could sit in very different regions
%   (high between) while each individually has a tight, low-noise
%   trajectory (low within), or vice versa.
%
%   dispOut = runGroupDispersionPermutationTest(Y, sessInfo, groupDefs, ...)
%
% INPUTS
%   Y         : [S x dimFull] MDS coordinates for ONE stream.
%   sessInfo  : matching sessInfo table (needs mouseId, sessWithin).
%   groupDefs : scalar struct, exactly 2 groups, e.g.
%                 groupDefs.fast = {'m1044','m1045','m1092','m1094'};
%                 groupDefs.slow = {'m1048','m1049','m1613','m1859','m1873'};
%
% NAME-VALUE ARGS
%   'dims'       : which columns of Y to use. Default: 1:size(Y,2).
%   'nFinalSess' : if provided, restricts to each animal's LAST N sessions
%                  only (matching the convergence analyses' late-session
%                  focus). Default: [] (use ALL sessions -- the more
%                  natural choice for a trajectory-wide dispersion
%                  question, as opposed to a late-session-only one).
%   'direction'  : 'slow_minus_fast' (default) or 'fast_minus_slow' --
%                  which group's dispersion is hypothesized to be LARGER,
%                  for the one-sided permutation p-value. Set based on
%                  your actual hypothesis (here: slow learners more
%                  dispersed), not chosen post-hoc to make a p-value
%                  smaller.
%   'nPerm'      : 1000 (default). IGNORED if exact enumeration is used
%                  (see below) -- kept only as a fallback/cap.
%   'useExact'   : true (default). If nchoosek(nAllMice, nGroup1) is
%                  small enough to be tractable (<= 5000 by default),
%                  enumerates ALL possible group partitions exactly
%                  instead of Monte Carlo sampling -- gives an exact,
%                  noise-free null and p-value, and is typically FASTER
%                  than a large Monte Carlo run for small cohorts (e.g.
%                  this project's 9-mouse, 4-vs-5 split: C(9,4)=126).
%                  Falls back to Monte Carlo with 'nPerm' draws if the
%                  exact count exceeds 'maxExact'.
%   'maxExact'   : 5000 (default) -- tractability cutoff for 'useExact'.
%   'permSeed'   : [] (default) or scalar seed (only relevant if Monte
%                  Carlo fallback is used).
%   'verbose'    : true (default).
%   'doSave'     : true (default). Requires 'saveDir' if true.
%   'saveDir', 'saveTag' : as in the other stats functions in this project.
%
% OUTPUT (dispOut)
%   .obsBetween, .obsWithin, .obsTotal   : observed TRUE-group values
%                  (group2 - group1, per 'direction'; group order is
%                  fieldnames(groupDefs) order, so 'slow_minus_fast'
%                  assumes groupDefs.fast is listed first)
%   .nullBetween, .nullWithin, .nullTotal : null distributions (exact or MC)
%   .pBetween, .pWithin, .pTotal          : one-sided empirical p-values
%   .nPartitions, .isExact, .direction, .dims, .nFinalSess
%
% See also: runConvergenceEllipsoidPermutationTest, runConvergenceRMANOVA_fastSlow

%% -------------------- parse options --------------------
p = inputParser;
p.addParameter('dims', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
p.addParameter('nFinalSess', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>=1));
p.addParameter('direction', 'slow_minus_fast', @(s) any(strcmpi(string(s), ["slow_minus_fast","fast_minus_slow"])));
p.addParameter('nPerm', 1000, @(x) isnumeric(x) && isscalar(x) && x>=1 && x==round(x));
p.addParameter('useExact', true, @(x) islogical(x) && isscalar(x));
p.addParameter('maxExact', 5000, @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('permSeed', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('verbose', true, @(x) islogical(x) && isscalar(x));
p.addParameter('doSave', true, @(x) islogical(x) && isscalar(x));
p.addParameter('saveDir', '', @(s) ischar(s) || isstring(s));
p.addParameter('saveTag', '', @(s) ischar(s) || isstring(s));
p.parse(varargin{:});
opt = p.Results;

if ~isempty(opt.permSeed)
    rng(opt.permSeed);
end

%% -------------------- validate + build fixed mouse universe --------------------
groupNames = fieldnames(groupDefs);
assert(numel(groupNames) == 2, 'runGroupDispersionPermutationTest requires exactly 2 groups.');

group1IdC = groupDefs.(groupNames{1});   % e.g. 'fast', listed first
group2IdC = groupDefs.(groupNames{2});   % e.g. 'slow'
overlap = intersect(group1IdC, group2IdC);
assert(isempty(overlap), 'Animal(s) %s appear in both groups.', strjoin(overlap, ', '));

allMiceIdC = [group1IdC(:)', group2IdC(:)'];
nAllMice   = numel(allMiceIdC);
n1         = numel(group1IdC);   % size of group listed first (e.g. fast)

assert(all(ismember({'mouseId','sessWithin'}, sessInfo.Properties.VariableNames)), ...
    'sessInfo must have mouseId and sessWithin columns.');

mouseId    = string(sessInfo.mouseId);
sessWithin = sessInfo.sessWithin;

missingMice = allMiceIdC(~ismember(string(allMiceIdC), unique(mouseId)));
if ~isempty(missingMice)
    error('The following animal(s) from groupDefs were not found in sessInfo.mouseId: %s', ...
        strjoin(missingMice, ', '));
end

if isempty(opt.dims)
    dims = 1:size(Y,2);
else
    dims = opt.dims(:)';
end
Yd = Y(:, dims);

%% -------------------- optionally restrict to each animal's last N sessions --------------------
keepRow = true(numel(mouseId), 1);
if ~isempty(opt.nFinalSess)
    for im = 1:nAllMice
        rowsIm = find(mouseId == allMiceIdC{im});
        maxSess = max(sessWithin(rowsIm));
        keepRow(rowsIm) = keepRow(rowsIm) & (sessWithin(rowsIm) > maxSess - opt.nFinalSess);
    end
end

mouseIdUse = mouseId(keepRow);
Yuse       = Yd(keepRow, :);

%% -------------------- core reusable subroutine --------------------
    function [totalDisp, betweenDisp, withinDisp] = local_dispersion(memberIdC)
        rowsG = ismember(mouseIdUse, string(memberIdC));
        Xg = Yuse(rowsG, :);
        gVec = mouseIdUse(rowsG);

        muG = mean(Xg, 1);
        totalDisp = mean(sum((Xg - muG).^2, 2));

        uAnimals = unique(gVec, 'stable');
        d2between = [];
        d2within  = [];
        for a = 1:numel(uAnimals)
            idxA = (gVec == uAnimals(a));
            Xa = Xg(idxA, :);
            muA = mean(Xa, 1);
            d2between = [d2between; repmat(sum((muA - muG).^2), sum(idxA), 1)]; %#ok<AGROW>
            d2within  = [d2within;  sum((Xa - muA).^2, 2)]; %#ok<AGROW>
        end
        betweenDisp = mean(d2between);
        withinDisp  = mean(d2within);
    end

    function [dTotal, dBetween, dWithin] = local_groupDiff(members1, members2)
        [t1, b1, w1] = local_dispersion(members1);
        [t2, b2, w2] = local_dispersion(members2);
        switch lower(string(opt.direction))
            case "slow_minus_fast"
                dTotal   = t2 - t1;   % group listed second minus group listed first
                dBetween = b2 - b1;
                dWithin  = w2 - w1;
            case "fast_minus_slow"
                dTotal   = t1 - t2;
                dBetween = b1 - b2;
                dWithin  = w1 - w2;
        end
    end

%% -------------------- observed --------------------
[obsTotal, obsBetween, obsWithin] = local_groupDiff(group1IdC, group2IdC);

fprintf('\n============================================================\n');
fprintf('Group dispersion permutation test (direction: %s)\n', opt.direction);
fprintf('Observed diff -- total: %.4f, between-animal: %.4f, within-animal: %.4f\n', ...
    obsTotal, obsBetween, obsWithin);

%% -------------------- null: exact enumeration or Monte Carlo --------------------
nPartitions = nchoosek(nAllMice, n1);
isExact = opt.useExact && (nPartitions <= opt.maxExact);

if isExact
    fprintf('Using EXACT enumeration: all C(%d,%d) = %d partitions.\n', nAllMice, n1, nPartitions);
    combos = nchoosek(1:nAllMice, n1);
    nDraws = size(combos, 1);
else
    fprintf('Using Monte Carlo: %d random partitions (exact count %d exceeds maxExact=%d).\n', ...
        opt.nPerm, nPartitions, opt.maxExact);
    nDraws = opt.nPerm;
end

nullTotal   = nan(nDraws, 1);
nullBetween = nan(nDraws, 1);
nullWithin  = nan(nDraws, 1);

tPerm = tic;
for pI = 1:nDraws
    if isExact
        idx1 = combos(pI, :);
    else
        permOrder = randperm(nAllMice);
        idx1 = permOrder(1:n1);
    end
    idx2 = setdiff(1:nAllMice, idx1);

    members1 = allMiceIdC(idx1);
    members2 = allMiceIdC(idx2);

    [nullTotal(pI), nullBetween(pI), nullWithin(pI)] = local_groupDiff(members1, members2);

    if opt.verbose && (mod(pI, 500) == 0 || pI == nDraws)
        fprintf('  %d/%d (%.1fs elapsed)\n', pI, nDraws, toc(tPerm));
    end
end

pTotal   = (1 + sum(nullTotal   >= obsTotal))   / (nDraws + 1);
pBetween = (1 + sum(nullBetween >= obsBetween)) / (nDraws + 1);
pWithin  = (1 + sum(nullWithin  >= obsWithin))  / (nDraws + 1);

fprintf('\nCompleted in %.1f sec.\n', toc(tPerm));
fprintf('  TOTAL          diff=%.4f | null mean=%.4f SD=%.4f | p=%.4f\n', obsTotal,   mean(nullTotal,'omitnan'),   std(nullTotal,'omitnan'),   pTotal);
fprintf('  BETWEEN-animal diff=%.4f | null mean=%.4f SD=%.4f | p=%.4f\n', obsBetween, mean(nullBetween,'omitnan'), std(nullBetween,'omitnan'), pBetween);
fprintf('  WITHIN-animal  diff=%.4f | null mean=%.4f SD=%.4f | p=%.4f\n', obsWithin,  mean(nullWithin,'omitnan'),  std(nullWithin,'omitnan'),  pWithin);
fprintf('============================================================\n');

%% -------------------- package output --------------------
dispOut = struct();
dispOut.obsTotal    = obsTotal;    dispOut.nullTotal   = nullTotal;   dispOut.pTotal   = pTotal;
dispOut.obsBetween  = obsBetween;  dispOut.nullBetween = nullBetween; dispOut.pBetween = pBetween;
dispOut.obsWithin   = obsWithin;   dispOut.nullWithin  = nullWithin;  dispOut.pWithin  = pWithin;
dispOut.nPartitions = nPartitions;
dispOut.isExact     = isExact;
dispOut.direction   = opt.direction;
dispOut.dims        = dims;
dispOut.nFinalSess  = opt.nFinalSess;
dispOut.group1IdC   = group1IdC;
dispOut.group2IdC   = group2IdC;

%% -------------------- optional save --------------------
if opt.doSave
    if strlength(strtrim(string(opt.saveDir))) == 0
        error('doSave is true but no ''saveDir'' was provided.');
    end
    outDir = char(string(opt.saveDir));
    if exist(outDir, 'dir') ~= 7
        mkdir(outDir);
    end
    tag = strtrim(char(string(opt.saveTag)));
    if isempty(tag)
        tag = 'unlabeled';
        warning('runGroupDispersionPermutationTest:noSaveTag', 'No ''saveTag'' provided -- using "unlabeled".');
    end
    dateStr  = char(datetime('today','Format','MMddyy'));
    saveName = sprintf('groupDispersionPermTest_%s_%s.mat', tag, dateStr);
    saveFullPath = fullfile(outDir, saveName);
    save(saveFullPath, 'dispOut', 'groupDefs');
    dispOut.save = struct('didSave', true, 'file', saveFullPath);
    fprintf('\nSaved to:\n%s\n', saveFullPath);
else
    dispOut.save = struct('didSave', false, 'file', '');
end

end