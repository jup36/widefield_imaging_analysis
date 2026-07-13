function [contrastTable, statsOut, h] = computeMotifClusterContrast_byGroup(csC_perAnimal, animalIDs, groupDefs, clusterA_idx, clusterB_idx, varargin)
%COMPUTEMOTIFCLUSTERCONTRAST_BYGROUP  Animal-level between/within-cluster cosine similarity, compared across groups.
%
% SYNOPSIS
%   [contrastTable, statsOut, h] = computeMotifClusterContrast_byGroup( ...
%       csC_perAnimal, animalIDs, groupDefs, clusterA_idx, clusterB_idx, ...)
%
% DESCRIPTION
%   Quantifies the "between-cluster contrast" visible in a beta cosine
%   similarity heatmap (e.g. a top-left supercluster vs. a bottom-right
%   cluster) at the ANIMAL level, so it can actually be statistically
%   compared across groups (e.g. fast vs. slow learners) rather than just
%   eyeballing the two groups' mean heatmaps.
%
%   For each animal i, using that animal's own session-averaged matrix
%   csC_perAnimal{i} (NOT a group-mean matrix):
%     betweenClusterMean(i) = mean of all cosine similarities between
%                              clusterA_idx and clusterB_idx motifs
%     withinA(i), withinB(i) = mean off-diagonal cosine similarity WITHIN
%                              clusterA_idx / clusterB_idx respectively
%                              (context metrics -- do within-cluster
%                              similarities differ too, or is it
%                              specifically the between-cluster term?)
%     segregationIndex(i)    = mean([withinA(i) withinB(i)]) - betweenClusterMean(i)
%                              (a single composite "modularity-like" score,
%                              provided for convenience/exploration)
%
%   These are collected into contrastTable (one row per animal), then
%   betweenClusterMean is compared between two specified groups (default:
%   the first two fields of groupDefs) via a dot plot with group means +/-
%   SEM, plus a rank-sum test and a label-permutation test (more
%   appropriate than a t-test given the very small n typical of this kind
%   of dataset).
%
% IMPORTANT CAVEAT
%   clusterA_idx/clusterB_idx are typically read off a heatmap that was
%   reordered by clustering the GRAND MEAN (i.e. pooling all groups). If
%   one group disproportionately drives that clustering, the cluster
%   boundaries themselves are not fully independent of the group
%   difference you are about to test. This function does not correct for
%   that -- it only computes and compares the contrast given clusters YOU
%   supply. Consider validating cluster stability independently (e.g. by
%   re-deriving clusters from one group only, or leave-one-group-out)
%   before treating the resulting p-value as fully independent evidence.
%
% INPUTS
%   csC_perAnimal : [nAnimals x 1] cell array of [K x K] animal-level mean
%                   cosine similarity matrices (as produced by
%                   analyzeMotifBetaCosineSimilarity_acrossAnimals.m).
%                   Natural (unreordered) motif indexing.
%   animalIDs     : [nAnimals x 1] cellstr of animal IDs, same order as
%                   csC_perAnimal.
%   groupDefs     : scalar struct, field name = group label, value =
%                   cellstr of animal IDs in that group (e.g.
%                   groupDefs.fast = {...}; groupDefs.slow = {...};).
%                   Must have at least 2 fields.
%   clusterA_idx  : numeric vector of motif indices (1..K, NATURAL motif
%                   numbering -- e.g. the actual motif numbers you'd read
%                   off a reordered heatmap's tick labels, not their
%                   display position) belonging to cluster A.
%   clusterB_idx  : numeric vector of motif indices belonging to cluster B.
%                   Must not overlap with clusterA_idx.
%
% NAME-VALUE ARGS
%   'groupA'        : name of first group to compare (must be a field of
%                      groupDefs). Default: fieldnames(groupDefs){1}.
%   'groupB'        : name of second group to compare. Default:
%                      fieldnames(groupDefs){2}.
%   'nPerm'          : number of permutations for the label-permutation
%                      test. Default: 10000.
%   'rngSeed'        : seed for the permutation test's random stream
%                      (uses a local RandStream -- does not alter your
%                      global rng state). Default: 1.
%   'groupColors'    : scalar struct, field name = group label, value =
%                      1x3 RGB, to color the dot plot. Any group not
%                      listed falls back to a built-in default (green for
%                      'fast', purple for 'slow', gray otherwise).
%   'figureScaleFactor' : figure size multiplier. Default: 1.2.
%   'visible'        : 'on'/'off'. Default: 'on'.
%   'figSaveDir', 'figSaveKeyword' : save-to-PDF options, same convention
%                      as other functions in this pipeline.
%
% OUTPUTS
%   contrastTable : table, one row per animal, columns: animalID, group,
%                   betweenClusterMean, withinA, withinB, segregationIndex.
%   statsOut      : struct with groupA/groupB names, n, mean, SEM per
%                   group, and pRanksum / pPermutation / obsDiff / nullDiff.
%   h             : struct with fig/ax handles and opt.
%
% EXAMPLE
%   groupDefs.fast = fast_learners;
%   groupDefs.slow = slow_learners;
%   clusterA_idx = [19 22 24 4 8 1 10 23 3 5 17 9 20];
%   clusterB_idx = [2 14 11 6 15 12];
%   [contrastTable, statsOut, h] = computeMotifClusterContrast_byGroup( ...
%       csC_perAnimal, animalIDs, groupDefs, clusterA_idx, clusterB_idx, ...
%       'figSaveDir', figSaveDir, 'figSaveKeyword', 'clusterContrast');
%
% See also: analyzeMotifBetaCosineSimilarity_acrossAnimals, plot_cosineSimilarity_groupDifference

% -------- parse args --------
grpNamesAll = fieldnames(groupDefs);
assert(numel(grpNamesAll) >= 2, 'groupDefs must have at least 2 group fields.');

p = inputParser;
p.addParameter('groupA', grpNamesAll{1}, @(s) ischar(s) || isstring(s));
p.addParameter('groupB', grpNamesAll{2}, @(s) ischar(s) || isstring(s));
p.addParameter('nPerm', 10000, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('rngSeed', 1, @(x) isnumeric(x) && isscalar(x));
p.addParameter('groupColors', struct(), @(x) isstruct(x) && isscalar(x));
p.addParameter('figureScaleFactor', 1.2, @(x) isnumeric(x) && x>0);
p.addParameter('visible', 'on', @(s) any(strcmpi(s, {'on','off'})));
p.addParameter('figSaveDir', {}, @(x) isempty(x) || ischar(x) || isstring(x) || iscell(x));
p.addParameter('figSaveKeyword', '', @(s) ischar(s) || isstring(s));
p.parse(varargin{:});
opt = p.Results;
opt.groupA = char(opt.groupA);
opt.groupB = char(opt.groupB);

assert(isfield(groupDefs, opt.groupA), 'groupDefs has no field "%s".', opt.groupA);
assert(isfield(groupDefs, opt.groupB), 'groupDefs has no field "%s".', opt.groupB);

% normalize figSaveDir
figSaveDir = opt.figSaveDir;
if iscell(figSaveDir)
    if isempty(figSaveDir), figSaveDir = ""; else, figSaveDir = string(figSaveDir{1}); end
else
    figSaveDir = string(figSaveDir);
end
figSaveKeyword = string(opt.figSaveKeyword);

% -------- sanity checks --------
nAnimals = numel(csC_perAnimal);
assert(numel(animalIDs)==nAnimals, 'animalIDs must match csC_perAnimal length.');

firstValid = find(~cellfun(@isempty, csC_perAnimal), 1, 'first');
assert(~isempty(firstValid), 'csC_perAnimal has no valid (non-empty) entries.');
K = size(csC_perAnimal{firstValid}, 1);

clusterA_idx = clusterA_idx(:)';
clusterB_idx = clusterB_idx(:)';
assert(all(clusterA_idx>=1 & clusterA_idx<=K), 'clusterA_idx must be within 1..%d.', K);
assert(all(clusterB_idx>=1 & clusterB_idx<=K), 'clusterB_idx must be within 1..%d.', K);
assert(isempty(intersect(clusterA_idx, clusterB_idx)), 'clusterA_idx and clusterB_idx must not overlap.');

% -------- 1) per-animal between/within cluster means --------
betweenClusterMean = nan(nAnimals,1);
withinA = nan(nAnimals,1);
withinB = nan(nAnimals,1);

for i = 1:nAnimals
    Ci = csC_perAnimal{i};
    if isempty(Ci)
        continue;
    end
    subAB = Ci(clusterA_idx, clusterB_idx);
    betweenClusterMean(i) = mean(subAB(:));
    withinA(i) = offdiag_mean_(Ci(clusterA_idx, clusterA_idx));
    withinB(i) = offdiag_mean_(Ci(clusterB_idx, clusterB_idx));
end
segregationIndex = mean([withinA, withinB], 2) - betweenClusterMean;

% -------- 2) group assignment + table --------
group = repmat({'unassigned'}, nAnimals, 1);
for gi = 1:numel(grpNamesAll)
    memberIDs = groupDefs.(grpNamesAll{gi});
    group(ismember(animalIDs, memberIDs)) = grpNamesAll(gi);
end

contrastTable = table(animalIDs(:), group(:), betweenClusterMean, withinA, withinB, segregationIndex, ...
    'VariableNames', {'animalID','group','betweenClusterMean','withinA','withinB','segregationIndex'});

% -------- 3) statistics: groupA vs groupB on betweenClusterMean --------
xA = contrastTable.betweenClusterMean(strcmp(contrastTable.group, opt.groupA) & ~isnan(contrastTable.betweenClusterMean));
xB = contrastTable.betweenClusterMean(strcmp(contrastTable.group, opt.groupB) & ~isnan(contrastTable.betweenClusterMean));

statsOut = struct();
statsOut.groupA = opt.groupA;
statsOut.groupB = opt.groupB;
statsOut.nA = numel(xA);
statsOut.nB = numel(xB);
statsOut.meanA = mean(xA);
statsOut.meanB = mean(xB);
statsOut.semA = std(xA) / sqrt(max(numel(xA),1));
statsOut.semB = std(xB) / sqrt(max(numel(xB),1));

if numel(xA) >= 1 && numel(xB) >= 1
    try
        statsOut.pRanksum = ranksum(xA, xB);
    catch ME
        statsOut.pRanksum = NaN;
        warning('computeMotifClusterContrast_byGroup:ranksumFailed', 'ranksum failed: %s', ME.message);
    end
    [statsOut.pPermutation, statsOut.obsDiff, statsOut.nullDiff] = ...
        permutation_test_twoSample_(xA, xB, opt.nPerm, opt.rngSeed);
else
    statsOut.pRanksum = NaN;
    statsOut.pPermutation = NaN;
    statsOut.obsDiff = NaN;
    statsOut.nullDiff = [];
    warning('computeMotifClusterContrast_byGroup:insufficientData', ...
        'Not enough data in "%s" and/or "%s" to run statistics.', opt.groupA, opt.groupB);
end

% -------- 4) plot --------
h = struct();
h.opt = opt;

h.fig = figure('Color','w', 'Visible', opt.visible);
set(h.fig, 'Units', 'normalized');
pos = get(h.fig, 'Position');
pos(3:4) = pos(3:4) * opt.figureScaleFactor;
set(h.fig, 'Position', pos);
h.ax = axes('Parent', h.fig); hold(h.ax, 'on');

colorA = resolve_group_color_(opt.groupA, opt.groupColors);
colorB = resolve_group_color_(opt.groupB, opt.groupColors);

jitterA = (rand(numel(xA),1) - 0.5) * 0.15;
jitterB = (rand(numel(xB),1) - 0.5) * 0.15;

scatter(h.ax, 1 + jitterA, xA, 60, 'MarkerFaceColor', colorA, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
scatter(h.ax, 2 + jitterB, xB, 60, 'MarkerFaceColor', colorB, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);

errorbar(h.ax, 1, statsOut.meanA, statsOut.semA, 'o', ...
    'Color', colorA*0.6, 'MarkerFaceColor', colorA, 'MarkerSize', 10, 'LineWidth', 2, 'CapSize', 12);
errorbar(h.ax, 2, statsOut.meanB, statsOut.semB, 'o', ...
    'Color', colorB*0.6, 'MarkerFaceColor', colorB, 'MarkerSize', 10, 'LineWidth', 2, 'CapSize', 12);

xlim(h.ax, [0.5 2.5]);
plot(h.ax, xlim(h.ax), [0 0], 'k--', 'LineWidth', 0.8);   % zero reference (anti-correlation threshold)

xticks(h.ax, [1 2]);
xticklabels(h.ax, {sprintf('%s (n=%d)', opt.groupA, statsOut.nA), sprintf('%s (n=%d)', opt.groupB, statsOut.nB)});
set(h.ax, 'TickLabelInterpreter', 'none');
ylabel(h.ax, 'Between-cluster mean cosine similarity');
title(h.ax, 'Cluster A \times Cluster B contrast', 'Interpreter', 'tex');

yl = ylim(h.ax);
text(h.ax, 1.5, yl(2) - 0.05*(yl(2)-yl(1)), ...
    sprintf('rank-sum p = %.3g\npermutation p = %.3g', statsOut.pRanksum, statsOut.pPermutation), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

box(h.ax, 'off');
grid(h.ax, 'on');
set(h.ax, 'TickDir', 'out');

% -------- save (optional) --------
if strlength(figSaveDir) > 0
    if ~isfolder(figSaveDir)
        mkdir(figSaveDir);
    end
    dateStr = char(datetime("today","Format","MMddyy"));
    parts = strings(0,1);
    parts(end+1,1) = "clusterContrast";
    if strlength(figSaveKeyword) > 0, parts(end+1,1) = figSaveKeyword; end
    parts(end+1,1) = dateStr;
    figSaveName = strjoin(parts, "_");
    print(h.fig, fullfile(figSaveDir, figSaveName), '-dpdf', '-painters', '-bestfit');
end

end % function


% ===== helper: mean of off-diagonal entries of a square submatrix =====
function m = offdiag_mean_(M)
n = size(M,1);
if n <= 1
    m = NaN;   % a single-motif "cluster" has no within-cluster pairs
    return;
end
mask = ~eye(n);
m = mean(M(mask));
end


% ===== helper: label-permutation two-sample test (uses a local RandStream --
%              does not alter the caller's global rng state) =====
function [pVal, obsDiff, nullDiff] = permutation_test_twoSample_(x, y, nPerm, seed)
x = x(:); y = y(:);
obsDiff = mean(x) - mean(y);
pooled = [x; y];
nX = numel(x);
n = numel(pooled);

s = RandStream('mt19937ar', 'Seed', seed);
nullDiff = nan(nPerm,1);
for pI = 1:nPerm
    idx = randperm(s, n);
    xp = pooled(idx(1:nX));
    yp = pooled(idx(nX+1:end));
    nullDiff(pI) = mean(xp) - mean(yp);
end
% +1 correction (Davison & Hinkley convention): avoids reporting p=0
pVal = (sum(abs(nullDiff) >= abs(obsDiff)) + 1) / (nPerm + 1);
end


% ===== helper: resolve a group's plot color, with sensible fallbacks =====
function c = resolve_group_color_(gName, groupColors)
if isstruct(groupColors) && isfield(groupColors, gName)
    c = groupColors.(gName);
    return;
end
switch lower(gName)
    case 'fast'
        c = hsv2rgb([0.36 0.70 0.65]);
    case 'slow'
        c = hsv2rgb([0.80 0.70 0.65]);
    otherwise
        c = [0.5 0.5 0.5];
end
end