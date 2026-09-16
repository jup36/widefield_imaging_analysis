function amongOut = runAmongDistancePermutationTest_blocked(Y, sessInfo, groupDefs, varargin)
%RUNAMONGDISTANCEPERMUTATIONTEST_BLOCKED
%   Same exact permutation test as runAmongDistancePermutationTest, but
%   with a different preprocessing step that AVOIDS discarding any
%   session data. Instead of right-aligning by each animal's final
%   session and clipping everyone to the SHORTEST animal's depth (which
%   throws away early sessions from longer-recorded animals -- e.g. a
%   17-session animal would lose its first 9 sessions if the shortest
%   animal only has 8), this function splits EACH animal's OWN sessions
%   into three RELATIVE, contiguous blocks -- Early / Middle / Late
%   thirds of THAT animal's own training -- and averages the MDS points
%   within each block. Every session from every animal contributes to
%   exactly one block's average; nothing is dropped.
%
%   Concretely, for an animal with N sessions, sessions 1..N are split
%   into 3 as-equal-as-possible contiguous chunks (e.g. N=8 -> chunks of
%   sizes 3,2,3), and that animal's "Early"/"Middle"/"Late" point is the
%   MEAN of its MDS coordinates across whichever sessions fall in that
%   chunk. This gives exactly ONE point per animal per block (3 points
%   total per animal) regardless of how many total sessions it has --
%   directly comparable across animals with very different session
%   counts, unlike raw session-index alignment.
%
%   The permutation-test logic downstream (exact enumeration of all
%   nchoosek(nAllMice,nGroup1) group-label permutations, contrast
%   statistic, one-sided empirical p-value) is otherwise IDENTICAL to
%   runAmongDistancePermutationTest -- only the preprocessing differs.
%
%   amongOut = runAmongDistancePermutationTest_blocked(Y, sessInfo, groupDefs, ...)
%
% INPUTS
%   Y         : [S x dimFull] MDS coordinates for ONE stream.
%   sessInfo  : matching sessInfo table (needs mouseId, sessWithin).
%   groupDefs : scalar struct, exactly 2 groups (e.g. fast/slow).
%
% NAME-VALUE ARGS
%   'dims'         : which columns of Y to use. Default: 1:size(Y,2).
%   'amongMetric'  : 'euclidean' (default) or 'mahal' -- SAME option as
%                    the unblocked version. If 'mahal', the whitening
%                    covariance is computed ONCE from the full set of
%                    block-averaged points, pooled across all mice and
%                    all 3 blocks -- a shared ruler, not built from
%                    either group, so this stays non-circular.
%   'contrastStat' : 'meanDiff' (default, averaged across all 3 blocks)
%                    or 'finalDiff' (Late block only).
%   'direction'    : 'slow_minus_fast' (default) or 'fast_minus_slow'.
%   'useExact'     : true (default). Same exact-enumeration logic as
%                    before (tractable here regardless, since group
%                    sizes haven't changed -- still C(9,4)=126 for this
%                    project's cohort).
%   'maxExact'     : 5000 (default).
%   'nPerm'        : 1000 (default, only used if exact is intractable).
%   'permSeed'     : [] (default) or scalar seed.
%   'doPlot'       : true (default) -- plots both groups' Early/Middle/
%                    Late trajectories with categorical x-tick labels.
%   'verbose'      : true (default).
%   'doSave', 'saveDir', 'saveTag' : as in the other stats functions here.
%
% OUTPUT (amongOut)
%   .blockLabels           : {"Early","Middle","Late"}
%   .amongY_group1, .amongY_group2 : [3 x 1] OBSERVED trajectories
%   .obsContrast, .nullContrast, .pValue
%   .nPartitions, .isExact, .direction, .contrastStat
%   .blockAssignmentTbl    : per-(mouse,session) block assignment, for
%                            sanity-checking the preprocessing itself
%                            (e.g. confirming no session was dropped and
%                            each animal's blocks are contiguous)

%% -------------------- parse options --------------------
p = inputParser;
p.addParameter('dims', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
p.addParameter('amongMetric', 'euclidean', @(s) any(strcmpi(string(s), ["euclidean","mahal"])));
p.addParameter('contrastStat', 'meanDiff', @(s) any(strcmpi(string(s), ["meanDiff","finalDiff"])));
p.addParameter('direction', 'slow_minus_fast', @(s) any(strcmpi(string(s), ["slow_minus_fast","fast_minus_slow"])));
p.addParameter('useExact', true, @(x) islogical(x) && isscalar(x));
p.addParameter('maxExact', 5000, @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('nPerm', 1000, @(x) isnumeric(x) && isscalar(x) && x>=1 && x==round(x));
p.addParameter('permSeed', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('doPlot', true, @(x) islogical(x) && isscalar(x));
p.addParameter('verbose', true, @(x) islogical(x) && isscalar(x));
p.addParameter('doSave', true, @(x) islogical(x) && isscalar(x));
p.addParameter('saveDir', '', @(s) ischar(s) || isstring(s));
p.addParameter('saveTag', '', @(s) ischar(s) || isstring(s));
p.parse(varargin{:});
opt = p.Results;

amongMetric  = lower(string(opt.amongMetric));
direction    = lower(string(opt.direction));
contrastStat = lower(string(opt.contrastStat));
blockLabels  = ["Early","Middle","Late"];
nBlocks      = 3;

if ~isempty(opt.permSeed)
    rng(opt.permSeed);
end

%% -------------------- validate + build fixed mouse universe --------------------
groupNames = fieldnames(groupDefs);
assert(numel(groupNames) == 2, 'runAmongDistancePermutationTest_blocked requires exactly 2 groups.');

group1IdC = groupDefs.(groupNames{1});
group2IdC = groupDefs.(groupNames{2});
overlap = intersect(group1IdC, group2IdC);
assert(isempty(overlap), 'Animal(s) %s appear in both groups.', strjoin(overlap, ', '));

allMiceIdC = [group1IdC(:)', group2IdC(:)'];
nAllMice   = numel(allMiceIdC);
n1         = numel(group1IdC);

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
d  = numel(dims);

%% -------------------- preprocessing: assign each session to Early/Middle/Late,
%% per animal, then average within each (animal, block) pair --------------------
blockAssign = nan(numel(mouseId), 1);   % 1/2/3, NaN for rows outside allMiceIdC

for im = 1:nAllMice
    rowsIm = find(mouseId == allMiceIdC{im});
    Nim = numel(rowsIm);
    assert(Nim >= nBlocks, ...
        'Animal "%s" has only %d session(s) -- need at least %d for a 3-block split.', ...
        allMiceIdC{im}, Nim, nBlocks);

    [~, ord] = sort(sessWithin(rowsIm), 'ascend');
    rowsSorted = rowsIm(ord);

    edges = round(linspace(0, Nim, nBlocks+1));   % e.g. N=8 -> [0 3 5 8] -> chunk sizes 3,2,3 (as-equal-as-possible)
    for b = 1:nBlocks
        rangeIdx = (edges(b)+1) : edges(b+1);
        blockAssign(rowsSorted(rangeIdx)) = b;
    end
end

% Sanity-check table: one row per session actually used, showing its
% animal, its own local session index, and which block it landed in.
blockAssignmentTbl = table( ...
    mouseId(~isnan(blockAssign)), sessWithin(~isnan(blockAssign)), blockAssign(~isnan(blockAssign)), ...
    'VariableNames', {'mouseId','sessWithin','block'});

% Per-(mouse,block) MEAN point.
blockPts = nan(nAllMice, nBlocks, d);
for im = 1:nAllMice
    for b = 1:nBlocks
        rowsB = (mouseId == allMiceIdC{im}) & (blockAssign == b);
        blockPts(im, b, :) = mean(Yd(rowsB, :), 1);
    end
end

if opt.verbose
    fprintf('\n============================================================\n');
    fprintf('Block preprocessing (Early/Middle/Late), per animal:\n');
    for im = 1:nAllMice
        rowsIm = (mouseId == allMiceIdC{im});
        nPerBlock = arrayfun(@(b) sum(blockAssign(rowsIm) == b), 1:nBlocks);
        fprintf('  %-8s : total %d sessions -> Early=%d, Middle=%d, Late=%d\n', ...
            allMiceIdC{im}, sum(rowsIm), nPerBlock(1), nPerBlock(2), nPerBlock(3));
    end
    fprintf('============================================================\n');
end

%% -------------------- shared, non-circular Mahalanobis ruler (if requested) --------------------
% Computed ONCE from ALL block-averaged points, pooled across every mouse
% and every block -- a shared metric, not built from either group.
W_mahal = [];
if amongMetric == "mahal"
    allBlockPts = reshape(permute(blockPts, [2 1 3]), nAllMice*nBlocks, d);
    SigmaGlobal = cov(allBlockPts, 'omitrows');
    SigmaGlobal = (SigmaGlobal + SigmaGlobal')/2 + 1e-10*eye(d);
    W_mahal = local_cholStable(SigmaGlobal);
end

%% -------------------- core reusable subroutine --------------------
    function amongTraj = local_amongTrajectory(memberIdC)
        memberMask = ismember(allMiceIdC, memberIdC);
        amongTraj = nan(nBlocks, 1);
        for b = 1:nBlocks
            pts = squeeze(blockPts(memberMask, b, :));
            if size(pts, 1) < 2
                continue;
            end
            switch amongMetric
                case "euclidean"
                    dists = pdist(pts, 'euclidean');
                case "mahal"
                    ptsW = pts / W_mahal;
                    dists = pdist(ptsW, 'euclidean');
            end
            amongTraj(b) = mean(dists);
        end
    end

    function c = local_contrast(traj1, traj2)
        switch direction
            case "slow_minus_fast"
                dvec = traj2 - traj1;
            case "fast_minus_slow"
                dvec = traj1 - traj2;
        end
        switch contrastStat
            case "meandiff"
                c = mean(dvec, 'omitnan');
            case "finaldiff"
                c = dvec(end);   % "Late" block
        end
    end

%% -------------------- observed --------------------
amongY_group1_obs = local_amongTrajectory(group1IdC);
amongY_group2_obs = local_amongTrajectory(group2IdC);
obsContrast = local_contrast(amongY_group1_obs, amongY_group2_obs);

fprintf('\nAmong-distance permutation test, BLOCKED preprocessing (metric: %s, contrast: %s, direction: %s)\n', ...
    amongMetric, contrastStat, direction);
fprintf('Observed contrast: %.4f\n', obsContrast);

%% -------------------- null: exact enumeration or Monte Carlo --------------------
nPartitions = nchoosek(nAllMice, n1);
isExact = opt.useExact && (nPartitions <= opt.maxExact);

if isExact
    fprintf('Using EXACT enumeration: all C(%d,%d) = %d partitions.\n', nAllMice, n1, nPartitions);
    combos = nchoosek(1:nAllMice, n1);
    nDraws = size(combos, 1);
else
    fprintf('Using Monte Carlo: %d random partitions.\n', opt.nPerm);
    nDraws = opt.nPerm;
end

nullContrast = nan(nDraws, 1);
tPerm = tic;
for pI = 1:nDraws
    if isExact
        idx1 = combos(pI, :);
    else
        permOrder = randperm(nAllMice);
        idx1 = permOrder(1:n1);
    end
    idx2 = setdiff(1:nAllMice, idx1);

    traj1 = local_amongTrajectory(allMiceIdC(idx1));
    traj2 = local_amongTrajectory(allMiceIdC(idx2));
    nullContrast(pI) = local_contrast(traj1, traj2);
end

pValue = (1 + sum(nullContrast >= obsContrast, 'omitnan')) / (nDraws + 1);

fprintf('Completed in %.1f sec.\n', toc(tPerm));
fprintf('Observed contrast = %.4f | null mean=%.4f SD=%.4f | p = %.4f\n', ...
    obsContrast, mean(nullContrast,'omitnan'), std(nullContrast,'omitnan'), pValue);
fprintf('============================================================\n');

%% -------------------- optional plot --------------------
if opt.doPlot
    figure('Color','w');
    hold on;
    plot(1:nBlocks, amongY_group1_obs, '-o', 'LineWidth', 2, 'DisplayName', groupNames{1});
    plot(1:nBlocks, amongY_group2_obs, '-o', 'LineWidth', 2, 'DisplayName', groupNames{2});
    set(gca, 'XTick', 1:nBlocks, 'XTickLabel', blockLabels);
    xlim([0.5, nBlocks+0.5]);
    ylabel(sprintf('Mean pairwise distance among group (%s)', amongMetric));
    title('Among-group dispersion, Early/Middle/Late blocks (all sessions used)');
    legend('Location','best');
    grid on; box off;
    hold off;
end

%% -------------------- package output --------------------
amongOut = struct();
amongOut.blockLabels        = blockLabels;
amongOut.amongY_group1      = amongY_group1_obs;
amongOut.amongY_group2      = amongY_group2_obs;
amongOut.groupNames         = groupNames;
amongOut.obsContrast        = obsContrast;
amongOut.nullContrast       = nullContrast;
amongOut.pValue             = pValue;
amongOut.nPartitions        = nPartitions;
amongOut.isExact            = isExact;
amongOut.direction          = char(direction);
amongOut.contrastStat       = char(contrastStat);
amongOut.amongMetric        = char(amongMetric);
amongOut.blockAssignmentTbl = blockAssignmentTbl;

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
        warning('runAmongDistancePermutationTest_blocked:noSaveTag', 'No ''saveTag'' provided -- using "unlabeled".');
    end
    dateStr  = char(datetime('today','Format','MMddyy'));
    saveName = sprintf('amongDistancePermTest_blocked_%s_%s.mat', tag, dateStr);
    saveFullPath = fullfile(outDir, saveName);
    save(saveFullPath, 'amongOut', 'groupDefs');
    amongOut.save = struct('didSave', true, 'file', saveFullPath);
    fprintf('\nSaved to:\n%s\n', saveFullPath);
else
    amongOut.save = struct('didSave', false, 'file', '');
end

end

%% ========================= local helper =========================
function R = local_cholStable(Sigma)
Sigma = (Sigma + Sigma')/2;
reg = 0;
for it = 1:10
    try
        R = chol(Sigma + reg*eye(size(Sigma,1)));
        return;
    catch
        if reg == 0
            reg = 1e-12;
        else
            reg = reg * 10;
        end
    end
end
error('Cholesky failed even after regularization. Sigma may be badly conditioned.');
end