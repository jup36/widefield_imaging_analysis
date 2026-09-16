function permOut = runConvergenceEllipsoidPermutationTest(Y, sessInfo, groupDefs, varargin)
%RUNCONVERGENCEELLIPSOIDPERMUTATIONTEST
%   Full circularity-control permutation test for the fast-learner
%   convergence analysis. This is the COMPLETE version of the mID-shuffle
%   null discussed previously -- that earlier version (inside
%   runConvergenceRMANOVA_fastSlow.m) held the reference ellipsoid FIXED
%   (built from the TRUE fast learners) and only shuffled which label got
%   attached to already-computed distances. That's a partial control: it
%   doesn't address the fact that the reference itself was built FROM the
%   fast learners, so they're guaranteed to look close to it by
%   construction to some degree.
%
%   THIS function closes that gap. Under each of nPerm permutations:
%     1) Randomly select numel(groupDefs.fast) mice (out of ALL mice
%        listed across groupDefs) to be "pseudo-fast" learners -- a
%        genuine relabeling of mouse IDs, not the true fast/slow split.
%     2) Rebuild the reference centroid (mu) and covariance (Sigma) from
%        ONLY the pseudo-fast mice's last 'nFinalSess' sessions (same
%        logic as computeFastLearnerEllipsoid, duplicated here
%        self-contained -- see NOTE ON SELF-CONTAINMENT below).
%     3) Recompute EVERY mouse's distance-to-this-pseudo-reference
%        (Mahalanobis by default, or Euclidean via 'refMode') for every
%        session.
%     4) Refit the SAME repeated-measures ANOVA design (final
%        nSessAnalyze sessions per mouse), using the PSEUDO-fast/pseudo-
%        slow labels as the Group factor -- i.e. fully self-referential,
%        matching the true analysis's structure exactly except for which
%        mice get called "fast."
%     5) Extract the same target F-statistic (Group main effect, or
%        Group:Session interaction).
%   The TRUE-label version (using the actual fast learners) is computed
%   ONCE through this identical code path (not imported from a separately
%   -computed result elsewhere), guaranteeing observed and null are
%   apples-to-apples. The empirical p-value is the fraction of pseudo-
%   relabelings whose F-statistic meets or exceeds the true one.
%
%   NOTE ON SELF-CONTAINMENT: this function does NOT call
%   computeFastLearnerEllipsoid or the local mouseId/sessWithin-inference
%   helpers embedded in your six-stream script -- those are LOCAL
%   functions defined inside that script file and are not callable from
%   any other file. The handful of lines of actual math they perform
%   (pool points -> mean/cov -> regularize; Mahalanobis distance) are
%   duplicated directly below instead, to avoid any risk of silently
%   calling the wrong/inaccessible copy (an issue that came up repeatedly
%   elsewhere in this project with duplicated functions across files).
%
%   permOut = runConvergenceEllipsoidPermutationTest(Y, sessInfo, groupDefs, ...)
%
% INPUTS
%   Y        : [S x dimFull] MDS coordinates for ONE stream (e.g. Y_cr),
%              i.e. the SAME Y that was fed into
%              plotMDSConvergenceAndAmongDistance for the real analysis.
%   sessInfo : the matching sessInfo table (must have columns mouseId,
%              sessWithin -- these come pre-populated with clean names
%              from buildSessionFeaturesFromXcorr_dir, so no inference
%              fallback logic is needed here, unlike
%              plotMDSConvergenceAndAmongDistance's more defensive
%              handling of arbitrary sessInfo tables).
%   groupDefs : scalar struct, field name = group label, value = cellstr
%              of animal IDs -- SAME structure used elsewhere, e.g.
%                groupDefs.fast = {'m1044','m1045','m1092','m1094'};
%                groupDefs.slow = {'m1048','m1049','m1613','m1859','m1873'};
%              Only 2 groups supported here (this test is specifically
%              about "is the fast-defining group special").
%
% NAME-VALUE ARGS
%   'dims'         : which columns of Y to use. Default: 1:size(Y,2)
%                    (matches plotMDSConvergenceAndAmongDistance's own
%                    default, and should match whatever dims your real
%                    convergence analysis used -- e.g. [1 2 3 4] if you
%                    passed that explicitly there).
%   'nFinalSess'   : how many of each pseudo-fast mouse's final sessions
%                    define the reference. Default: 3 (matches
%                    computeFastLearnerEllipsoid's usage throughout this
%                    project).
%   'confLevel'    : passed through for consistency/bookkeeping only --
%                    NOT actually used in the F-statistic computation
%                    (that only needs mu/Sigma, not a confidence radius).
%                    Default: 0.95.
%   'refMode'      : 'mahal' (default, matches this specific request) or
%                    'euclidean'. Use 'euclidean' if the stream you're
%                    testing used 'refMode','euclidean' in its real
%                    convergence analysis (Hit and combinedAll did, per
%                    this project's earlier choices) -- for an
%                    apples-to-apples null, this MUST match what the real
%                    analysis used.
%   'nSessAnalyze' : final-session window depth for the RM-ANOVA. Default:
%                    [], meaning auto = minimum session count among all
%                    mice in groupDefs (same convention as
%                    runConvergenceRMANOVA_fastSlow).
%   'permTarget'   : 'group' (default) or 'interaction' -- which
%                    F-statistic to build a null distribution for. See
%                    runConvergenceRMANOVA_fastSlow's docstring for the
%                    full explanation of this choice; same meaning here.
%   'nPerm'        : 1000 (default).
%   'permSeed'     : [] (default, non-reproducible) or a scalar seed.
%   'verbose'      : true (default).
%   'doSave'       : true (default). Requires 'saveDir' if true.
%   'saveDir'      : folder to save into.
%   'saveTag'      : short label for the saved filename (e.g. 'cr').
%
% OUTPUT (permOut)
%   .obsF          : observed F-statistic, TRUE fast/slow labels, using
%                    THIS function's own (self-contained) computation
%                    path -- should closely match (allowing for the
%                    refMode/dims settings actually used) whatever
%                    runConvergenceRMANOVA_fastSlow reported for the
%                    equivalent target, as a sanity cross-check.
%   .nullF         : [nPerm x 1] null distribution.
%   .pEmpirical    : (1 + sum(nullF >= obsF)) / (nPerm + 1)
%   .permTarget, .nSessAnalyze, .nFast, .allMiceIdC, .dims, .refMode
%   .pseudoFastSets : {nPerm x 1} cellstr, which mice were pseudo-fast on
%                    each permutation -- kept for diagnostics (e.g.
%                    checking the null isn't dominated by a handful of
%                    mice that happen to look "convergent" under almost
%                    any grouping).
%
% See also: runConvergenceRMANOVA_fastSlow, computeFastLearnerEllipsoid

%% -------------------- parse options --------------------
p = inputParser;
p.addParameter('dims', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
p.addParameter('nFinalSess', 3, @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('confLevel', 0.95, @(x) isnumeric(x) && isscalar(x) && x>0 && x<1); %#ok<INUSA> % bookkeeping only, see docstring
p.addParameter('refMode', 'mahal', @(s) any(strcmpi(string(s), ["mahal","euclidean"])));
p.addParameter('nSessAnalyze', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>=2 && x==round(x)));
p.addParameter('permTarget', 'group', @(s) any(strcmpi(string(s), ["group","interaction"])));
p.addParameter('nPerm', 1000, @(x) isnumeric(x) && isscalar(x) && x>=1 && x==round(x));
p.addParameter('permSeed', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('verbose', true, @(x) islogical(x) && isscalar(x));
p.addParameter('doSave', true, @(x) islogical(x) && isscalar(x));
p.addParameter('saveDir', '', @(s) ischar(s) || isstring(s));
p.addParameter('saveTag', '', @(s) ischar(s) || isstring(s));
p.parse(varargin{:});
opt = p.Results;

refMode    = lower(string(opt.refMode));
permTarget = lower(string(opt.permTarget));

if ~isempty(opt.permSeed)
    rng(opt.permSeed);
end

%% -------------------- validate + build fixed mouse universe --------------------
groupNames = fieldnames(groupDefs);
assert(numel(groupNames) == 2, ...
    'runConvergenceEllipsoidPermutationTest is designed for exactly 2 groups (e.g. fast/slow); got %d.', numel(groupNames));

trueFastIdC = groupDefs.(groupNames{1});
trueSlowIdC = groupDefs.(groupNames{2});
overlap = intersect(trueFastIdC, trueSlowIdC);
assert(isempty(overlap), 'Animal(s) %s appear in both groups.', strjoin(overlap, ', '));

allMiceIdC = [trueFastIdC(:)', trueSlowIdC(:)'];
nAllMice   = numel(allMiceIdC);
nFast      = numel(trueFastIdC);

assert(all(ismember('mouseId', sessInfo.Properties.VariableNames)) && ...
       all(ismember('sessWithin', sessInfo.Properties.VariableNames)), ...
    'sessInfo must have mouseId and sessWithin columns (as produced by buildSessionFeaturesFromXcorr_dir).');

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

%% -------------------- fixed per-mouse session bookkeeping (independent of grouping) --------------------
perMouseMax = zeros(nAllMice, 1);
for im = 1:nAllMice
    perMouseMax(im) = max(sessWithin(mouseId == allMiceIdC{im}));
end
minDepth = min(perMouseMax);   % max sessWithin == depth, since sessWithin starts at 1

if isempty(opt.nSessAnalyze)
    nSessAnalyze = minDepth;
else
    nSessAnalyze = opt.nSessAnalyze;
    if nSessAnalyze > minDepth
        error('Requested nSessAnalyze=%d exceeds the minimum available depth (%d).', nSessAnalyze, minDepth);
    end
end

sessFromEnd = nan(numel(mouseId), 1);
for im = 1:nAllMice
    rowsIm = (mouseId == allMiceIdC{im});
    sessFromEnd(rowsIm) = perMouseMax(im) - sessWithin(rowsIm);
end

withinDesign = table((1:nSessAnalyze)', 'VariableNames', {'Session'});
withinDesign.Session = categorical(withinDesign.Session);
varNames = arrayfun(@(k) sprintf('Session%d', k), 1:nSessAnalyze, 'UniformOutput', false);
rmFormula = sprintf('%s-%s ~ Group', varNames{1}, varNames{end});

fprintf('\n============================================================\n');
fprintf('Full circularity-control permutation test\n');
fprintf('refMode=%s, dims=%s, nFinalSess=%d, nSessAnalyze=%d, nFast=%d/%d mice, target=%s, nPerm=%d\n', ...
    refMode, mat2str(dims), opt.nFinalSess, nSessAnalyze, nFast, nAllMice, permTarget, opt.nPerm);
fprintf('============================================================\n');

%% -------------------- core reusable subroutine (self-contained) --------------------
    function Fstat = local_computeFstat(fastSetIdC)
        % 1) Build the pseudo (or true) reference from ONLY fastSetIdC's
        % last nFinalSess sessions, pooled directly (matching
        % computeFastLearnerEllipsoid's logic exactly, duplicated here).
        idxPts = [];
        for k = 1:numel(fastSetIdC)
            iMouse = find(mouseId == fastSetIdC{k});
            [~, ord] = sort(sessWithin(iMouse), 'ascend');
            iMouse = iMouse(ord);
            take = max(1, numel(iMouse)-opt.nFinalSess+1) : numel(iMouse);
            idxPts = [idxPts; iMouse(take)]; %#ok<AGROW>
        end
        idxPts = unique(idxPts, 'stable');

        Xref  = Yd(idxPts, :);
        muRef = mean(Xref, 1);

        SigmaRef = cov(Xref);
        SigmaRef = (SigmaRef + SigmaRef')/2 + 1e-8*eye(d);

        % 2) Distance from EVERY session (all mice) to this reference.
        switch refMode
            case "euclidean"
                diffAll = Yd - muRef;
                distAll = sqrt(sum(diffAll.^2, 2));
            case "mahal"
                R = local_cholStable(SigmaRef);
                Z = (Yd - muRef) / R;
                distAll = sqrt(sum(Z.^2, 2));
        end

        % 3) Group label per row, SELF-REFERENTIALLY tied to fastSetIdC.
        groupRow = repmat("slow", numel(mouseId), 1);
        groupRow(ismember(mouseId, string(fastSetIdC))) = "fast";

        % 4) Pivot to wide format, one row per mouse in allMiceIdC.
        wideVals = nan(nAllMice, nSessAnalyze);
        groupOf  = strings(nAllMice, 1);
        for im = 1:nAllMice
            rowsIm = (mouseId == allMiceIdC{im}) & (sessFromEnd < nSessAnalyze);
            sessPosIm = nSessAnalyze - sessFromEnd(rowsIm);
            [~, ord] = sort(sessPosIm);
            distIm = distAll(rowsIm);
            wideVals(im, :) = distIm(ord)';
            groupOf(im) = groupRow(find(rowsIm, 1, 'first'));
        end

        wideTablePerm = array2table(wideVals, 'VariableNames', varNames);
        wideTablePerm.Group = categorical(groupOf);
        wideTablePerm = wideTablePerm(:, ['Group', varNames]);

        rmPerm = fitrm(wideTablePerm, rmFormula, 'WithinDesign', withinDesign);

        switch permTarget
            case "group"
                ranovaBetweenPerm = ranova(rmPerm, 'WithinModel', '1');
                Fstat = ranovaBetweenPerm{'Group', 'F'};
            case "interaction"
                ranovaPerm = ranova(rmPerm);
                Fstat = ranovaPerm{'Group:Session', 'F'};
        end
    end

%% -------------------- observed (TRUE labels), through the SAME code path --------------------
obsF = local_computeFstat(trueFastIdC);
fprintf('Observed F (true fast/slow labels): %.4f\n', obsF);

%% -------------------- permutation null --------------------
nullF = nan(opt.nPerm, 1);
pseudoFastSets = cell(opt.nPerm, 1);

tPerm = tic;
for pI = 1:opt.nPerm
    permOrder = randperm(nAllMice);
    pseudoFastIdC = allMiceIdC(permOrder(1:nFast));
    pseudoFastSets{pI} = pseudoFastIdC;

    nullF(pI) = local_computeFstat(pseudoFastIdC);

    if opt.verbose && (mod(pI, 100) == 0 || pI == opt.nPerm)
        fprintf('  permutation %d/%d (%.1fs elapsed)\n', pI, opt.nPerm, toc(tPerm));
    end
end

pEmpirical = (1 + sum(nullF >= obsF)) / (opt.nPerm + 1);

fprintf('\nPermutation test complete in %.1f sec.\n', toc(tPerm));
fprintf('Observed F = %.4f | null F: mean=%.4f, SD=%.4f, 95th pct=%.4f | empirical p = %.4f (n=%d perms)\n', ...
    obsF, mean(nullF,'omitnan'), std(nullF,'omitnan'), prctile(nullF,95), pEmpirical, opt.nPerm);

%% -------------------- package output --------------------
permOut = struct();
permOut.obsF           = obsF;
permOut.nullF          = nullF;
permOut.pEmpirical     = pEmpirical;
permOut.permTarget     = char(permTarget);
permOut.nSessAnalyze   = nSessAnalyze;
permOut.nFast          = nFast;
permOut.allMiceIdC     = allMiceIdC;
permOut.trueFastIdC    = trueFastIdC;
permOut.dims           = dims;
permOut.refMode        = char(refMode);
permOut.nFinalSess     = opt.nFinalSess;
permOut.pseudoFastSets = pseudoFastSets;

%% -------------------- optional save to disk --------------------
if opt.doSave
    if strlength(strtrim(string(opt.saveDir))) == 0
        error(['doSave is true (the default) but no ''saveDir'' was provided. ' ...
               'Pass ''saveDir'',''<path>'' (and ideally ''saveTag'',''<label>''), or ''doSave'',false.']);
    end
    outDir = char(string(opt.saveDir));
    if exist(outDir, 'dir') ~= 7
        mkdir(outDir);
    end
    tag = strtrim(char(string(opt.saveTag)));
    if isempty(tag)
        tag = 'unlabeled';
        warning('runConvergenceEllipsoidPermutationTest:noSaveTag', ...
            'No ''saveTag'' provided -- saving with tag "unlabeled".');
    end
    dateStr  = char(datetime('today','Format','MMddyy'));
    saveName = sprintf('convergenceEllipsoidPermTest_%s_%s_n%d_%s.mat', tag, permTarget, nSessAnalyze, dateStr);
    saveFullPath = fullfile(outDir, saveName);
    save(saveFullPath, 'permOut', 'groupDefs');
    permOut.save = struct('didSave', true, 'file', saveFullPath);
    fprintf('\nSaved permutation test output to:\n%s\n', saveFullPath);
else
    permOut.save = struct('didSave', false, 'file', '');
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