%% ========================================================================
%  SYNOPSIS
%  ========================================================================
%  For ONE OR MORE target motifs, plots the mean positive-lag residual
%  xcorr in a chosen DIRECTION relative to those motifs, split into two
%  subplots by whether the other motif is Go-preferring or NoGo-
%  preferring, one point per animal (session-averaged), grouped fast vs.
%  slow learners.
%
%  targetMotif CAN BE AN ARRAY (e.g. targetMotif = goIdx, to ask "how does
%  the WHOLE Go-preferring ensemble drive/get driven by other motifs,"
%  rather than just one specific motif). When targetMotif has more than
%  one element, the SAME per-target extraction below is computed
%  INDEPENDENTLY for each target motif (each with its own self-exclusion
%  from its own preference set), and those per-target scalars are then
%  averaged across all targets, PER ANIMAL -- i.e. a two-stage average
%  (per-target, then across-target), matching the per-session-then-
%  per-animal averaging convention used everywhere else in this project.
%  This is NOT the same as pooling every target's raw xcorr entries into
%  one giant vector before averaging -- each target motif contributes
%  exactly one number to the across-target average, regardless of how
%  many Go-/NoGo-preferring partners it individually has.
%
%  DIRECTION CONVENTION (established/verified elsewhere in this project):
%    xcorrPosLagMat(row, col): COLUMN motif leads, ROW motif follows.
%
%    direction = 'into' : target motif(s) are the FOLLOWER. Fix ROW =
%                targetMotif, vary COLUMN across other motifs (the
%                leaders). Answers "what drives INTO the target motif(s)?"
%                -- this is A(targetMotif, :), a ROW trace per target.
%    direction = 'from' : target motif(s) are the LEADER. Fix COLUMN =
%                targetMotif, vary ROW across other motifs (the
%                followers). Answers "what does the target motif(s) drive?"
%                -- this is A(:, targetMotif), a COLUMN trace per target.
%
%  Reuses averageMatAcrossSessions and plotJitteredGroups verbatim from
%  the prior significance-testing script -- no new averaging/plotting
%  mechanics introduced here.
%  ========================================================================

%% -------------------- whereabouts --------------------
filePath_xcorrRez = compatiblepath("Z:\Rodent Data\dualImaging_parkj\collectData\xcorr_withinTrialTimeShuffle\xcorrResidualPosLag_motifH_collect_timeshuffle_n1000_080526.mat");
load(filePath_xcorrRez, 'headerC', 'mIdC', 'xcorrRezC');

% -------- Go-/NoGo-preferring motif labels (bin-wise-difference definition) --------
goNogoDir = fullfile(compatiblepath('/Volumes/buschman/Rodent Data/dualImaging_parkj'), ...
    'collectData', 'motifGLM_goNogoPreference');
goNogoFiles = dir(fullfile(goNogoDir, 'motifGoNogoPreference_fromGLM_beta_permNull_strictAndDiffBins_*.mat'));
assert(~isempty(goNogoFiles), 'No Go/NoGo preference result found under %s.', goNogoDir);
[~, mostRecentIdx_goNogo] = max([goNogoFiles.datenum]);
goNogoPath = fullfile(goNogoDir, goNogoFiles(mostRecentIdx_goNogo).name);
fprintf('Using Go/NoGo preference result:\n%s\n', goNogoPath);
load(goNogoPath, 'resultT_strict');

goIdx   = find(resultT_strict.isGoPreferring)';
nogoIdx = find(resultT_strict.isNoGoPreferring)';
% load(goNogoPath, 'resultT_diffBins');
%  
% goIdx   = find(resultT_diffBins.isGoPreferring_diff)';
% nogoIdx = find(resultT_diffBins.isNoGoPreferring_diff)';

assert(~any(resultT_strict.isGoPreferring & resultT_strict.isNoGoPreferring), ...
    'A motif is labeled both Go- and NoGo-preferring simultaneously -- check labeling logic before proceeding.');

fast_learners = {'m1044', 'm1045', 'm1092', 'm1094'};
slow_learners = {'m1048', 'm1049', 'm1613', 'm1859', 'm1873'};

%% -------------------- user parameters --------------------
targetMotif = goIdx;                          % SCALAR or ARRAY of motifs of interest, e.g. targetMotif = goIdx
direction   = 'into';                         % 'into' (target = follower, row trace) or 'from' (target = leader, column trace)
fieldName   = 'xcorrPosLagMat_cr_residual';   % 'xcorrPosLagMat_hit_residual' for Hit trials instead
trialTag    = 'CR trials (residual)';         % just for titles/labels -- keep in sync with fieldName above

fastColor = [0.85 0.33 0.10];   % orange, matches this project's convention for fast learners elsewhere
slowColor = [0.00 0.45 0.74];   % blue,   matches this project's convention for slow learners elsewhere

figSaveDir = "Z:\Rodent Data\dualImaging_parkj\collectFigure\motifTDR\xcorr_mds\motif_posLag_xcorrMatrix";
if exist(figSaveDir, 'dir') ~= 7
    mkdir(figSaveDir);
end

assert(any(strcmpi(direction, {'into','from'})), 'direction must be ''into'' or ''from''.');

%% -------------------- validate EVERY target motif, up front --------------------
% Each target motif needs its own preference-class check and its own
% self-exclusion set (computed once here, reused inside the animal loop
% below -- no need to redo the ismember/assert checks per animal, since
% these only depend on the motif identity, not the animal).
targetMotif = unique(targetMotif(:)', 'stable');
nTargets = numel(targetMotif);

isTargetGo_perTarget   = false(1, nTargets);
targetLabel_perTarget  = strings(1, nTargets);

for ti = 1:nTargets
    t = targetMotif(ti);
    isTGo   = ismember(t, goIdx);
    isTNoGo = ismember(t, nogoIdx);
    assert(isTGo || isTNoGo, ...
        'targetMotif element %d is neither Go- nor NoGo-preferring per resultT_strict -- check targetMotif or the preference labeling.', t);
    assert(~(isTGo && isTNoGo), ...
        'targetMotif element %d is labeled BOTH Go- and NoGo-preferring -- check labeling logic.', t);
    isTargetGo_perTarget(ti) = isTGo;
    targetLabel_perTarget(ti) = ternary_local(isTGo, "Go", "NoGo");
end

fprintf('Target motif(s) (n=%d): %s\n', nTargets, mat2str(targetMotif));
fprintf('  Preference per target: %s\n', strjoin(targetLabel_perTarget, ', '));
fprintf('Direction: %s\n', direction);

%% -------------------- per-animal, session-averaged extraction --------------------
allAnimals = [fast_learners, slow_learners];
mIdCol = mIdC(:,1);

goVal_perAnimal   = nan(numel(allAnimals), 1);
nogoVal_perAnimal = nan(numel(allAnimals), 1);
animalHasData     = false(numel(allAnimals), 1);

for a = 1:numel(allAnimals)
    thisID = allAnimals{a};
    rowIdx = find(strcmp(mIdCol, thisID), 1, 'first');
    if isempty(rowIdx)
        warning('Animal "%s" not found in mIdC -- skipping.', thisID);
        continue;
    end

    Manimal = averageMatAcrossSessions(xcorrRezC, rowIdx, @(entry) entry.obs.(fieldName));
    if isempty(Manimal)
        warning('Animal "%s" has no valid sessions for field "%s" -- skipping.', thisID, fieldName);
        continue;
    end

    % ---- inner loop: compute the SAME per-target extraction independently
    % for each target motif (own self-exclusion per target), then average
    % across targets -- two-stage average, matching this project's
    % established per-session-then-per-animal averaging convention. ----
    goVal_perTarget   = nan(1, nTargets);
    nogoVal_perTarget = nan(1, nTargets);

    for ti = 1:nTargets
        t = targetMotif(ti);

        % Self-exclusion is PER TARGET (a different target motif may or
        % may not itself be Go- or NoGo-preferring, so the set it gets
        % excluded from differs target to target).
        goIdx_noSelf_t   = goIdx(goIdx ~= t);
        nogoIdx_noSelf_t = nogoIdx(nogoIdx ~= t);

        switch lower(direction)
            case 'into'
                % target = follower: ROW = t, vary COLUMN
                traceVec = Manimal(t, :);
            case 'from'
                % target = leader: COLUMN = t, vary ROW
                traceVec = Manimal(:, t)';
        end

        goVal_perTarget(ti)   = mean(traceVec(goIdx_noSelf_t),   'omitnan');
        nogoVal_perTarget(ti) = mean(traceVec(nogoIdx_noSelf_t), 'omitnan');
    end

    goVal_perAnimal(a)   = mean(goVal_perTarget,   'omitnan');
    nogoVal_perAnimal(a) = mean(nogoVal_perTarget, 'omitnan');
    animalHasData(a) = true;
end

if any(~animalHasData)
    fprintf('WARNING: %d/%d animals had no valid data and were excluded from the plot: %s\n', ...
        sum(~animalHasData), numel(allAnimals), strjoin(allAnimals(~animalHasData), ', '));
end

%% -------------------- split into fast/slow for plotting --------------------
nFast = numel(fast_learners);
isFast = false(numel(allAnimals), 1);
isFast(1:nFast) = true;   % allAnimals = [fast_learners, slow_learners], by construction above

goFastVals   = goVal_perAnimal(isFast & animalHasData);
goSlowVals   = goVal_perAnimal(~isFast & animalHasData);
nogoFastVals = nogoVal_perAnimal(isFast & animalHasData);
nogoSlowVals = nogoVal_perAnimal(~isFast & animalHasData);

%% -------------------- plot: two subplots, target->Go and target->NoGo --------------------
if nTargets == 1
    targetDisp = sprintf('motif %d', targetMotif);
else
    targetDisp = sprintf('%d motifs', nTargets);
end

switch lower(direction)
    case 'into'
        arrowGo   = sprintf('Go-preferring \\rightarrow %s', targetDisp);
        arrowNoGo = sprintf('NoGo-preferring \\rightarrow %s', targetDisp);
        yLabelStr = sprintf('Mean pos-lag xcorr INTO %s', targetDisp);
    case 'from'
        arrowGo   = sprintf('%s \\rightarrow Go-preferring', targetDisp);
        arrowNoGo = sprintf('%s \\rightarrow NoGo-preferring', targetDisp);
        yLabelStr = sprintf('Mean pos-lag xcorr FROM %s', targetDisp);
end

fig = figure('Color','w', 'Position', [200 200 800 420]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing','compact', 'Padding','compact');

axGo = nexttile(tl, 1);
plotJitteredGroups(axGo, goFastVals, goSlowVals, fastColor, slowColor);
ylabel(axGo, yLabelStr);
title(axGo, arrowGo);

axNoGo = nexttile(tl, 2);
plotJitteredGroups(axNoGo, nogoFastVals, nogoSlowVals, fastColor, slowColor);
ylabel(axNoGo, yLabelStr);
title(axNoGo, arrowNoGo);

% Match y-limits across both panels for direct visual comparability.
yl1 = ylim(axGo); yl2 = ylim(axNoGo);
yl = [min(yl1(1), yl2(1)), max(yl1(2), yl2(2))];
ylim(axGo, yl); ylim(axNoGo, yl);

sgtitle(tl, sprintf('Session-averaged, per-animal (%s) %s -- %s', direction, targetDisp, trialTag), ...
    'Interpreter', 'tex');

%% -------------------- save --------------------
dateStr = string(datetime('today','Format','MMddyy'));
if nTargets == 1
    targetTag = sprintf('motif%d', targetMotif);
else
    targetTag = sprintf('nTargets%d', nTargets);
end
outFile = fullfile(figSaveDir, sprintf('%s_%s_goVsNogo_%s_%s.pdf', targetTag, direction, ...
    matlab.lang.makeValidName(trialTag), dateStr));
set(fig, 'InvertHardcopy', 'off');
print(fig, outFile, '-dpdf', '-painters', '-bestfit');
fprintf('\nSaved figure to:\n%s\n', outFile);

%% ========================================================================
%  LOCAL FUNCTIONS (verbatim from the significance-testing script, plus
%  one small helper for the target-preference print statement above)
%  ========================================================================
function out = ternary_local(cond, valTrue, valFalse)
if cond
    out = valTrue;
else
    out = valFalse;
end
end


function plotJitteredGroups(ax, fastVals, slowVals, fastColor, slowColor)
% Jittered strip plot: fast at x=1, slow at x=2, small random horizontal
% jitter for visibility, group mean +/- SEM overlaid as a black errorbar.
jitterWidth = 0.15;

xFast = 1 + (rand(numel(fastVals), 1) - 0.5) * jitterWidth * 2;
xSlow = 2 + (rand(numel(slowVals), 1) - 0.5) * jitterWidth * 2;

hold(ax, 'on');
scatter(ax, xFast, fastVals, 40, fastColor, 'filled', 'MarkerFaceAlpha', 0.7);
scatter(ax, xSlow, slowVals, 40, slowColor, 'filled', 'MarkerFaceAlpha', 0.7);

meanFast = mean(fastVals, 'omitnan'); semFast = std(fastVals, 'omitnan') / sqrt(sum(isfinite(fastVals)));
meanSlow = mean(slowVals, 'omitnan'); semSlow = std(slowVals, 'omitnan') / sqrt(sum(isfinite(slowVals)));

errorbar(ax, 1, meanFast, semFast, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6, 'LineWidth', 1.3, 'CapSize', 8);
errorbar(ax, 2, meanSlow, semSlow, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6, 'LineWidth', 1.3, 'CapSize', 8);

set(ax, 'XTick', [1 2], 'XTickLabel', {'Fast', 'Slow'});
xlim(ax, [0.5 2.5]);
grid(ax, 'on'); box(ax, 'on');
hold(ax, 'off');
end


function Manimal = averageMatAcrossSessions(xcorrRezC, rowIdx, extractFn)
% Generalized version: extractFn is a function handle taking one
% session's result struct and returning a [K x K] matrix.
stack = [];
nSessions = size(xcorrRezC, 2);
for j = 1:nSessions
    entry = xcorrRezC{rowIdx, j};
    if isempty(entry) || ~isstruct(entry)
        continue;
    end
    try
        M = extractFn(entry);
    catch
        continue;   % missing field -- skip
    end
    if isempty(M)
        continue;
    end
    stack = cat(3, stack, M);
end

if isempty(stack)
    Manimal = [];
else
    Manimal = mean(stack, 3, 'omitnan');
end
end