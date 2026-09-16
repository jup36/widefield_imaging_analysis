function statsOut = runConvergenceRMANOVA_fastSlow(convergenceTable, groupDefs, varargin)
%RUNCONVERGENCERMANOVA_FASTSLOW
%   Repeated-measures ANOVA on MDS-convergence distance over each animal's
%   FINAL N sessions (N = the minimum session count among the mice being
%   analyzed, so every subject contributes a complete row -- this is what
%   makes a classic RM-ANOVA valid despite animals having unequal total
%   session counts). Tests the main effect of Group (e.g. fast vs. slow
%   learners), the main effect of Session (position within the final-N
%   window), and the Group x Session interaction, plus per-session
%   post-hoc Group comparisons.
%
%   statsOut = runConvergenceRMANOVA_fastSlow(convergenceTable, groupDefs, ...)
%
% INPUTS
%   convergenceTable : table with (at minimum) columns mouseId,
%                      sessFromEnd, and the response variable (default
%                      'distToRef') -- exactly what
%                      plotMDSConvergenceAndAmongDistance.m now returns as
%                      h.convergenceTable. sessFromEnd = 0 must mean "that
%                      mouse's own last session" (see that function).
%   groupDefs        : scalar struct, field name = group label, value =
%                      cellstr of animal IDs, e.g.
%                        groupDefs.fast = {'m1044','m1045','m1092','m1094'};
%                        groupDefs.slow = {'m1048','m1049','m1613','m1859','m1873'};
%                      Designed for exactly 2 groups (fast/slow); more
%                      groups will still run (multcompare handles >2
%                      levels) but the framing below assumes 2.
%
% NAME-VALUE ARGS
%   'nSessAnalyze'  : how many final sessions to include (default: [],
%                     meaning auto = the minimum session count among all
%                     mice listed in groupDefs -- i.e. "however many
%                     sessions the shortest-recorded mouse has"). If you
%                     supply a value larger than that minimum, this
%                     errors with a clear message (some mouse would have
%                     missing data at that window depth) rather than
%                     silently truncating or imputing.
%   'valueVar'      : which column of convergenceTable to analyze.
%                     Default: 'distToRef'.
%   'comparisonType': passed to multcompare's 'ComparisonType'. Default:
%                     'tukey-kramer' (appropriate for unequal group
%                     sizes, e.g. 4 fast vs. 5 slow learners here).
%   'verbose'       : true (default). Prints the ranova table, Mauchly's
%                     sphericity test, and the per-session post-hoc table.
%   'doSave'        : true (default). Saves statsOut (+ groupDefs) to a
%                     .mat file. Set false to skip saving entirely.
%   'saveDir'       : folder to save into. REQUIRED if doSave is true (no
%                     hardcoded default -- this function has no way to
%                     know your machine's path conventions). Created if
%                     it doesn't already exist.
%   'saveTag'       : a short label included in the saved filename, e.g.
%                     'cr', 'hit', 'combinedCorrect' -- meant to mirror
%                     whatever you're naming the output variable in your
%                     own script (e.g. statsOut_cr -> 'saveTag','cr'),
%                     since MATLAB cannot read that variable name back
%                     out of the caller automatically (see note above
%                     the save block, below). Defaults to 'unlabeled'
%                     with a warning if omitted -- provide one for a
%                     filename you can actually tell apart from others.
%
%   Saved filename: convergenceRMANOVA_<saveTag>_<valueVar>_n<nSessAnalyze>_<MMddyy>.mat
%   e.g. convergenceRMANOVA_cr_distToRef_n8_080526.mat
%
% OUTPUT (statsOut)
%   .nSessAnalyze      : the window depth actually used
%   .mouseList         : animal IDs included, in the order used to build
%                        the wide table (grouped: all of group 1, then
%                        group 2, ... matching fieldnames(groupDefs) order)
%   .groupOf           : cellstr, same order as mouseList, group label
%                        per animal
%   .wideTable         : [nMice x nSessAnalyze] table actually fed to
%                        fitrm (columns Session1..SessionN), plus Group
%   .rm                : the fitted repeatedMeasuresModel object (from
%                        fitrm) -- keep this around if you want to run
%                        additional contrasts/comparisons yourself later
%   .ranovaTbl         : output of ranova(rm) -- Session main effect and
%                        Group:Session interaction ONLY (both effects
%                        that involve the within-subject factor). This
%                        does NOT include the pure Group main effect --
%                        see .ranovaBetweenTbl for that.
%   .ranovaBetweenTbl  : output of ranova(rm,'WithinModel','1') -- the
%                        pure BETWEEN-subjects Group main effect, i.e.
%                        is Group different on average across all
%                        nSessAnalyze sessions, ignoring session
%                        structure entirely. THIS is the classical "main
%                        effect of Group" row (labeled 'Group' in the
%                        table). Not automatically part of ranovaTbl --
%                        requires this separate call.
%   .mauchlyTbl        : output of mauchly(rm) -- sphericity test, applies
%                        to ranovaTbl's within-subject effects only (not
%                        meaningful for ranovaBetweenTbl, which has no
%                        within-subject structure). If its p-value < 0.05,
%                        sphericity is violated and you should report
%                        ranovaTbl's pValueGG (Greenhouse-Geisser
%                        corrected) column instead of the uncorrected
%                        pValue column for the Session and Group:Session
%                        rows.
%   .multCompareTbl    : output of multcompare(rm,'Group','By','Session')
%                        -- pairwise Group comparisons AT EACH Session
%                        level, i.e. exactly "statistical comparison
%                        between groups per session."
%
% IMPORTANT CAVEATS (read before reporting these results)
%   1) DATA LOSS: any session beyond the shortest mouse's final-N window
%      is simply excluded from this analysis. If (say) m1873 has 17
%      sessions and the group minimum is 8, only m1873's LAST 8 sessions
%      are used here -- its first 9 sessions contribute nothing to this
%      test (though they still appear in the plotted convergence figure).
%      This is the necessary trade-off for a complete-design classic
%      RM-ANOVA; a linear mixed-effects model (fitlme, with Session as a
%      continuous or categorical fixed effect and Animal as a random
%      effect) would use ALL the data without this truncation, at the
%      cost of a somewhat less standard/familiar test. Worth considering
%      as a complementary analysis if reviewers ask about the truncation.
%   2) SPHERICITY: classic RM-ANOVA assumes sphericity of the within-
%      subject covariance (Mauchly's test, returned in .mauchlyTbl). With
%      only 4-5 subjects per group this test itself has low power to
%      detect violations -- consider reporting the Greenhouse-Geisser
%      corrected p-value (ranovaTbl.pValueGG) regardless, as a more
%      conservative default.
%   3) UNBALANCED GROUPS: fitrm/ranova handle unequal group sizes (e.g.
%      4 fast vs. 5 slow) without difficulty, but statistical power for
%      the Group effect is inherently limited by the smaller group.
%
% EXAMPLE
%   groupDefs.fast = {'m1044','m1045','m1092','m1094'};
%   groupDefs.slow = {'m1048','m1049','m1613','m1859','m1873'};
%   statsOut = runConvergenceRMANOVA_fastSlow(h_conv_allM_hit.convergenceTable, groupDefs);
%
% See also: plotMDSConvergenceAndAmongDistance, fitrm, ranova, mauchly, multcompare

%% -------------------- parse options --------------------
p = inputParser;
p.addParameter('nSessAnalyze', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>=2 && x==round(x)));
p.addParameter('valueVar', 'distToRef', @(s) ischar(s) || isstring(s));
p.addParameter('comparisonType', 'tukey-kramer', @(s) ischar(s) || isstring(s));
p.addParameter('verbose', true, @(x) islogical(x) && isscalar(x));

% -------- NEW: optional save-to-disk --------
p.addParameter('doSave', true, @(x) islogical(x) && isscalar(x));
p.addParameter('saveDir', '', @(s) ischar(s) || isstring(s));
p.addParameter('saveTag', '', @(s) ischar(s) || isstring(s));
p.parse(varargin{:});
opt = p.Results;

valueVar = char(opt.valueVar);

%% -------------------- validate inputs --------------------
requiredCols = {'mouseId', 'sessFromEnd', valueVar};
missingCols = requiredCols(~ismember(requiredCols, convergenceTable.Properties.VariableNames));
if ~isempty(missingCols)
    error(['convergenceTable is missing required column(s): %s. ' ...
           'Make sure you are passing h.convergenceTable from the UPDATED ' ...
           'plotMDSConvergenceAndAmongDistance.m (older versions did not ' ...
           'compute sessFromEnd).'], strjoin(missingCols, ', '));
end

groupNames = fieldnames(groupDefs);
assert(numel(groupNames) >= 2, 'groupDefs must have at least 2 group fields.');

for gi = 1:numel(groupNames)
    for gj = gi+1:numel(groupNames)
        ov = intersect(groupDefs.(groupNames{gi}), groupDefs.(groupNames{gj}));
        assert(isempty(ov), ...
            'Animal(s) %s appear in both group "%s" and group "%s".', ...
            strjoin(ov, ', '), groupNames{gi}, groupNames{gj});
    end
end

%% -------------------- build mouseList / groupOf, in group order --------------------
mouseList = {};
groupOf   = {};
for gi = 1:numel(groupNames)
    members = groupDefs.(groupNames{gi});
    mouseList = [mouseList, members(:)']; %#ok<AGROW>
    groupOf   = [groupOf,   repmat(groupNames(gi), 1, numel(members))]; %#ok<AGROW>
end

convMouseId = string(convergenceTable.mouseId);
missingMice = mouseList(~ismember(string(mouseList), unique(convMouseId)));
if ~isempty(missingMice)
    error('The following animal(s) from groupDefs were not found in convergenceTable.mouseId: %s', ...
        strjoin(missingMice, ', '));
end

nMice = numel(mouseList);

%% -------------------- determine each mouse's available depth --------------------
% depth(mouse) = number of sessions available = max(sessFromEnd)+1 for
% that mouse (sessFromEnd=0 is the last session, so the count of distinct
% sessFromEnd values -- 0,1,2,...,max -- equals max+1, assuming no gaps,
% which holds by construction since sessFromEnd is derived from sessWithin,
% itself a gap-free valid-session counter).
depthPerMouse = nan(nMice, 1);
for i = 1:nMice
    rows_i = convMouseId == string(mouseList{i});
    depthPerMouse(i) = max(convergenceTable.sessFromEnd(rows_i)) + 1;
end

minDepth = min(depthPerMouse);

if isempty(opt.nSessAnalyze)
    nSessAnalyze = minDepth;
else
    nSessAnalyze = opt.nSessAnalyze;
    if nSessAnalyze > minDepth
        error(['Requested nSessAnalyze=%d exceeds the minimum available depth (%d sessions, ' ...
               'from animal "%s"). Some animal(s) would have missing data at that window size. ' ...
               'Use nSessAnalyze <= %d, or omit it to auto-use the minimum.'], ...
               nSessAnalyze, minDepth, mouseList{depthPerMouse == minDepth}, minDepth);
    end
end

fprintf('\n============================================================\n');
fprintf('RM-ANOVA on convergence distance (%s), final %d session(s) per animal\n', valueVar, nSessAnalyze);
fprintf('Per-animal available depth: %s\n', mat2str(depthPerMouse'));
fprintf('Mice included: %s\n', strjoin(mouseList, ', '));
fprintf('Group assignment: %s\n', strjoin(groupOf, ', '));
fprintf('============================================================\n');

%% -------------------- pivot to wide format --------------------
varNames = arrayfun(@(k) sprintf('Session%d', k), 1:nSessAnalyze, 'UniformOutput', false);

wideVals = nan(nMice, nSessAnalyze);
for i = 1:nMice
    rows_i = convergenceTable(convMouseId == string(mouseList{i}), :);
    rows_i = rows_i(rows_i.sessFromEnd < nSessAnalyze, :);

    assert(height(rows_i) == nSessAnalyze, ...
        ['Mouse "%s" has %d rows in the analysis window (expected %d) -- check for gaps ' ...
         'in sessFromEnd (this should not happen if convergenceTable came from the current ' ...
         'plotMDSConvergenceAndAmongDistance.m unmodified).'], ...
        mouseList{i}, height(rows_i), nSessAnalyze);

    % sessPos: 1 = earliest session IN THE WINDOW, nSessAnalyze = that
    % mouse's own final session (ascending in calendar/session time).
    sessPos = nSessAnalyze - rows_i.sessFromEnd;
    [~, ord] = sort(sessPos);
    wideVals(i, :) = rows_i.(valueVar)(ord)';
end

wideTable = array2table(wideVals, 'VariableNames', varNames);
wideTable.Group = categorical(groupOf(:));
wideTable = wideTable(:, ['Group', varNames]);
wideTable.Properties.RowNames = mouseList;

%% -------------------- fit repeated-measures model --------------------
withinDesign = table((1:nSessAnalyze)', 'VariableNames', {'Session'});
withinDesign.Session = categorical(withinDesign.Session);

rmFormula = sprintf('%s-%s ~ Group', varNames{1}, varNames{end});
rm = fitrm(wideTable, rmFormula, 'WithinDesign', withinDesign);

ranovaTbl  = ranova(rm);
mauchlyTbl = mauchly(rm);

% Pure between-subjects Group main effect (averaged across all
% nSessAnalyze sessions) -- NOT automatically included in ranova(rm)'s
% default output above, which only reports effects involving the
% within-subject factor (Session main effect, Group:Session interaction).
% Collapsing the within-model to an intercept-only design ('WithinModel','1')
% runs the test on each animal's MEAN across all analyzed sessions, which
% is the classical "is Group different on average, ignoring session"
% question -- this was missing from the original version of this
% function despite being one of the three explicitly requested effects
% (Group, Session, Group x Session).
ranovaBetweenTbl = ranova(rm, 'WithinModel', '1');

multCompareTbl = multcompare(rm, 'Group', 'By', 'Session', 'ComparisonType', char(opt.comparisonType));

%% -------------------- verbose output --------------------
if opt.verbose
    fprintf('\n--- Wide-format data fed to fitrm ---\n');
    disp(wideTable);

    fprintf('\n--- Repeated-measures ANOVA, WITHIN-subject effects (ranova) ---\n');
    fprintf('(Session main effect, Group:Session interaction -- NOT the pure Group main effect; see below)\n');
    disp(ranovaTbl);
    fprintf(['NOTE: check mauchlyTbl below before trusting the uncorrected pValue column above -- ' ...
             'if sphericity is violated (p<0.05), use pValueGG (Greenhouse-Geisser corrected) instead.\n']);

    fprintf('\n--- Repeated-measures ANOVA, BETWEEN-subjects Group main effect ---\n');
    fprintf('(test on each animal''s MEAN across all %d analyzed sessions -- this IS the "main effect of Group")\n', nSessAnalyze);
    disp(ranovaBetweenTbl);

    fprintf('\n--- Mauchly''s test of sphericity ---\n');
    disp(mauchlyTbl);

    fprintf('\n--- Post-hoc: Group comparison at each Session (%s correction) ---\n', char(opt.comparisonType));
    disp(multCompareTbl);
end

%% -------------------- package output --------------------
statsOut = struct();
statsOut.nSessAnalyze   = nSessAnalyze;
statsOut.mouseList      = mouseList;
statsOut.groupOf        = groupOf;
statsOut.depthPerMouse  = depthPerMouse;
statsOut.wideTable      = wideTable;
statsOut.rm             = rm;
statsOut.ranovaTbl        = ranovaTbl;
statsOut.ranovaBetweenTbl = ranovaBetweenTbl;   % pure Group main effect -- see note above
statsOut.mauchlyTbl       = mauchlyTbl;
statsOut.multCompareTbl = multCompareTbl;
statsOut.valueVar       = valueVar;
statsOut.comparisonType = char(opt.comparisonType);

%% -------------------- NEW: optional save to disk --------------------
% MATLAB cannot see what variable name the CALLER assigns this function's
% return value to (e.g. "statsOut_cr = runConvergenceRMANOVA_fastSlow(...)")
% -- inputname() only resolves INPUT argument names, and the assignment
% itself only happens after this function has already returned, so there
% is no way to auto-derive a filename from it. Hence the explicit
% 'saveTag' name-value pair below.
if opt.doSave
    if strlength(strtrim(string(opt.saveDir))) == 0
        error(['doSave is true (the default) but no ''saveDir'' was provided. ' ...
               'Pass ''saveDir'', ''<path>'' (and ideally ''saveTag'', ''<label>'' for a ' ...
               'discernible filename, e.g. ''cr'' or ''combinedCorrect''), or pass ''doSave'', false ' ...
               'to skip saving.']);
    end

    outDir = char(string(opt.saveDir));
    if exist(outDir, 'dir') ~= 7
        mkdir(outDir);
    end

    tag = strtrim(char(string(opt.saveTag)));
    if isempty(tag)
        tag = 'unlabeled';
        warning('runConvergenceRMANOVA_fastSlow:noSaveTag', ...
            ['No ''saveTag'' provided -- saving with tag "unlabeled". Pass e.g. ''saveTag'',''cr'' ' ...
             'so the saved filename actually tells you which analysis it is.']);
    end

    dateStr  = char(datetime('today','Format','MMddyy'));
    saveName = sprintf('convergenceRMANOVA_%s_%s_n%d_%s.mat', tag, valueVar, nSessAnalyze, dateStr);
    saveFullPath = fullfile(outDir, saveName);

    save(saveFullPath, 'statsOut', 'groupDefs');

    statsOut.save = struct('didSave', true, 'file', saveFullPath, 'saveDir', outDir, 'saveTag', tag, 'dateStr', dateStr);
    fprintf('\nSaved stats output to:\n%s\n', saveFullPath);
else
    statsOut.save = struct('didSave', false, 'file', '', 'saveDir', '', 'saveTag', '', 'dateStr', '');
end

end