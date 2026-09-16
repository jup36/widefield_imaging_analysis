function [xcorrRezC, headerC, mIdC, missingC, nullT] = collectMotifXcorrResidualPosLagShuffleResults(filePathBase, varargin)
% COLLECTMOTIFXCORRRESIDUALPOSLAGSHUFFLERESULTS
%   Reassemble per-session results saved by
%   MotifXcorrResidualPosLagShuffle_Spock_ArrayTask (one array task per
%   session, run via Data_GNG_Scotty_MotifXcorrRsdPosLagShuff_ArraySubmit)
%   into the same xcorrRezC / headerC / mIdC animal-by-session cell array
%   shape used elsewhere in this pipeline.
%
%   'doPSTHSubtraction' (default true) selects WHICH variant to collect --
%   true (default, UNCHANGED naming/behavior from every prior version of
%   this collector) looks for "xcorrResidualPosLagShuffle" files
%   (PSTH-subtracted); false looks for "xcorrRawTracePosLagShuffle" files
%   (raw-trace, no PSTH subtraction) instead. The two variants live in the
%   same Matfiles folders but can never cross-match, since their base
%   keywords share no substring. Output FIELD NAMES inside each collected
%   session struct are IDENTICAL regardless of which variant was collected
%   (S.params.doPSTHSubtraction inside each entry records which one it
%   actually is) -- this is what lets every downstream script stay
%   unchanged regardless of which collection gets loaded.
%
%   TWO NULL DISTRIBUTIONS
%   ----------------------
%   The analysis function now computes two independent shuffle nulls,
%   toggled by doTimeShuffle (within-trial circshift/permute) and
%   doTrialShuffle (across-trial permutation). Because this collector
%   stores each session's result struct WHOLE, both nulls are carried into
%   xcorrRezC automatically -- no field-by-field plumbing is needed here.
%   Downstream scripts reach them as:
%       S.shuf.xcorrPosLagMat_<stream>_residual                 (within-trial null)
%       S.shuf.xcorrPosLagMat_<stream>_residual_trialShuffle    (trial-shuffle null)
%   and likewise mean_/std_/z_/p_ and S.shuf.curve.* with the same suffix.
%
%   What this collector DOES now do is AUDIT null availability per session,
%   because results predating the doTrialShuffle addition have no
%   _trialShuffle fields at all. Collecting old and new sessions together
%   would otherwise produce an xcorrRezC where some entries carry the
%   trial-shuffle null and others silently don't -- a downstream loop would
%   then error partway through, or skip those sessions without saying so.
%   The returned nullT table and the printed summary make any such mix
%   explicit BEFORE it reaches an analysis script. 'requireConsistentNulls'
%   (default false) escalates a mix from warning to error.
%
%   USAGE
%   [xcorrRezC, headerC, mIdC, missingC, nullT] = collectMotifXcorrResidualPosLagShuffleResults( ...
%       filePathBase, ...
%       'doPSTHSubtraction',      true, ...   % false for the raw-trace variant
%       'submissionRecord',       '<path>', ...% optional; auto-finds latest matching variant if omitted
%       'requireConsistentNulls', false, ...  % true = error (not warn) if sessions disagree on available nulls
%       'saveCollection',         true, ...   % optional; saves combined .mat when done
%       'fileSaveDir',            '<path>');  % optional; where to save the collection
%
%   OUTPUTS
%   xcorrRezC : {nAnimals x nSessions} cell array of per-session result
%               structs (same fields regardless of doPSTHSubtraction:
%               .meta, .params, .obs, .shuf)
%   headerC   : {nAnimals x nSessions} cell array of session headers
%   mIdC      : {nAnimals x nSessions} cell array of mouse IDs
%   missingC  : table of sessions whose saved result file could not be
%               found/loaded -- check this before trusting xcorrRezC is
%               complete. Does NOT include sessions that loaded fine but
%               had a trial type skipped for thin trial counts -- see the
%               printed skipped-trial-type summary and
%               result.meta.skipped_* for those.
%   nullT     : per-session table of which nulls are present
%               (hasTimeShuffle / hasTrialShuffle, plus the params flags
%               each session recorded). Also saved into the collection.

%% -------------------- Parse inputs --------------------
p = inputParser;
addRequired(p, 'filePathBase', @(x) ischar(x) || isstring(x));
addParameter(p, 'doPSTHSubtraction', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'requireConsistentNulls', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'submissionRecord', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'saveCollection', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'fileSaveDir', '', @(x) ischar(x) || isstring(x));
parse(p, filePathBase, varargin{:});

filePathBase           = char(p.Results.filePathBase);
doPSTHSubtraction      = p.Results.doPSTHSubtraction;
requireConsistentNulls = p.Results.requireConsistentNulls;
submissionRecord       = char(p.Results.submissionRecord);
saveCollection         = p.Results.saveCollection;
fileSaveDir            = char(p.Results.fileSaveDir);

% -------- naming base, conditional on doPSTHSubtraction (mirrors the
% submission wrapper's convention exactly -- must stay in sync with it) --------
if doPSTHSubtraction
    baseKeyword = 'xcorrResidualPosLagShuffle';   % UNCHANGED from every prior version
else
    baseKeyword = 'xcorrRawTracePosLagShuffle';   % no substring overlap with the above
end

manifestDir = fullfile(filePathBase, 'collectData', 'xcorr_withinTrialTimeShuffle', 'arrayManifests');

if isempty(fileSaveDir)
    fileSaveDir = fullfile(filePathBase, 'collectData', 'xcorr_withinTrialTimeShuffle');
end

%% -------------------- Locate submission record --------------------
% CAUTION: when auto-selecting (submissionRecord not supplied), this picks
% the file with the most recent filesystem modification time (datenum).
% That's usually fine, but if this manifest folder now contains submission
% records from more than one nShuffle run (e.g. n100 and n1000), or from
% BOTH doPSTHSubtraction variants (which now live side by side under
% different base keywords), and filePathBase is a network-mounted drive
% (SMB mtimes can occasionally be stale/skewed), it's worth passing
% 'submissionRecord' explicitly instead of trusting auto-detection --
% especially since picking the wrong record would silently collect a
% different run's files without erroring (nShuffle and saveKeyword below
% are both derived from whichever record gets selected). The
% doPSTHSubtraction-based baseKeyword filter below at least guarantees
% auto-detection can never cross the two PSTH-subtraction variants, even
% if it picks an unintended nShuffle/date among files of the SAME variant.
% The auto-selected path is always printed below; confirm it matches the
% run you intend before trusting the output.
if isempty(submissionRecord)
    recFiles = dir(fullfile(manifestDir, [baseKeyword '_submissionRecord_*.mat']));
    if isempty(recFiles)
        error('No submission record found under %s matching baseKeyword "%s" (doPSTHSubtraction=%d). Pass one explicitly via ''submissionRecord''.', ...
            manifestDir, baseKeyword, doPSTHSubtraction);
    end
    [~, mostRecentIdx] = max([recFiles.datenum]);
    submissionRecord = fullfile(manifestDir, recFiles(mostRecentIdx).name);
    fprintf('Using most recent submission record (doPSTHSubtraction=%d):\n%s\n', doPSTHSubtraction, submissionRecord);
end

rec = load(submissionRecord, 'headerC', 'mIdC', 'nSessions', 'funcNV', 'doPSTHSubtraction', 'baseKeyword');

% Cross-check: if an explicit submissionRecord was supplied, confirm it
% actually matches the requested doPSTHSubtraction -- otherwise saveKeyword
% below gets built with the WRONG baseKeyword, and every session would
% silently land in missingC ("result not found") instead of surfacing a
% clear error about the actual mismatch. Older records predating this
% field are allowed through with a warning rather than blocked.
if isfield(rec, 'doPSTHSubtraction')
    if rec.doPSTHSubtraction ~= doPSTHSubtraction
        error(['Mismatch: the loaded submission record was built with doPSTHSubtraction=%d ' ...
               '(baseKeyword="%s"), but this call requested doPSTHSubtraction=%d (baseKeyword="%s"). ' ...
               'Pass the matching submissionRecord, or correct the doPSTHSubtraction argument.'], ...
               rec.doPSTHSubtraction, rec.baseKeyword, doPSTHSubtraction, baseKeyword);
    end
else
    warning(['Loaded submission record predates the doPSTHSubtraction field -- cannot verify it matches ' ...
             'the requested doPSTHSubtraction=%d. Proceeding, but confirm this record is the intended variant.'], ...
             doPSTHSubtraction);
end

headerC_flat = rec.headerC;
mIdC_flat    = rec.mIdC;
nSessions    = rec.nSessions;

nShuffleIdx = find(strcmpi(rec.funcNV, 'nShuffle'));
if isempty(nShuffleIdx)
    error('Could not find nShuffle in submission record''s funcNV.');
end
nShuffle = rec.funcNV{nShuffleIdx + 1};

% What the submission record ASKED for, where recorded -- compared against
% what each session actually contains, further below.
requested_doTimeShuffle  = getNVdefault(rec.funcNV, 'doTimeShuffle',  NaN);
requested_doTrialShuffle = getNVdefault(rec.funcNV, 'doTrialShuffle', NaN);

fprintf('\n============================================================\n');
fprintf('Collecting %d session results (doPSTHSubtraction=%d, baseKeyword=%s)\n', ...
    nSessions, doPSTHSubtraction, baseKeyword);
fprintf('nShuffle: %d\n', nShuffle);
fprintf('Submission record requested: doTimeShuffle=%s, doTrialShuffle=%s\n', ...
    flagStr(requested_doTimeShuffle), flagStr(requested_doTrialShuffle));
fprintf('============================================================\n');

%% -------------------- Group sessions by animal (preserve j/jj shape) --------------------
uniqueAnimals = unique(mIdC_flat, 'stable');
nAnimals = numel(uniqueAnimals);

xcorrRezC = cell(nAnimals, 0);
headerC   = cell(nAnimals, 0);
mIdC      = cell(nAnimals, 0);

sessionCounter = zeros(nAnimals, 1);

missingRows = struct('mId', {}, 'header', {}, 'reason', {});

% Skip tracking now covers all four trial types and the four pooled
% streams, matching what the analysis function actually reports -- the
% previous hit/cr-only version could not flag a session whose FA or Miss
% count sank a pooled stream.
skippedRows = struct('mId', {}, 'header', {}, ...
    'nTrials_hit', {}, 'nTrials_cr', {}, 'nTrials_fa', {}, 'nTrials_miss', {}, ...
    'skipped_hit', {}, 'skipped_cr', {}, 'skipped_fa', {}, 'skipped_miss', {}, ...
    'skipped_combinedCorrect', {}, 'skipped_combinedAll', {}, ...
    'skipped_allGo', {}, 'skipped_allNoGo', {});

% Null-availability audit, one row per successfully loaded session.
nullRows = struct('mId', {}, 'header', {}, ...
    'param_doTimeShuffle', {}, 'param_doTrialShuffle', {}, ...
    'hasTimeShuffle', {}, 'hasTrialShuffle', {}, 'nShuffle_session', {});

%% -------------------- Load each session's saved result --------------------
for i = 1:nSessions
    mId    = mIdC_flat{i};
    header = headerC_flat{i};

    animalIdx = find(strcmp(uniqueAnimals, mId), 1, 'first');
    sessionCounter(animalIdx) = sessionCounter(animalIdx) + 1;
    jj = sessionCounter(animalIdx);

    headerC{animalIdx, jj} = header;
    mIdC{animalIdx, jj}    = mId;

    animalFolderC = find_keyword_containing_folder(filePathBase, mId, 'recursive', false);
    if isempty(animalFolderC)
        missingRows(end+1) = struct('mId', mId, 'header', header, ...
            'reason', 'animal folder not found locally'); %#ok<AGROW>
        continue;
    end
    animalFolder = animalFolderC{1};

    sessionFolderC = find_keyword_containing_folder(animalFolder, header, 'recursive', false);
    if isempty(sessionFolderC)
        missingRows(end+1) = struct('mId', mId, 'header', header, ...
            'reason', 'session folder not found locally'); %#ok<AGROW>
        continue;
    end

    taskFolderC = find_keyword_containing_folder(sessionFolderC{1}, 'task', 'recursive', false);
    if isempty(taskFolderC)
        missingRows(end+1) = struct('mId', mId, 'header', header, ...
            'reason', 'task folder not found locally'); %#ok<AGROW>
        continue;
    end

    filePath_mat = cell2mat(find_keyword_containing_folder(taskFolderC{1}, 'Matfiles', 'recursive', false));
    if isempty(filePath_mat)
        missingRows(end+1) = struct('mId', mId, 'header', header, ...
            'reason', 'Matfiles folder not found locally'); %#ok<AGROW>
        continue;
    end

    % Match by header + baseKeyword + nShuffle (NOT exact date stamp) --
    % robust to sessions (re)submitted on a different day, and the
    % baseKeyword filter guarantees this can never cross-match the other
    % doPSTHSubtraction variant's files even though they live in the same
    % Matfiles folder.
    saveKeyword = sprintf('%s_%s_n%d_', header, baseKeyword, nShuffle);
    resultFileC = find_keyword_containing_files(filePath_mat, saveKeyword, 'recursive', false);

    if isempty(resultFileC)
        missingRows(end+1) = struct('mId', mId, 'header', header, ...
            'reason', 'saved result .mat not found (task may not have run/finished yet, or failed)'); %#ok<AGROW>
        continue;
    end

    if numel(resultFileC) > 1
        % NOTE: with two nulls now in play, a stale pre-doTrialShuffle file
        % and a fresh one can both match this keyword. Most-recent wins,
        % which is what you want, but the warning is worth reading.
        fileInfo = cellfun(@(f) dir(f), resultFileC);
        [~, mostRecentIdx] = max([fileInfo.datenum]);
        warning('Multiple saved result files found for %s; using most recent: %s', ...
            header, resultFileC{mostRecentIdx});
        resultFile = resultFileC{mostRecentIdx};
    else
        resultFile = resultFileC{1};
    end

    try
        S = load(resultFile, 'result');
        xcorrRezC{animalIdx, jj} = S.result;
        fprintf('[%3d/%3d] Loaded %s  %s\n', i, nSessions, mId, header);

        % ---- skip tracking (all trial types + pooled streams) ----
        if isfield(S.result, 'meta')
            m = S.result.meta;
            skipFlags = [getfielddefault(m, 'skipped_hit',  false), ...
                         getfielddefault(m, 'skipped_cr',   false), ...
                         getfielddefault(m, 'skipped_fa',   false), ...
                         getfielddefault(m, 'skipped_miss', false), ...
                         getfielddefault(m, 'skipped_combinedCorrect', false), ...
                         getfielddefault(m, 'skipped_combinedAll',     false), ...
                         getfielddefault(m, 'skipped_allGo',   false), ...
                         getfielddefault(m, 'skipped_allNoGo', false)];
            if any(skipFlags)
                skippedRows(end+1) = struct('mId', mId, 'header', header, ...
                    'nTrials_hit',  getfielddefault(m, 'nTrials_hit',  NaN), ...
                    'nTrials_cr',   getfielddefault(m, 'nTrials_cr',   NaN), ...
                    'nTrials_fa',   getfielddefault(m, 'nTrials_fa',   NaN), ...
                    'nTrials_miss', getfielddefault(m, 'nTrials_miss', NaN), ...
                    'skipped_hit',  skipFlags(1), 'skipped_cr',   skipFlags(2), ...
                    'skipped_fa',   skipFlags(3), 'skipped_miss', skipFlags(4), ...
                    'skipped_combinedCorrect', skipFlags(5), ...
                    'skipped_combinedAll',     skipFlags(6), ...
                    'skipped_allGo',   skipFlags(7), ...
                    'skipped_allNoGo', skipFlags(8)); %#ok<AGROW>
            end
        end

        % ---- null-availability audit ----
        % Presence is checked on the FIELDS, not just the params flags: a
        % session could have been submitted with doTrialShuffle=true but
        % produced by an older function build that ignored it. The Hit
        % stream is the probe, since it is present in every complete run.
        prm = getfielddefault(S.result, 'params', struct());
        shf = getfielddefault(S.result, 'shuf',   struct());

        nullRows(end+1) = struct('mId', mId, 'header', header, ...
            'param_doTimeShuffle',  getfielddefault(prm, 'doTimeShuffle',  NaN), ...
            'param_doTrialShuffle', getfielddefault(prm, 'doTrialShuffle', NaN), ...
            'hasTimeShuffle',  isfield(shf, 'xcorrPosLagMat_hit_residual'), ...
            'hasTrialShuffle', isfield(shf, 'xcorrPosLagMat_hit_residual_trialShuffle'), ...
            'nShuffle_session', getfielddefault(prm, 'nShuffle', NaN)); %#ok<AGROW>

        % Sanity check: confirm the loaded session's own recorded
        % doPSTHSubtraction flag matches what we intended to collect --
        % catches a mismatched/corrupted file rather than silently mixing
        % variants into one collection.
        if isfield(prm, 'doPSTHSubtraction') && prm.doPSTHSubtraction ~= doPSTHSubtraction
            warning('[%s %s] Loaded file''s doPSTHSubtraction=%d does not match requested %d -- check for a naming/collection mismatch.', ...
                mId, header, prm.doPSTHSubtraction, doPSTHSubtraction);
        end
    catch ME
        missingRows(end+1) = struct('mId', mId, 'header', header, ...
            'reason', sprintf('failed to load: %s', ME.message)); %#ok<AGROW>
    end
end

if isempty(missingRows)
    missingC = table('Size', [0 3], 'VariableTypes', {'cell','cell','cell'}, ...
        'VariableNames', {'mId','header','reason'});
else
    missingC = struct2table(missingRows, 'AsArray', true);
end

if isempty(nullRows)
    nullT = table();
else
    nullT = struct2table(nullRows, 'AsArray', true);
end

%% -------------------- Report --------------------
nLoaded = sum(~cellfun(@isempty, xcorrRezC(:)));
fprintf('\n============================================================\n');
fprintf('Loaded %d/%d sessions successfully.\n', nLoaded, nSessions);
if ~isempty(missingRows)
    fprintf('%d session(s) missing/failed -- see returned missingC table:\n', numel(missingRows));
    disp(missingC);
else
    fprintf('No missing sessions.\n');
end

% -------- null availability --------
fprintf('\n---- Shuffle null availability ----\n');
if isempty(nullT)
    fprintf('No sessions loaded -- nothing to audit.\n');
else
    nTime  = sum(nullT.hasTimeShuffle);
    nTrial = sum(nullT.hasTrialShuffle);
    nBoth  = sum(nullT.hasTimeShuffle & nullT.hasTrialShuffle);
    nNone  = sum(~nullT.hasTimeShuffle & ~nullT.hasTrialShuffle);

    fprintf('  within-trial null present : %d/%d sessions\n', nTime,  height(nullT));
    fprintf('  trial-shuffle null present: %d/%d sessions\n', nTrial, height(nullT));
    fprintf('  both nulls present        : %d/%d sessions\n', nBoth,  height(nullT));
    if nNone > 0
        fprintf('  NEITHER null present      : %d session(s) (observed-only runs)\n', nNone);
    end

    mixedTime  = nTime  > 0 && nTime  < height(nullT);
    mixedTrial = nTrial > 0 && nTrial < height(nullT);

    if mixedTime || mixedTrial
        badT = nullT(~nullT.hasTimeShuffle | ~nullT.hasTrialShuffle, ...
            {'mId','header','hasTimeShuffle','hasTrialShuffle'});
        msg = sprintf(['Sessions in this collection DISAGREE on which shuffle nulls are present ' ...
            '(within-trial: %d/%d, trial-shuffle: %d/%d). Results predating the doTrialShuffle ' ...
            'addition have no _trialShuffle fields. A downstream loop over sessions will error ' ...
            'or silently skip on the incomplete ones. Sessions lacking a null are listed below.'], ...
            nTime, height(nullT), nTrial, height(nullT));
        if requireConsistentNulls
            disp(badT);
            error('%s\nCalled with requireConsistentNulls=true -- aborting rather than saving a mixed collection.', msg);
        else
            warning('%s', msg);
            disp(badT);
        end
    else
        fprintf('  Consistent across all loaded sessions.\n');
    end

    % nShuffle should match the record; a mismatch means a stale file slipped through
    nSh = unique(nullT.nShuffle_session(isfinite(nullT.nShuffle_session)));
    if numel(nSh) > 1
        warning('Sessions report differing nShuffle values (%s) -- a stale result file may have been matched.', ...
            mat2str(nSh(:)'));
    elseif ~isempty(nSh) && nSh ~= nShuffle
        warning('Sessions report nShuffle=%d but the submission record says %d.', nSh, nShuffle);
    end
end

% -------- skipped trial types --------
if ~isempty(skippedRows)
    skippedT = struct2table(skippedRows, 'AsArray', true);
    fprintf('\n%d loaded session(s) had at least one trial type or pooled stream SKIPPED:\n', ...
        numel(skippedRows));
    disp(skippedT);
    fprintf(['NOTE: these sessions ARE included in xcorrRezC (file loaded fine), but the skipped ' ...
             'stream''s matrices are NaN. Check result.meta.skipped_* before pooling across sessions.\n']);
else
    skippedT = table();
    fprintf('\nNo loaded sessions had a skipped trial type or pooled stream.\n');
end
fprintf('============================================================\n');

%% -------------------- Save combined collection --------------------
if saveCollection
    if exist(fileSaveDir, 'dir') ~= 7
        mkdir(fileSaveDir);
    end
    dateStr = string(datetime('today','Format','MMddyy'));
    % Output collection filename baseKeyword-tagged, mirroring the
    % submission side -- true (default) keeps the EXACT prior naming
    % ("xcorrResidualPosLag_motifH_collect_..."), so nothing about
    % already-collected PSTH-subtracted results changes; false gets a
    % distinctly different, non-overlapping name. Deliberately NOT tagged
    % with which nulls are present: that varies per session and belongs in
    % nullT, not in a filename that downstream load paths depend on.
    if doPSTHSubtraction
        collectionBase = "xcorrResidualPosLag_motifH_collect_timeshuffle";
    else
        collectionBase = "xcorrRawTracePosLag_motifH_collect_timeshuffle";
    end
    saveName = collectionBase + "_n" + string(nShuffle) + "_" + dateStr + ".mat";
    saveFullPath = fullfile(fileSaveDir, saveName);

    % nullT saved alongside so the audit travels with the data -- a script
    % loading this collection months from now can check which nulls it has
    % without re-deriving it from the structs.
    save(saveFullPath, "xcorrRezC", "mIdC", "headerC", "missingC", "nullT", ...
        "skippedT", "doPSTHSubtraction", '-v7.3');

    fprintf('\nSaved combined collection:\n%s\n', saveFullPath);
end

end

%% ========================================================================
function v = getfielddefault(s, fieldName, defaultVal)
if isstruct(s) && isfield(s, fieldName)
    v = s.(fieldName);
else
    v = defaultVal;
end
end

%% ========================================================================
function v = getNVdefault(nvCell, name, defaultVal)
% Pull a value out of a {'name', value, ...} cell, or return a default.
idx = find(strcmpi(nvCell, name), 1, 'first');
if isempty(idx) || idx == numel(nvCell)
    v = defaultVal;
else
    v = nvCell{idx + 1};
end
end

%% ========================================================================
function s = flagStr(v)
if isnumeric(v) && isnan(v)
    s = 'not recorded';
else
    s = sprintf('%d', v);
end
end