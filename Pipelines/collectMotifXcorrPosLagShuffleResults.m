function [xcorrRezC, headerC, mIdC, missingC] = collectMotifXcorrPosLagShuffleResults(filePathBase, varargin)
% COLLECTMOTIFXCORRPOSLAGSHUFFLERESULTS
%   Reassemble per-session results saved by MotifXcorrPosLagShuffle_Spock_ArrayTask
%   (one array task per session, run via
%   Data_Dual_GNG_Scotty_MotifXcorrPosLagShuffle_ArraySubmit) into the
%   same xcorrRezC / headerC / mIdC animal-by-session cell array shape
%   that the original serial batch script produced in-memory, so the
%   example plots at the bottom of that script (post-lag xcorr heatmap,
%   motif-pair xcorr with shuffle CI, etc.) work unchanged.
%
%   Rather than converting bucket paths (e.g. /jukebox/buschman/...) back
%   to local Windows paths -- which has proven fragile throughout this
%   whole pipeline (mount-point mismatches, symlink surprises) -- this
%   re-locates each session's Matfiles folder locally the same way the
%   submitter's own discovery loop did: by mId, then header, via
%   find_keyword_containing_folder. That's the one lookup mechanism we've
%   already confirmed works correctly end to end.
%
%   USAGE
%   [xcorrRezC, headerC, mIdC, missingC] = collectMotifXcorrPosLagShuffleResults( ...
%       filePathBase, ...
%       'submissionRecord', '<path>', ...   % optional; auto-finds latest if omitted
%       'saveCollection',   true, ...       % optional; saves combined .mat when done
%       'fileSaveDir',      '<path>');      % optional; where to save the collection
%
%   OUTPUTS
%   xcorrRezC : {nAnimals x nSessions} cell array of per-session result
%               structs (same fields as returned by
%               motifH_perTrial_crosscorr_posLag_withinTrialTimeShuffle_func:
%               .meta, .params, .obs, .shuf)
%   headerC   : {nAnimals x nSessions} cell array of session headers
%   mIdC      : {nAnimals x nSessions} cell array of mouse IDs
%   missingC  : table of sessions whose saved result file could not be
%               found (e.g. still running, or the array task failed) --
%               check this before trusting xcorrRezC is complete

%% -------------------- Parse inputs --------------------
p = inputParser;
addRequired(p, 'filePathBase', @(x) ischar(x) || isstring(x));
addParameter(p, 'submissionRecord', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'saveCollection', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'fileSaveDir', '', @(x) ischar(x) || isstring(x));
parse(p, filePathBase, varargin{:});

filePathBase      = char(p.Results.filePathBase);
submissionRecord  = char(p.Results.submissionRecord);
saveCollection    = p.Results.saveCollection;
fileSaveDir       = char(p.Results.fileSaveDir);

manifestDir = fullfile(filePathBase, 'collectData', 'xcorr_withinTrialTimeShuffle', 'arrayManifests');

if isempty(fileSaveDir)
    fileSaveDir = fullfile(filePathBase, 'collectData', 'xcorr_withinTrialTimeShuffle');
end

%% -------------------- Locate submission record --------------------
if isempty(submissionRecord)
    recFiles = dir(fullfile(manifestDir, 'xcorrPosLagShuffle_submissionRecord_*.mat'));
    if isempty(recFiles)
        error('No submission record found under %s. Pass one explicitly via ''submissionRecord''.', manifestDir);
    end
    [~, mostRecentIdx] = max([recFiles.datenum]);
    submissionRecord = fullfile(manifestDir, recFiles(mostRecentIdx).name);
    fprintf('Using most recent submission record:\n%s\n', submissionRecord);
end

rec = load(submissionRecord, 'headerC', 'mIdC', 'nSessions', 'funcNV');

headerC_flat = rec.headerC;
mIdC_flat    = rec.mIdC;
nSessions    = rec.nSessions;

% Extract nShuffle from funcNV (name-value cell array) for constructing
% the expected save filename pattern.
nShuffleIdx = find(strcmpi(rec.funcNV, 'nShuffle'));
if isempty(nShuffleIdx)
    error('Could not find nShuffle in submission record''s funcNV.');
end
nShuffle = rec.funcNV{nShuffleIdx + 1};

fprintf('\n============================================================\n');
fprintf('Collecting %d session results\n', nSessions);
fprintf('nShuffle: %d\n', nShuffle);
fprintf('============================================================\n');

%% -------------------- Group sessions by animal (preserve j/jj shape) --------------------
uniqueAnimals = unique(mIdC_flat, 'stable');
nAnimals = numel(uniqueAnimals);

xcorrRezC = cell(nAnimals, 0);
headerC   = cell(nAnimals, 0);
mIdC      = cell(nAnimals, 0);

sessionCounter = zeros(nAnimals, 1);

missingRows = struct('mId', {}, 'header', {}, 'reason', {});

%% -------------------- Load each session's saved result --------------------
for i = 1:nSessions
    mId    = mIdC_flat{i};
    header = headerC_flat{i};

    animalIdx = find(strcmp(uniqueAnimals, mId), 1, 'first');
    sessionCounter(animalIdx) = sessionCounter(animalIdx) + 1;
    jj = sessionCounter(animalIdx);

    headerC{animalIdx, jj} = header;
    mIdC{animalIdx, jj}    = mId;

    % Re-locate this session's Matfiles folder locally (mirrors the
    % submitter's own discovery loop -- the one lookup mechanism we've
    % already confirmed works correctly).
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

    % Match by the stable part of the filename (header + nShuffle count),
    % not the exact date stamp -- robust to sessions having been
    % (re)submitted on a different day than others.
    saveKeyword = sprintf('%s_xcorrPosLagShuffle_n%d_', header, nShuffle);
    resultFileC = find_keyword_containing_files(filePath_mat, saveKeyword, 'recursive', false);

    if isempty(resultFileC)
        missingRows(end+1) = struct('mId', mId, 'header', header, ...
            'reason', 'saved result .mat not found (task may not have run/finished yet, or failed)'); %#ok<AGROW>
        continue;
    end

    if numel(resultFileC) > 1
        % More than one match (e.g. re-run on a different date) -- take
        % the most recently modified one and warn.
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
fprintf('============================================================\n');

%% -------------------- Save combined collection --------------------
if saveCollection
    if exist(fileSaveDir, 'dir') ~= 7
        mkdir(fileSaveDir);
    end
    dateStr  = string(datetime('today','Format','MMddyy'));
    saveName = "xcorrPosLag_motifH_collect_timeshuffle_" + dateStr + ".mat";
    saveFullPath = fullfile(fileSaveDir, saveName);

    save(saveFullPath, "xcorrRezC", "mIdC", "headerC", "missingC", '-v7.3');
    fprintf('\nSaved combined collection:\n%s\n', saveFullPath);
end

end
