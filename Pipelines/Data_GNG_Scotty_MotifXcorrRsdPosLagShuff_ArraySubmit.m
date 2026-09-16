function Data_GNG_Scotty_MotifXcorrRsdPosLagShuff_ArraySubmit(filePathBase, fileKeyword, varargin)
% Data_GNG_Scotty_MotifXcorrRsdPosLagShuff_ArraySubmit
%
% Discover all valid sessions under filePathBase, build one SLURM array
% manifest entry per session, and submit a single array job. Each array
% task runs motifH_perTrial_xcorr_residual_posLag_withinTrialTimeShuff_func
% for its session and SAVES the result into that session's own Matfiles
% folder.
%
% NEW: 'doPSTHSubtraction' (default true) is passed straight through to
% the core function -- true = PSTH-subtracted ("residual") xcorr, exactly
% the prior behavior; false = xcorr on RAW per-trial traces, everything
% else identical. See that function's header for the full rationale.
%
% NAMING (IMPORTANT -- backward compatibility):
%   doPSTHSubtraction = true  (default): saved filenames, manifest name,
%     job stem, and submission-record name are UNCHANGED from every prior
%     version of this script -- this run is indistinguishable on disk from
%     any already-completed PSTH-subtracted Scotty run, and will collide/
%     overwrite correctly (same as always) if resubmitted with the same
%     nShuffle+date.
%   doPSTHSubtraction = false (NEW): base keyword changes from
%     "xcorrResidualPosLagShuffle" to "xcorrRawTracePosLagShuffle"
%     throughout (saved filenames, manifest, job stem, submission record)
%     -- chosen specifically so there is NO substring overlap with the
%     PSTH-subtracted naming (avoids any risk of the collector's
%     substring-based file matching cross-matching the two variants).
%
% Required helper on Scotty MATLAB path:
%   MotifXcorrResidualPosLagShuffle_Spock_ArrayTask.m (UNCHANGED -- it
%   already passes funcNV{:} through generically, and save_fn is fully
%   constructed here before being handed to it via the manifest, so no
%   changes are needed there for this new parameter to take effect).
%
% USAGE
%   Data_GNG_Scotty_MotifXcorrRsdPosLagShuff_ArraySubmit( ...
%       filePathBase, fileKeyword, ...
%       'fileKeywordBeh',        '_alignedPupilOrofacial.mat', ...
%       'excludeAnimals',        {'m1893','m1897'}, ...
%       'xcorrLagWindowSec',     1.0, ...
%       'posLagPoolWindowSec',   0.5, ...
%       'doTimeShuffle',         true, ...
%       'nShuffle',              1000, ...
%       'shuffleMethod',         'circshift', ...
%       'rngSeed',               1, ...
%       'zscoreHs',              true, ...
%       'useSymmetry',           true, ...
%       'doFisherZ',             true, ...
%       'clipR',                 0.999, ...
%       'minCorrectTrials',      10, ...
%       'doPSTHSubtraction',     true, ...   % NEW -- false for the raw-trace variant
%       'showProgress',          true, ...
%       'progressEvery',         10, ...
%       'sbatch_time',           3600, ...
%       'sbatch_memory',         16, ...
%       'sbatch_cpus',           12, ...
%       'array_throttle',        10);

%% -------------------- Parse inputs --------------------
p = inputParser;

addRequired(p, 'filePathBase', @(x) ischar(x) || isstring(x));
addRequired(p, 'fileKeyword',  @(x) ischar(x) || isstring(x));

addParameter(p, 'fileKeywordBeh', '_alignedPupilOrofacial.mat', @(x) ischar(x) || isstring(x));
addParameter(p, 'excludeAnimals', {'m1893','m1897'}, @iscell);

addParameter(p, 'xcorrLagWindowSec',   1.0,  @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'posLagPoolWindowSec', 0.5,  @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'doTimeShuffle',       true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'nShuffle',            100,  @(x) isnumeric(x) && isscalar(x) && x>=1);
addParameter(p, 'shuffleMethod',       'circshift', @(s) ischar(s) || isstring(s));
addParameter(p, 'rngSeed',             1,    @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'zscoreHs',            true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'useSymmetry',         true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'doFisherZ',           true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'clipR',               0.999, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'minCorrectTrials',    10,   @(x) isnumeric(x) && isscalar(x) && x>=1);

% NEW: PSTH-subtraction toggle, passed straight through to the core function.
addParameter(p, 'doPSTHSubtraction',   true, @(x) islogical(x) && isscalar(x));

addParameter(p, 'showProgress',        true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'progressEvery',       10,   @(x) isnumeric(x) && isscalar(x));

addParameter(p, 'sbatch_time',    3600, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'sbatch_memory',  16,   @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'sbatch_cpus',    12,   @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'array_throttle', 10,   @(x) isnumeric(x) && isscalar(x) && x>0);

parse(p, filePathBase, fileKeyword, varargin{:});

filePathBase    = char(p.Results.filePathBase);
fileKeyword     = char(p.Results.fileKeyword);
fileKeywordBeh  = char(p.Results.fileKeywordBeh);
excludeAnimals  = p.Results.excludeAnimals;

sbatch_time     = p.Results.sbatch_time;
sbatch_memory   = p.Results.sbatch_memory;
sbatch_cpus     = p.Results.sbatch_cpus;
array_throttle  = p.Results.array_throttle;

doPSTHSubtraction = p.Results.doPSTHSubtraction;

% -------- naming base, conditional on doPSTHSubtraction (see header note) --------
if doPSTHSubtraction
    baseKeyword = 'xcorrResidualPosLagShuffle';   % UNCHANGED from every prior version
else
    baseKeyword = 'xcorrRawTracePosLagShuffle';   % NEW, no substring overlap with the above
end

funcNV = { ...
    'xcorrLagWindowSec',   p.Results.xcorrLagWindowSec, ...
    'posLagPoolWindowSec', p.Results.posLagPoolWindowSec, ...
    'doTimeShuffle',       p.Results.doTimeShuffle, ...
    'nShuffle',            p.Results.nShuffle, ...
    'shuffleMethod',       p.Results.shuffleMethod, ...
    'rngSeed',             p.Results.rngSeed, ...
    'zscoreHs',            p.Results.zscoreHs, ...
    'useSymmetry',         p.Results.useSymmetry, ...
    'doFisherZ',           p.Results.doFisherZ, ...
    'clipR',               p.Results.clipR, ...
    'minCorrectTrials',    p.Results.minCorrectTrials, ...
    'doPSTHSubtraction',   doPSTHSubtraction, ...
    'showProgress',        p.Results.showProgress, ...
    'progressEvery',       p.Results.progressEvery};

fprintf('\n============================================================\n');
fprintf('Motif xcorr pos-lag time-shuffle ARRAY submission\n');
fprintf('filePathBase: %s\n', filePathBase);
fprintf('fileKeyword:  %s\n', fileKeyword);
fprintf('nShuffle:     %d\n', p.Results.nShuffle);
fprintf('minCorrectTrials: %d\n', p.Results.minCorrectTrials);
fprintf('doPSTHSubtraction: %d  (baseKeyword: %s)\n', doPSTHSubtraction, baseKeyword);
fprintf('sbatch_time/memory/cpus: %d / %d / %d\n', sbatch_time, sbatch_memory, sbatch_cpus);
fprintf('============================================================\n');

%% -------------------- Discover sessions --------------------
mListC_jRGECO = GrabFiles_sort_trials('m*jRGECO', 0, {filePathBase});
mListC = mListC_jRGECO;

mListI = cellfun(@(a) all(~contains(a, excludeAnimals)), mListC);
mListC = mListC(mListI);

fprintf('\nFound %d animal folders after exclusion.\n', numel(mListC));
disp(mListC(:));

sessionRows = struct('mId', {}, 'header', {}, 'filePath', {}, 'filePath_mat', {});

for j = 1:numel(mListC)
    mId = cell2mat(regexp(mListC{j}, 'm\d{4}', 'match'));

    filePathSessions = find_keyword_containing_folder(mListC{j}, mId, 'recursive', false);

    for jj = 1:numel(filePathSessions)
        filePathC = find_keyword_containing_folder(filePathSessions{jj}, 'task', 'recursive', false);
        if isempty(filePathC); continue; end
        filePath = filePathC{1};

        filePath_mat = cell2mat(find_keyword_containing_folder(filePath, 'Matfiles', 'recursive', false));
        if isempty(filePath_mat); continue; end

        filePath_refit = find_keyword_containing_files(filePath_mat, fileKeyword, 'recursive', true);
        if isempty(filePath_refit); continue; end

        filePath_beh = find_keyword_containing_files(filePath_mat, fileKeywordBeh, 'recursive', false);
        if isempty(filePath_beh); continue; end

        header = extract_date_animalID_header(filePath);

        sessionRows(end+1) = struct( ...
            'mId',          mId, ...
            'header',       header, ...
            'filePath',     filePath, ...
            'filePath_mat', filePath_mat); %#ok<AGROW>
    end
end

nSessions = numel(sessionRows);
fprintf('\nDiscovered %d valid sessions with curated motif + behavioral files.\n', nSessions);

if nSessions == 0
    error('No valid sessions found under %s with keyword %s.', filePathBase, fileKeyword);
end

for s = 1:nSessions
    fprintf('  [%d] %s  %s\n', s, sessionRows(s).mId, sessionRows(s).header);
end

%% -------------------- Connect to Scotty --------------------
keyFile = fullfile(getenv('USERPROFILE'), '.ssh', 'id_ed25519_scotty_matlab');
s_conn = ssh2_command_scotty('connect', 'scotty', keyFile);

%% -------------------- Build manifest --------------------
dateStr = string(datetime('today','Format','MMddyy'));

filePath_bucket     = cell(nSessions, 1);
filePath_mat_bucket = cell(nSessions, 1);
save_fn_bucket       = cell(nSessions, 1);
mIdC                = cell(nSessions, 1);
headerC             = cell(nSessions, 1);

for s = 1:nSessions
    filePath_bucket{s}     = ConvertWinToBucketPath(sessionRows(s).filePath);
    filePath_mat_bucket{s} = ConvertWinToBucketPath(sessionRows(s).filePath_mat);
    mIdC{s}    = sessionRows(s).mId;
    headerC{s} = sessionRows(s).header;

    saveName = sprintf('%s_%s_n%d_%s.mat', ...
        sessionRows(s).header, baseKeyword, p.Results.nShuffle, dateStr);
    save_fn_bucket{s} = ConvertWinToBucketPath(fullfile(sessionRows(s).filePath_mat, saveName));
end

manifestDir = fullfile(filePathBase, 'collectData', 'xcorr_withinTrialTimeShuffle', 'arrayManifests');
if exist(manifestDir, 'dir') ~= 7
    mkdir(manifestDir);
end

manifestName = sprintf('%s_arrayManifest_n%d_%s%s', ...
    baseKeyword, p.Results.nShuffle, dateStr, fileKeyword);
manifest_fn = fullfile(manifestDir, manifestName);

save(manifest_fn, ...
    'filePath_bucket', ...
    'filePath_mat_bucket', ...
    'save_fn_bucket', ...
    'mIdC', ...
    'headerC', ...
    'fileKeyword', ...
    'funcNV', ...
    'nSessions', ...
    '-v7.3');

fprintf('\nSaved array manifest:\n%s\n', manifest_fn);

manifest_fn_bucket = ConvertWinToBucketPath(manifest_fn);

%% -------------------- Write + patch sbatch array script --------------------
if doPSTHSubtraction
    jobTag = 'xcorrResidShuf';
else
    jobTag = 'xcorrRawTraceShuf';
end
jobStem = sanitize_for_slurm(sprintf('%s_n%d_%s_ARRAY', jobTag, p.Results.nShuffle, dateStr));

script_name = WriteBashScriptWinScotty(jobStem, ...
    'MotifXcorrResidualPosLagShuffle_Spock_ArrayTask', ...
    {manifest_fn_bucket}, ...
    {"'%s'"}, ...
    'sbatch_time', sbatch_time, ...
    'sbatch_memory', sbatch_memory, ...
    'sbatch_name', jobStem, ...
    'sbatch_path', ...
      "/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/Spock/");

arrayDirective  = sprintf('#SBATCH --array=1-%d%%%d', nSessions, array_throttle);
cpuDirective    = sprintf('#SBATCH --cpus-per-task=%d', sbatch_cpus);
outputDirective = sprintf('#SBATCH --output=out/%s_%%A_%%a.out', jobStem);

dynamicScriptsDir     = '/jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts';
scriptFullPath_esc    = sprintf('%s/%s', dynamicScriptsDir, script_name);

fprintf('\nEnsuring remote out/ directory exists...\n');
respMkdir = ssh2_command_scotty(s_conn, sprintf('mkdir -p %s/out', dynamicScriptsDir));
if isfield(respMkdir, 'command_result') && ~isempty(respMkdir.command_result)
    disp(respMkdir.command_result);
end

fprintf('\nFetching generated script for local patching:\n%s\n', scriptFullPath_esc);
respCat = ssh2_command_scotty(s_conn, sprintf('cat %s', scriptFullPath_esc));

if ~isfield(respCat, 'command_result') || isempty(respCat.command_result)
    error('Could not read back generated script from Scotty: %s', scriptFullPath_esc);
end

origLines = respCat.command_result(:);
if ~iscell(origLines)
    origLines = cellstr(origLines);
end

shebangIdx = find(startsWith(strtrim(origLines), '#!'), 1, 'first');
if isempty(shebangIdx)
    error(['Could not locate a shebang ("#!...") line in the fetched script -- ' ...
           'cannot safely distinguish real content from SSH banner output. ' ...
           'Raw response was:\n%s'], strjoin(origLines, newline));
end
if shebangIdx > 1
    fprintf('Discarding %d banner/MOTD line(s) prepended by the SSH wrapper before the script content.\n', shebangIdx - 1);
end
origLines = origLines(shebangIdx:end);

newLines = origLines;
newLines = upsertSbatchDirective(newLines, '#SBATCH --array=',           arrayDirective);
newLines = upsertSbatchDirective(newLines, '#SBATCH --cpus-per-task=',   cpuDirective);
newLines = upsertSbatchDirective(newLines, '#SBATCH --output=',          outputDirective);

patchedContent = strjoin(newLines, newline);
if ~endsWith(patchedContent, newline)
    patchedContent = [char(patchedContent), newline];
end

b64Content = matlab.net.base64encode(uint8(char(patchedContent)));
writeCmd = sprintf('printf ''%%s'' ''%s'' | base64 -d > %s', b64Content, scriptFullPath_esc);

fprintf('Writing patched script back to Scotty (base64 round-trip)...\n');
respWrite = ssh2_command_scotty(s_conn, writeCmd);
if isfield(respWrite, 'command_result') && ~isempty(respWrite.command_result)
    disp(respWrite.command_result);
end

respVerify = ssh2_command_scotty(s_conn, sprintf('cat %s', scriptFullPath_esc));
if ~isfield(respVerify, 'command_result') || isempty(respVerify.command_result)
    error('Could not re-read script from Scotty to verify the patch.');
end
verifyLines = respVerify.command_result(:);
if ~iscell(verifyLines)
    verifyLines = cellstr(verifyLines);
end

hasArray  = any(contains(verifyLines, '#SBATCH --array='));
hasCpus   = any(contains(verifyLines, sprintf('#SBATCH --cpus-per-task=%d', sbatch_cpus)));
hasOutput = any(contains(verifyLines, '%A_%a.out'));

if ~(hasArray && hasCpus && hasOutput)
    fprintf(2, '\nPatch verification FAILED. Current remote script contents:\n');
    disp(strjoin(verifyLines, newline));
    error(['Patched sbatch script is missing expected directives ' ...
           '(array=%d, cpus=%d, output=%d). Aborting before submission.'], ...
           hasArray, hasCpus, hasOutput);
end

fprintf('Verified patched script contains --array, --cpus-per-task=%d, and --output=...%%A_%%a.out.\n', ...
    sbatch_cpus);

%% -------------------- Submit array job --------------------
remoteCmd = ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts ; ' ...
    sprintf('sbatch %s', script_name)];

resp = ssh2_command_scotty(s_conn, remoteCmd);
disp(resp.command_result{end});

jidLine = resp.command_result(contains(resp.command_result, "Submitted batch job"));
if isempty(jidLine)
    error('No array job ID found in sbatch response.');
end

array_job_id = extractAfter(jidLine, "Submitted batch job ");
array_job_id = regexprep(array_job_id, '[^0-9]', '');

fprintf('\nSubmitted ARRAY job %s for %d sessions (doPSTHSubtraction=%d).\n', array_job_id, nSessions, doPSTHSubtraction);
fprintf('Array throttle: %d concurrent tasks max.\n', array_throttle);
fprintf('Array directive: %s\n', arrayDirective);
fprintf('CPU directive:   %s\n', cpuDirective);

%% -------------------- Save local submission record --------------------
swarm_id = cell(nSessions, 1);
for s = 1:nSessions
    swarm_id{s} = sprintf('%s_%d', array_job_id, s);
end

recordName = sprintf('%s_submissionRecord_n%d_%s.mat', baseKeyword, p.Results.nShuffle, dateStr);
save(fullfile(manifestDir, recordName), ...
    'array_job_id', 'array_throttle', 'arrayDirective', 'cpuDirective', ...
    'script_name', 'manifest_fn', 'manifest_fn_bucket', ...
    'save_fn_bucket', 'filePath_bucket', 'filePath_mat_bucket', ...
    'mIdC', 'headerC', 'nSessions', 'swarm_id', 'funcNV', 'fileKeyword', ...
    'doPSTHSubtraction', 'baseKeyword', ...
    'sbatch_time', 'sbatch_memory', 'sbatch_cpus', '-v7.3');

fprintf('\nSaved submission record:\n%s\n', fullfile(manifestDir, recordName));

clearvars s_conn

end

%% ========================================================================
function linesOut = upsertSbatchDirective(linesIn, keyPrefix, fullLine)
matchIdx = find(startsWith(strtrim(linesIn), keyPrefix), 1, 'first');

if ~isempty(matchIdx)
    linesOut = linesIn;
    linesOut{matchIdx} = fullLine;
    return;
end

firstSbatchIdx = find(startsWith(strtrim(linesIn), '#SBATCH'), 1, 'first');
if isempty(firstSbatchIdx)
    error('upsertSbatchDirective: no existing "#SBATCH" line found to insert %s after.', keyPrefix);
end

linesOut = [linesIn(1:firstSbatchIdx); {fullLine}; linesIn(firstSbatchIdx+1:end)];
end
