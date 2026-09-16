function Data_Dual_GNG_Scotty_MotifXcorrPosLagShuffle_ArraySubmit(filePathBase, fileKeyword, varargin)
% Data_Dual_GNG_Scotty_MotifXcorrPosLagShuffle_ArraySubmit
%
% Discover all valid sessions under filePathBase (same discovery logic as
% the serial batch runner), build one SLURM array manifest entry per
% session, and submit a single array job. Each array task runs
% motifH_perTrial_crosscorr_posLag_withinTrialTimeShuffle_func for its
% session and SAVES the result into that session's own Matfiles folder
% (rather than returning it to this function).
%
% Required helper on Scotty MATLAB path:
%   MotifXcorrPosLagShuffle_Spock_ArrayTask.m
%
% USAGE
%   Data_Dual_GNG_Scotty_MotifXcorrPosLagShuffle_ArraySubmit( ...
%       filePathBase, fileKeyword, ...
%       'fileKeywordBeh',        '_alignedPupilOrofacial.mat', ...
%       'excludeAnimals',        {'m1893','m1897'}, ...
%       'xcorrLagWindowSec',     1.0, ...
%       'posLagPoolWindowSec',   0.5, ...
%       'doTimeShuffle',         true, ...
%       'nShuffle',              100, ...
%       'shuffleMethod',         'circshift', ...
%       'rngSeed',               1, ...
%       'zscoreHs',              true, ...
%       'useSymmetry',           true, ...
%       'doFisherZ',             true, ...
%       'clipR',                 0.999, ...
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

% Pass-through params for motifH_perTrial_crosscorr_posLag_withinTrialTimeShuffle_func
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
addParameter(p, 'showProgress',        true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'progressEvery',       10,   @(x) isnumeric(x) && isscalar(x));

% SLURM resource params
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

% Bundle the function's own name-value pairs so the array task can just
% do: motifH_perTrial_crosscorr_posLag_withinTrialTimeShuffle_func(filePath, fileKeyword, funcNV{:})
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
    'showProgress',        p.Results.showProgress, ...
    'progressEvery',       p.Results.progressEvery};

fprintf('\n============================================================\n');
fprintf('Motif xcorr pos-lag time-shuffle ARRAY submission\n');
fprintf('filePathBase: %s\n', filePathBase);
fprintf('fileKeyword:  %s\n', fileKeyword);
fprintf('nShuffle:     %d\n', p.Results.nShuffle);
fprintf('sbatch_time/memory/cpus: %d / %d / %d\n', sbatch_time, sbatch_memory, sbatch_cpus);
fprintf('============================================================\n');

%% -------------------- Discover sessions (mirrors serial batch runner) --------------------
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

% All array tasks share the same fileKeyword + funcNV; store once.
for s = 1:nSessions
    filePath_bucket{s}     = ConvertWinToBucketPath(sessionRows(s).filePath);
    filePath_mat_bucket{s} = ConvertWinToBucketPath(sessionRows(s).filePath_mat);
    mIdC{s}    = sessionRows(s).mId;
    headerC{s} = sessionRows(s).header;

    saveName = sprintf('%s_xcorrPosLagShuffle_n%d_%s.mat', ...
        sessionRows(s).header, p.Results.nShuffle, dateStr);
    save_fn_bucket{s} = ConvertWinToBucketPath(fullfile(sessionRows(s).filePath_mat, saveName));
end

manifestDir = fullfile(filePathBase, 'collectData', 'xcorr_withinTrialTimeShuffle', 'arrayManifests');
if exist(manifestDir, 'dir') ~= 7
    mkdir(manifestDir);
end

manifestName = sprintf('xcorrPosLagShuffle_arrayManifest_n%d_%s%s', ...
    p.Results.nShuffle, dateStr, fileKeyword);
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
jobStem = sanitize_for_slurm(sprintf('xcorrShuf_n%d_%s_ARRAY', p.Results.nShuffle, dateStr));

script_name = WriteBashScriptWinScotty(jobStem, ...
    'MotifXcorrPosLagShuffle_Spock_ArrayTask', ...
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
scriptFullPath_esc    = sprintf('%s/%s', dynamicScriptsDir, script_name); % backslash-escaped spaces, for use unquoted in remote shell commands

fprintf('\nEnsuring remote out/ directory exists...\n');
respMkdir = ssh2_command_scotty(s_conn, sprintf('mkdir -p %s/out', dynamicScriptsDir));
if isfield(respMkdir, 'command_result') && ~isempty(respMkdir.command_result)
    disp(respMkdir.command_result);
end

%% ---- Read the generated script back, edit locally, write back, verify ----
% Doing the edit in MATLAB (rather than blind remote awk/sed) avoids the
% multi-layer shell-quoting fragility that silently dropped the --array,
% --cpus-per-task, and --output edits in a prior attempt (verified by
% inspecting the resulting sbatch script -- none of the three directives
% had actually changed, meaning the remote patch command never took
% effect, even though it looked syntactically fine on the MATLAB side).
fprintf('\nFetching generated script for local patching:\n%s\n', scriptFullPath_esc);
respCat = ssh2_command_scotty(s_conn, sprintf('cat %s', scriptFullPath_esc));

if ~isfield(respCat, 'command_result') || isempty(respCat.command_result)
    error('Could not read back generated script from Scotty: %s', scriptFullPath_esc);
end

origLines = respCat.command_result(:);
if ~iscell(origLines)
    origLines = cellstr(origLines);
end

% ssh2_command_scotty prepends a banner/MOTD line to every command's
% output (confirmed via diagnostic testing -- e.g. "PLEASE DELETE
% UNUSED/UNNEEDED DATA..." shows up ahead of actual `cat` content on
% every call). That banner is NOT part of the script file; if left in,
% it would land as a literal first line of the rewritten script, ahead
% of the shebang, breaking it. Discard everything before the real
% shebang line rather than trusting command_result to be pure file
% content.
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

% Insert/replace each directive using a shared helper that always checks
% for an existing "#SBATCH --key=" line first and replaces it in place;
% only inserts a new line (right after the first #SBATCH line) if no
% matching directive exists yet. This matters because
% WriteBashScriptWinScotty's default template already emits its own
% "#SBATCH --cpus-per-task=1" line -- blindly inserting a second
% --cpus-per-task line alongside it would leave two conflicting
% directives in the file, and whichever one SLURM resolves as "last wins"
% could silently reintroduce the exact single-core bug this patch exists
% to fix. --array has no default line in the observed template, so for
% it this reduces to a plain insert, but running it through the same
% duplicate-safe path costs nothing and protects against a future
% template change adding one.
newLines = origLines;
newLines = upsertSbatchDirective(newLines, '#SBATCH --array=',           arrayDirective);
newLines = upsertSbatchDirective(newLines, '#SBATCH --cpus-per-task=',   cpuDirective);
newLines = upsertSbatchDirective(newLines, '#SBATCH --output=',          outputDirective);

patchedContent = strjoin(newLines, newline);
if ~endsWith(patchedContent, newline)
    patchedContent = [char(patchedContent), newline]; % concatenation, NOT '+' (which does numeric char-code addition)
end

% Write back via base64 encode/decode round-trip rather than a heredoc.
% Diagnostic testing showed multi-line commands (heredocs) sent through
% this SSH wrapper on Windows do not reliably preserve embedded newlines
% ("here-document ... delimited by end-of-file" error) -- the remote side
% never saw them as separate lines. Base64 avoids the problem entirely by
% sending the whole payload as a single line containing only
% [A-Za-z0-9+/=], which a prior round-trip test confirmed survives intact
% (including literal $VARS and embedded quotes) end to end.
b64Content = matlab.net.base64encode(uint8(char(patchedContent)));
writeCmd = sprintf('printf ''%%s'' ''%s'' | base64 -d > %s', b64Content, scriptFullPath_esc);

fprintf('Writing patched script back to Scotty (base64 round-trip)...\n');
respWrite = ssh2_command_scotty(s_conn, writeCmd);
if isfield(respWrite, 'command_result') && ~isempty(respWrite.command_result)
    disp(respWrite.command_result);
end

% Verify: re-fetch and assert the three directives actually landed. Fail
% loudly here rather than silently submitting a broken sbatch script.
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

fprintf('\nSubmitted ARRAY job %s for %d sessions.\n', array_job_id, nSessions);
fprintf('Array throttle: %d concurrent tasks max.\n', array_throttle);
fprintf('Array directive: %s\n', arrayDirective);
fprintf('CPU directive:   %s\n', cpuDirective);

%% -------------------- Save local submission record --------------------
swarm_id = cell(nSessions, 1);
for s = 1:nSessions
    swarm_id{s} = sprintf('%s_%d', array_job_id, s);
end

recordName = sprintf('xcorrPosLagShuffle_submissionRecord_n%d_%s.mat', p.Results.nShuffle, dateStr);
save(fullfile(manifestDir, recordName), ...
    'array_job_id', 'array_throttle', 'arrayDirective', 'cpuDirective', ...
    'script_name', 'manifest_fn', 'manifest_fn_bucket', ...
    'save_fn_bucket', 'filePath_bucket', 'filePath_mat_bucket', ...
    'mIdC', 'headerC', 'nSessions', 'swarm_id', 'funcNV', 'fileKeyword', ...
    'sbatch_time', 'sbatch_memory', 'sbatch_cpus', '-v7.3');

fprintf('\nSaved submission record:\n%s\n', fullfile(manifestDir, recordName));

clearvars s_conn

end

%% ========================================================================
function linesOut = upsertSbatchDirective(linesIn, keyPrefix, fullLine)
% UPSERTSBATCHDIRECTIVE
%   If a line starting with keyPrefix (e.g. '#SBATCH --cpus-per-task=')
%   already exists in linesIn, replace that line in place with fullLine.
%   Otherwise, insert fullLine as a new line immediately after the first
%   '#SBATCH' line. Either way, the result has exactly one line matching
%   keyPrefix -- never zero, never two. linesIn/linesOut are column
%   cellstr arrays of script lines.

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