function Data_GNG_Scotty_GlobalDAMotifPixelXcorr_ArraySubmit(filePathBase, fileKeyword, fileW, varargin)
%DATA_GNG_SCOTTY_GLOBALDAMOTIFPIXELXCORR_ARRAYSUBMIT
%
% Discover all valid sessions under filePathBase, build one SLURM array
% manifest entry per session, and submit a single array job. Each array
% task runs globalDA_motifPixel_xcorr_func for its session and saves the
% result into that session's own Matfiles folder.
%
% Structure follows Data_GNG_Scotty_MotifXcorrRsdPosLagShuff_ArraySubmit
% (same discovery loop, ConvertWinToBucketPath, WriteBashScriptWinScotty +
% base64 patch round-trip, submission record).
%
% DISCOVERY REQUIRES ALL THREE PER-SESSION FILES
%   <header>_green_dff_smCollect.mat    (dffsmCell)
%   <header>_green_tbytDat_dff.mat      (tbytDat with frameT)
%   refit file matching fileKeyword     (hC, tbytDat)
% A session missing any of them is reported and excluded rather than
% consuming an array slot and failing on the node. The basis file (fileW,
% W_basis + nanpxs) is COMMON to every session and is passed once, not
% discovered per session.
%
% WHY frameT: <header>_green_dff_smCollect.mat is written WITHOUT -append
% by dffPostprocess_auditory_gng_dual, so anything appended to it later
% (e.g. DAglobalC) is destroyed whenever the postprocess is re-run. The
% compute function therefore takes DA frame times from frameT in the green
% tbytDat_dff file, which is written in the same pass as dffsmCell.
%
% USAGE
%   Data_GNG_Scotty_GlobalDAMotifPixelXcorr_ArraySubmit( ...
%       filePathBase, ...
%       '_red_dff_combined.mat', ...                 % refit keyword
%       'Z:\...\clusterW_output_..._L10K10_062526.mat', ...
%       'excludeAnimals',   {'m1893','m1897'}, ...
%       'hRow',             1, ...
%       'alignDt',          0.05, ...
%       'maxLagSec',        2.0, ...
%       'pixelSelect',      'footprint', ...
%       'nShuffleTrial',    1000, ...
%       'doCircShuffle',    true, ...
%       'nShuffleCirc',     200, ...
%       'sbatch_time',      480, ...
%       'sbatch_memory',    32, ...
%       'sbatch_cpus',      12, ...
%       'array_throttle',   8);

%% -------------------- Parse inputs --------------------
p = inputParser;
addRequired(p, 'filePathBase', @(x) ischar(x) || isstring(x));
addRequired(p, 'fileKeyword',  @(x) ischar(x) || isstring(x));
addRequired(p, 'fileW',        @(x) ischar(x) || isstring(x));

addParameter(p, 'fileKeywordDA',     '_green_dff_smCollect.mat', @(x) ischar(x) || isstring(x));
addParameter(p, 'fileKeywordTbytDA', '_green_tbytDat_dff.mat',   @(x) ischar(x) || isstring(x));
addParameter(p, 'excludeAnimals',    {'m1893','m1897'}, @iscell);

% forwarded to globalDA_motifPixel_xcorr_func
addParameter(p, 'hRow',              1,    @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'tRow',              3,    @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'alignWin',          [-0.9 5], @(x) isnumeric(x) && numel(x)==2);
addParameter(p, 'alignDt',           0.05, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'maxLagSec',         2.0,  @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'pixelSelect',       'footprint', @(s) ischar(s) || isstring(s));
addParameter(p, 'topFrac',           0.20, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'doPSTHSubtraction', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'doFisherZ',         true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'clipR',             0.999, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'minCorrectTrials',  10,   @(x) isnumeric(x) && isscalar(x) && x>=1);
addParameter(p, 'nShuffleTrial',     1000, @(x) isnumeric(x) && isscalar(x) && x>=1);
addParameter(p, 'doCircShuffle',     true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'nShuffleCirc',      200,  @(x) isnumeric(x) && isscalar(x) && x>=1);
addParameter(p, 'circMotifs',        [],   @(x) isempty(x) || isnumeric(x));
addParameter(p, 'rngSeed',           1,    @(x) isempty(x) || (isnumeric(x) && isscalar(x)));

addParameter(p, 'sbatch_time',    480, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'sbatch_memory',  32,  @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'sbatch_cpus',    12,  @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'array_throttle', 8,   @(x) isnumeric(x) && isscalar(x) && x>0);

parse(p, filePathBase, fileKeyword, fileW, varargin{:});

filePathBase      = char(p.Results.filePathBase);
fileKeyword       = char(p.Results.fileKeyword);
fileW             = char(p.Results.fileW);
fileKeywordDA     = char(p.Results.fileKeywordDA);
fileKeywordTbytDA = char(p.Results.fileKeywordTbytDA);
excludeAnimals    = p.Results.excludeAnimals;

sbatch_time    = p.Results.sbatch_time;
sbatch_memory  = p.Results.sbatch_memory;
sbatch_cpus    = p.Results.sbatch_cpus;
array_throttle = p.Results.array_throttle;

assert(exist(fileW, 'file') == 2, 'Common basis file not found:\n  %s', fileW);

% Naming splits on doPSTHSubtraction so the two variants can be written
% into the SAME per-session Matfiles folders without colliding. The two
% keywords deliberately share no substring ("Rsd" vs "Raw"), matching the
% convention used for the motif-motif variants -- so any future
% substring-based collector can never cross-match them.
if p.Results.doPSTHSubtraction
    baseKeyword = 'globalDAMotifPixelXcorrRsd';
else
    baseKeyword = 'globalDAMotifPixelXcorrRaw';
end

funcNV = { ...
    'hRow',              p.Results.hRow, ...
    'tRow',              p.Results.tRow, ...
    'alignWin',          p.Results.alignWin, ...
    'alignDt',           p.Results.alignDt, ...
    'maxLagSec',         p.Results.maxLagSec, ...
    'pixelSelect',       char(string(p.Results.pixelSelect)), ...
    'topFrac',           p.Results.topFrac, ...
    'doPSTHSubtraction', p.Results.doPSTHSubtraction, ...
    'doFisherZ',         p.Results.doFisherZ, ...
    'clipR',             p.Results.clipR, ...
    'minCorrectTrials',  p.Results.minCorrectTrials, ...
    'nShuffleTrial',     p.Results.nShuffleTrial, ...
    'doCircShuffle',     p.Results.doCircShuffle, ...
    'nShuffleCirc',      p.Results.nShuffleCirc, ...
    'circMotifs',        p.Results.circMotifs, ...
    'rngSeed',           p.Results.rngSeed, ...
    'saveKeyword',       baseKeyword, ...
    'useParfor',         true, ...
    'verbose',           true};

fprintf('\n============================================================\n');
fprintf('Global DA x motif PIXEL-SPACE xcorr ARRAY submission\n');
fprintf('filePathBase : %s\n', filePathBase);
fprintf('refit keyword: %s\n', fileKeyword);
fprintf('basis (W)    : %s\n', fileW);
fprintf('pixelSelect  : %s | alignDt %g s | maxLag %g s\n', ...
    char(string(p.Results.pixelSelect)), p.Results.alignDt, p.Results.maxLagSec);
fprintf('nShuffleTrial: %d | circShuffle %d (n=%d)\n', ...
    p.Results.nShuffleTrial, p.Results.doCircShuffle, p.Results.nShuffleCirc);
fprintf('doPSTHSubtraction: %d  (baseKeyword: %s)\n', p.Results.doPSTHSubtraction, baseKeyword);
fprintf('sbatch time/mem/cpus: %d / %d / %d\n', sbatch_time, sbatch_memory, sbatch_cpus);
fprintf('============================================================\n');

%% -------------------- Discover sessions --------------------
mListC = GrabFiles_sort_trials('m*jRGECO', 0, {filePathBase});
mListI = cellfun(@(a) all(~contains(a, excludeAnimals)), mListC);
mListC = mListC(mListI);

fprintf('\nFound %d animal folders after exclusion.\n', numel(mListC));

sessionRows = struct('mId', {}, 'header', {}, 'filePath', {}, 'filePath_mat', {}, ...
    'fileDA', {}, 'fileTbytDA', {}, 'fileH', {});
skipRows = struct('mId', {}, 'header', {}, 'reason', {});

for j = 1:numel(mListC)
    mId = cell2mat(regexp(mListC{j}, 'm\d{4}', 'match'));
    filePathSessions = find_keyword_containing_folder(mListC{j}, mId, 'recursive', false);

    for jj = 1:numel(filePathSessions)
        filePathC = find_keyword_containing_folder(filePathSessions{jj}, 'task', 'recursive', false);
        if isempty(filePathC); continue; end
        filePath = filePathC{1};

        filePath_mat = cell2mat(find_keyword_containing_folder(filePath, 'Matfiles', 'recursive', false));
        if isempty(filePath_mat); continue; end

        header = extract_date_animalID_header(filePath);

        % All three per-session inputs must exist. Recording WHICH one is
        % missing is the point -- a bare "skipped" would leave you grepping
        % the tree by hand later.
        fDA = find_keyword_containing_files(filePath_mat, fileKeywordDA, 'recursive', false);
        if isempty(fDA)
            skipRows(end+1) = struct('mId', mId, 'header', header, 'reason', ...
                sprintf('no %s', fileKeywordDA)); %#ok<AGROW>
            continue;
        end

        fTB = find_keyword_containing_files(filePath_mat, fileKeywordTbytDA, 'recursive', false);
        if isempty(fTB)
            skipRows(end+1) = struct('mId', mId, 'header', header, 'reason', ...
                sprintf('no %s (run dffPostprocess first)', fileKeywordTbytDA)); %#ok<AGROW>
            continue;
        end

        fH = find_keyword_containing_files(filePath_mat, fileKeyword, 'recursive', true);
        if isempty(fH)
            skipRows(end+1) = struct('mId', mId, 'header', header, 'reason', ...
                sprintf('no refit file matching %s', fileKeyword)); %#ok<AGROW>
            continue;
        end

        sessionRows(end+1) = struct( ...
            'mId', mId, 'header', header, ...
            'filePath', filePath, 'filePath_mat', filePath_mat, ...
            'fileDA', fDA{1}, 'fileTbytDA', fTB{1}, 'fileH', fH{1}); %#ok<AGROW>
    end
end

nSessions = numel(sessionRows);
fprintf('\nDiscovered %d sessions with all three required files.\n', nSessions);
for s = 1:nSessions
    fprintf('  [%d] %s  %s\n', s, sessionRows(s).mId, sessionRows(s).header);
end

if ~isempty(skipRows)
    fprintf('\n%d session(s) EXCLUDED:\n', numel(skipRows));
    for s = 1:numel(skipRows)
        fprintf('  %s  %s  -- %s\n', skipRows(s).mId, skipRows(s).header, skipRows(s).reason);
    end
end

if nSessions == 0
    error('No sessions had all three required files under %s.', filePathBase);
end

%% -------------------- Connect to Scotty --------------------
keyFile = fullfile(getenv('USERPROFILE'), '.ssh', 'id_ed25519_scotty_matlab');
s_conn = ssh2_command_scotty('connect', 'scotty', keyFile);

%% -------------------- Build manifest --------------------
dateStr = string(datetime('today','Format','MMddyy'));

fileDA_bucket     = cell(nSessions, 1);
fileTbytDA_bucket = cell(nSessions, 1);
fileH_bucket      = cell(nSessions, 1);
saveDir_bucket    = cell(nSessions, 1);
mIdC              = cell(nSessions, 1);
headerC           = cell(nSessions, 1);

for s = 1:nSessions
    fileDA_bucket{s}     = ConvertWinToBucketPath(sessionRows(s).fileDA);
    fileTbytDA_bucket{s} = ConvertWinToBucketPath(sessionRows(s).fileTbytDA);
    fileH_bucket{s}      = ConvertWinToBucketPath(sessionRows(s).fileH);
    saveDir_bucket{s}    = ConvertWinToBucketPath(sessionRows(s).filePath_mat);
    mIdC{s}              = sessionRows(s).mId;
    headerC{s}           = sessionRows(s).header;
end

fileW_bucket = ConvertWinToBucketPath(fileW);   % common, single path

manifestDir = fullfile(filePathBase, 'collectData', 'globalDA_motifPixelXcorr', 'arrayManifests');
if exist(manifestDir, 'dir') ~= 7
    mkdir(manifestDir);
end

manifestName = sprintf('%s_arrayManifest_nT%d_%s.mat', ...
    baseKeyword, p.Results.nShuffleTrial, dateStr);
manifest_fn = fullfile(manifestDir, manifestName);

% fileW_bucket is stored under the name the array task expects (fileW_bucket),
% as a single char rather than a per-session cell.
save(manifest_fn, ...
    'fileDA_bucket', 'fileTbytDA_bucket', 'fileH_bucket', 'fileW_bucket', ...
    'saveDir_bucket', 'mIdC', 'headerC', 'funcNV', 'nSessions', ...
    'fileKeyword', 'fileKeywordDA', 'fileKeywordTbytDA', '-v7.3');

fprintf('\nSaved array manifest:\n%s\n', manifest_fn);
manifest_fn_bucket = ConvertWinToBucketPath(manifest_fn);

%% -------------------- Write + patch sbatch array script --------------------
if p.Results.doPSTHSubtraction, jobTag = 'DAmotifPixRsd'; else, jobTag = 'DAmotifPixRaw'; end
jobStem = sanitize_for_slurm(sprintf('%s_nT%d_%s_ARRAY', jobTag, p.Results.nShuffleTrial, dateStr));

script_name = WriteBashScriptWinScotty(jobStem, ...
    'GlobalDAMotifPixelXcorr_Spock_ArrayTask', ...
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

dynamicScriptsDir  = '/jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts';
scriptFullPath_esc = sprintf('%s/%s', dynamicScriptsDir, script_name);

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
if ~iscell(origLines), origLines = cellstr(origLines); end

shebangIdx = find(startsWith(strtrim(origLines), '#!'), 1, 'first');
if isempty(shebangIdx)
    error(['Could not locate a shebang ("#!...") line in the fetched script -- ' ...
           'cannot safely distinguish real content from SSH banner output. ' ...
           'Raw response was:\n%s'], strjoin(origLines, newline));
end
if shebangIdx > 1
    fprintf('Discarding %d banner/MOTD line(s) before the script content.\n', shebangIdx - 1);
end
origLines = origLines(shebangIdx:end);

newLines = origLines;
newLines = upsertSbatchDirective(newLines, '#SBATCH --array=',         arrayDirective);
newLines = upsertSbatchDirective(newLines, '#SBATCH --cpus-per-task=', cpuDirective);
newLines = upsertSbatchDirective(newLines, '#SBATCH --output=',        outputDirective);

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
if ~iscell(verifyLines), verifyLines = cellstr(verifyLines); end

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
fprintf('Verified patched script contains --array, --cpus-per-task=%d, and --output=...%%A_%%a.out.\n', sbatch_cpus);

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

recordName = sprintf('%s_submissionRecord_nT%d_%s.mat', baseKeyword, p.Results.nShuffleTrial, dateStr);
save(fullfile(manifestDir, recordName), ...
    'array_job_id', 'array_throttle', 'arrayDirective', 'cpuDirective', ...
    'script_name', 'manifest_fn', 'manifest_fn_bucket', ...
    'fileDA_bucket', 'fileTbytDA_bucket', 'fileH_bucket', 'fileW_bucket', ...
    'saveDir_bucket', 'mIdC', 'headerC', 'nSessions', 'swarm_id', ...
    'funcNV', 'fileKeyword', 'baseKeyword', 'skipRows', 'fileW', ...
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
