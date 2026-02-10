function S = motifH_trialTime_pca_func(filePath, fileKeyword, varargin)
% MOTIFH_TRIALTIME_PCA_FUNC
%   Per-session helper that loads aligned motif H (tbytDat_hAligned),
%   stacks trials (Hs.Y3: N x K x T), and converts it into a sample×motif
%   matrix X where samples are (trial,timebin) rows.
%
%   This function is designed for a 2-pass workflow:
%     PASS 1: call this per session -> store S (contains X + bookkeeping)
%     PASS 2: build global/expert PCA bases across sessions, then project.
%
% OUTPUT (struct S):
%   S.meta.header, S.meta.mouseId, S.meta.filePath
%   S.trI (trial type masks)
%   S.H (Hs info)
%   S.X (sample x K), S.rowInfo (trialIdx, timeIdx, timeSec)
%   S.pca.withinSession (optional): coeff, score, latent, explained, mu
%
% Name–Value
%   'doZscore'          : z-score H inside stack_trials_H (default true)
%   'useTrials'         : 'all'|'go'|'nogo' (default 'all')
%   'timeWindowSec'     : [] or [t0 t1] in seconds relative to Hs.winCtrs (default [])
%   'doWithinSessPCA'   : true/false (default false)
%   'nPC'               : # PCs to keep for within-session PCA (default 10)
%   'pcaAlgorithm'      : 'svd' (default) or 'eig' (MATLAB pca uses svd by default)
%   'signFix'           : 'none'|'maxabs' (default 'maxabs') for within-session PCA
%   'nanPolicy'         : 'dropRows' (default) | 'zeroFill'
%
% NOTE:
%   Global/expert-anchored PCA projections are added in PASS 2 by a separate function.

% -------------------- parse --------------------
p = inputParser;
p.addParameter('doZscore', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('useTrials', 'all', @(s) any(strcmpi(string(s), ["all","go","nogo"])));
p.addParameter('timeWindowSec', [], @(x) isempty(x) || (isnumeric(x)&&numel(x)==2));
p.addParameter('doWithinSessPCA', false, @(x)islogical(x)&&isscalar(x));
p.addParameter('nPC', 10, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('pcaAlgorithm', 'svd', @(s)ischar(s)||isstring(s));
p.addParameter('signFix', 'maxabs', @(s) any(strcmpi(string(s), ["none","maxabs"])));
p.addParameter('nanPolicy', 'dropRows', @(s) any(strcmpi(string(s), ["dropRows","zeroFill"])));
p.parse(varargin{:});
opt = p.Results;

% -------------------- 0) Grab files --------------------
header      = extract_date_animalID_header(filePath);
keyword_beh = '_alignedPupilOrofacial.mat';

filePath_matfiles = cell2mat(GrabFiles_sort_trials('Matfiles', 0, {filePath}));
filePath_H        = cell2mat(GrabFiles_sort_trials([header '*' fileKeyword], 0, {filePath_matfiles}));
filePath_B        = cell2mat(GrabFiles_sort_trials([header '*' keyword_beh], 0, {filePath_matfiles}));

load(filePath_H, 'tbytDat_hAligned');
load(filePath_B, 'tbytDat');

trI = trialTypeInfoAuditoryGngTbytDat(tbytDat);
assert(size(tbytDat_hAligned, 2) == numel(tbytDat), 'Mismatch in trial counts.');

% infer mouseId
tok = regexp(string(header), '(m\d{3,5})', 'tokens', 'once');
if isempty(tok), mouseId = "unknown";
else, mouseId = string(tok{1});
end

% -------------------- 1) Stack trials --------------------
Hs = stack_trials_H(tbytDat_hAligned, 'zscore', opt.doZscore); % Hs.Y3: N x K x T
[N, K, T] = size(Hs.Y3);
tvec = Hs.winCtrs(:)'; % 1 x T

% time selection
tMask = true(1,T);
if ~isempty(opt.timeWindowSec)
    tMask = (tvec >= opt.timeWindowSec(1)) & (tvec <= opt.timeWindowSec(2));
end
tSel = find(tMask);
Tsel = numel(tSel);

% trial selection
useTrials = lower(string(opt.useTrials));
switch useTrials
    case "go"
        trMask = trI.goI(:);
    case "nogo"
        trMask = trI.nogoI(:);
    otherwise
        trMask = true(N,1);
end
trSel = find(trMask);
Nsel = numel(trSel);

% extract selected cube: [Nsel x K x Tsel]
Y = Hs.Y3(trSel, :, tSel);

% -------------------- 2) Build sample×motif matrix X --------------------
% samples = trial-time rows
% Xraw is TIME-MAJOR (time blocks, trials inside each block)
Xraw = reshape(permute(Y, [1 3 2]), [], K); % [Nsel*Tsel x K]

% row bookkeeping MUST match TIME-MAJOR order
trialIdx = repmat(trSel(:), Tsel, 1);      % [Nsel*Tsel x 1] trials repeat per time bin
timeIdx  = repelem(tSel(:), Nsel, 1);      % [Nsel*Tsel x 1] each time bin repeats Nsel times
timeSec  = tvec(timeIdx).';

rowInfoFull = table(trialIdx, timeIdx, timeSec, ...
    'VariableNames', {'trialIdx','timeIdx','timeSec'});

% missingness mask on the ORIGINAL values (what you want to restore later)
wasMissing = any(~isfinite(Xraw), 2);

switch lower(string(opt.nanPolicy))
    case "zerofill"
        X = Xraw;
        X(~isfinite(X)) = 0;          % only for PCA/projection stability
        goodRow = true(size(X,1),1);  % keep rectangular grid
        rowInfo = rowInfoFull;

    otherwise % "dropRows"
        goodRow = all(isfinite(Xraw), 2);
        X = Xraw(goodRow,:);
        rowInfo = rowInfoFull(goodRow,:);
        wasMissing = wasMissing(goodRow);
end

% -------------------- 3) Package outputs --------------------
S = struct();
S.meta = struct('header', header, 'mouseId', char(mouseId), 'filePath', filePath);
S.params = opt;
S.trI = trI;

S.H = struct();
S.H.winCtrs   = Hs.winCtrs;
S.H.winBounds = Hs.winBounds;
S.H.stepSec   = Hs.params.Step;
S.H.K         = K;

S.X = X;                 % sample x K
S.rowInfo = rowInfo;     % maps rows back to (trial,time)

% also keep the selected indices (useful later)
S.select = struct('trialSel', trSel, 'timeSel', tSel, 'useTrials', char(useTrials));

% -------------------- store --------------------
S.Xraw = Xraw;                % (optional but useful for debugging)
S.X    = X;                   % what is used for PCA/projection
S.rowInfo = rowInfo;
S.rowMask = struct();
S.rowMask.goodRow    = goodRow;
S.rowMask.wasMissing = wasMissing;   % <--- THIS is your “restore NaNs” key
S.grid = struct('Nsel', Nsel, 'Tsel', Tsel, 'K', K, 'tSel', tSel, 'trSel', trSel);

% -------------------- 4) Optional within-session PCA --------------------
S.pca = struct();
S.pca.withinSession = [];

if opt.doWithinSessPCA
    nPC = min([opt.nPC, size(X,2), size(X,1)-1]);
    [coeff, score, latent, ~, explained, mu] = pca(X, ...
        'Algorithm', char(opt.pcaAlgorithm), ...
        'NumComponents', nPC);

    % sign convention (important for session-to-session comparability)
    coeff = local_fixPCSigns(coeff, char(opt.signFix));

    S.pca.withinSession = struct( ...
        'coeff', coeff, ...
        'score', score, ...
        'latent', latent, ...
        'explained', explained, ...
        'mu', mu, ...
        'nPC', nPC);
end

end

% ========================= helper =========================
function coeff = local_fixPCSigns(coeff, mode)
% Fix arbitrary sign of PCA components for reproducibility.
% 'maxabs': flip so the largest-magnitude loading is positive.
mode = lower(string(mode));
if mode == "none", return; end

for pc = 1:size(coeff,2)
    v = coeff(:,pc);
    [~, idx] = max(abs(v));
    if v(idx) < 0
        coeff(:,pc) = -coeff(:,pc);
    end
end
end
