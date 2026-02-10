function pcaRezC = applyAcrossSessionPCABases(pcaRezC, hdrC_all, hdrC, varargin)
% APPLYACROSSSESSIONPCABASES
%   PASS 2: Build shared PCA bases and project every session onto them.
%
%   Builds:
%     (1) global PCA basis from ALL sessions concatenated
%     (2) per-mouse expert PCA basis from that mouse's expert session (hdrC_all)
%     (3) expert-global PCA basis from expert sessions across selected mice (hdrC)
%
%   Then projects each session's trial_time samples (rows) onto each basis.
%   If Pass1 used nanPolicy='zeroFill', we restore NaNs in PC scores using
%   S.rowMask.wasMissing, and also provide trial×time×PC trajectories.
%
% Inputs
%   pcaRezC  : animal x session cell array; each entry is struct S from motifH_trialTime_pca_func
%   hdrC_all : cellstr, expert header per qualifying mouse (one per mouse)
%   hdrC     : cellstr, expert headers for selected fast-learners (subset of hdrC_all)
%
% Name–Value
%   'nPC'     : # PCs for shared bases (default 10)
%   'signFix' : 'maxabs' (default) | 'none'  (flip PC signs for reproducibility)
%
% Output
%   pcaRezC updated; for each session S:
%     S.pca.globalBasis,            S.pca.globalScore_rows,            S.pca.globalScore_traj
%     S.pca.expertGlobalBasis,      S.pca.expertGlobalScore_rows,      S.pca.expertGlobalScore_traj
%     S.pca.expertWithinMouseBasis, S.pca.expertWithinMouseScore_rows, S.pca.expertWithinMouseScore_traj

% -------------------- parse --------------------
p = inputParser;
p.addParameter('nPC', 10, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('signFix', 'maxabs', @(s) any(strcmpi(string(s),["maxabs","none"])));
p.parse(varargin{:});
opt = p.Results;

nPC = opt.nPC;
signFix = char(opt.signFix);

% -------------------- 1) GLOBAL basis (all sessions) --------------------
Xall = {};
sessRef = {};

for a = 1:size(pcaRezC,1)
    for s = 1:size(pcaRezC,2)
        S = pcaRezC{a,s};

        if isempty(S) || ~isstruct(S), continue; end
        if ~isfield(S,'X') || isempty(S.X), continue; end
        if ~isfield(S,'meta') || ~isfield(S.meta,'header'), continue; end

        Xtrain = local_getTrainRows(S);
        if isempty(Xtrain), continue; end

        Xall{end+1,1}   = Xtrain;        %#ok<AGROW>
        sessRef{end+1,1}= S.meta.header; %#ok<AGROW>
    end
end

if isempty(Xall)
    error('applyAcrossSessionPCABases:NoData', 'No valid session X found in pcaRezC.');
end

Xglob = vertcat(Xall{:});
globBasis = local_fitPCABasis(Xglob, nPC, signFix);

% -------------------- 2) EXPERT-GLOBAL basis (expert sessions across mice) --------------------
XexpG = local_collectByHeaders(pcaRezC, hdrC);
expGlobBasis = local_fitPCABasis(XexpG, nPC, signFix);

% -------------------- 3) EXPERT-within-mouse bases (one expert per mouse) --------------------
expWithinMap = containers.Map('KeyType','char','ValueType','any');

for i = 1:numel(hdrC_all)
    hdr = string(hdrC_all{i});
    if strlength(hdr)==0, continue; end

    [Sexp, ok] = local_findSessionByHeader(pcaRezC, hdr);
    if ~ok || isempty(Sexp) || ~isstruct(Sexp), continue; end
    if ~isfield(Sexp,'X') || isempty(Sexp.X), continue; end

    tok = regexp(hdr, '(m\d{3,5})', 'tokens', 'once');
    if isempty(tok), continue; end
    mid = char(string(tok{1}));

    Xtrain = local_getTrainRows(Sexp);
    if isempty(Xtrain), continue; end

    expWithinMap(mid) = local_fitPCABasis(Xtrain, nPC, signFix);
end

% -------------------- 4) Project every session onto each basis --------------------
for a = 1:size(pcaRezC,1)
    for s = 1:size(pcaRezC,2)

        S = pcaRezC{a,s};
        if isempty(S) || ~isstruct(S), continue; end
        if ~isfield(S,'X') || isempty(S.X), continue; end

        % Ensure PCA field exists
        if ~isfield(S,'pca') || isempty(S.pca), S.pca = struct(); end

        X = S.X;  % rows = (trial,time), cols = motifs

        % Missing-row mask (from original Xraw NaNs before zerofill)
        miss = [];
        if isfield(S,'rowMask') && isfield(S.rowMask,'wasMissing') && ~isempty(S.rowMask.wasMissing)
            miss = logical(S.rowMask.wasMissing(:));
            if numel(miss) ~= size(X,1)
                % If something is inconsistent, ignore restoration rather than crash
                miss = [];
            end
        end

        % ---- global projection ----
        S.pca.globalBasis = globBasis;
        scoreRows = local_project(X, globBasis);
        if ~isempty(miss), scoreRows(miss,:) = NaN; end

        S.pca.globalScore_rows = scoreRows;
        S.pca.globalScore_traj = local_scoresToTrialTimeCube(S, scoreRows);

        % ---- expert-global projection ----
        S.pca.expertGlobalBasis = expGlobBasis;
        scoreRows = local_project(X, expGlobBasis);
        if ~isempty(miss), scoreRows(miss,:) = NaN; end

        S.pca.expertGlobalScore_rows = scoreRows;
        S.pca.expertGlobalScore_traj = local_scoresToTrialTimeCube(S, scoreRows);

        % ---- expert-within-mouse projection ----
        mid = '';
        if isfield(S,'meta') && isfield(S.meta,'mouseId') && ~isempty(S.meta.mouseId)
            mid = S.meta.mouseId;
        else
            % try infer from header if missing
            if isfield(S,'meta') && isfield(S.meta,'header')
                tok = regexp(string(S.meta.header), '(m\d{3,5})', 'tokens', 'once');
                if ~isempty(tok), mid = char(string(tok{1})); end
            end
        end

        if ~isempty(mid) && isKey(expWithinMap, mid)
            b = expWithinMap(mid);

            S.pca.expertWithinMouseBasis = b;
            scoreRows = local_project(X, b);
            if ~isempty(miss), scoreRows(miss,:) = NaN; end

            S.pca.expertWithinMouseScore_rows = scoreRows;
            S.pca.expertWithinMouseScore_traj = local_scoresToTrialTimeCube(S, scoreRows);
        else
            S.pca.expertWithinMouseBasis = [];
            S.pca.expertWithinMouseScore_rows = [];
            S.pca.expertWithinMouseScore_traj = [];
        end

        pcaRezC{a,s} = S;
    end
end

end

% ========================= helpers =========================
function basis = local_fitPCABasis(X, nPC, signFix)
% Fit PCA basis on X (rows=samples, cols=motifs). Uses MATLAB pca().
if isempty(X)
    error('local_fitPCABasis:EmptyX','Empty training matrix.');
end
nPC = min([nPC, size(X,2), size(X,1)-1]);
[coeff, ~, latent, ~, explained, mu] = pca(X, 'NumComponents', nPC);

coeff = local_fixPCSigns(coeff, signFix);

basis = struct( ...
    'coeff', coeff, ...
    'mu', mu, ...
    'latent', latent, ...
    'explained', explained, ...
    'nPC', nPC);
end

function score = local_project(X, basis)
% Project X onto a PCA basis.
Xc = X - basis.mu;
score = Xc * basis.coeff; % [nRow x nPC]
end

function X = local_collectByHeaders(pcaRezC, hdrList)
% Concatenate training rows from sessions matching hdrList.
Xcell = {};
for i = 1:numel(hdrList)
    hdr = string(hdrList{i});
    if strlength(hdr)==0, continue; end

    [S, ok] = local_findSessionByHeader(pcaRezC, hdr);
    if ~ok || isempty(S) || ~isstruct(S), continue; end
    if ~isfield(S,'X') || isempty(S.X), continue; end

    Xtrain = local_getTrainRows(S);
    if isempty(Xtrain), continue; end

    Xcell{end+1,1} = Xtrain; %#ok<AGROW>
end
if isempty(Xcell)
    error('applyAcrossSessionPCABases:NoExpertSessions', ...
        'No expert sessions found in pcaRezC for provided header list.');
end
X = vertcat(Xcell{:});
end

function [S, ok] = local_findSessionByHeader(pcaRezC, header)
% Find the session struct S in pcaRezC whose meta.header matches.
ok = false;
S = [];
header = string(header);

for a = 1:size(pcaRezC,1)
    for s = 1:size(pcaRezC,2)
        Si = pcaRezC{a,s};
        if isempty(Si) || ~isstruct(Si), continue; end
        if ~isfield(Si,'meta') || ~isfield(Si.meta,'header'), continue; end
        if string(Si.meta.header) == header
            S = Si;
            ok = true;
            return;
        end
    end
end
end

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

function score3 = local_scoresToTrialTimeCube(S, scoreRows)
% Map row-wise scores back to [Nsel x Tsel x nPC] cube.
%
% Requires S.grid.Nsel, S.grid.Tsel, and S.rowInfo with trialIdx/timeIdx,
% plus S.select.trialSel and S.select.timeSel.

if ~isfield(S,'grid') || ~isfield(S.grid,'Nsel') || ~isfield(S.grid,'Tsel')
    error('local_scoresToTrialTimeCube:MissingGrid','S.grid missing required fields.');
end
if ~isfield(S,'rowInfo') || isempty(S.rowInfo)
    error('local_scoresToTrialTimeCube:MissingRowInfo','S.rowInfo is missing/empty.');
end
if ~isfield(S,'select') || ~isfield(S.select,'trialSel') || ~isfield(S.select,'timeSel')
    error('local_scoresToTrialTimeCube:MissingSelect','S.select missing required fields.');
end

Nsel = S.grid.Nsel;
Tsel = S.grid.Tsel;
nPC  = size(scoreRows,2);

score3 = nan(Nsel, Tsel, nPC);

% rowInfo stores absolute indices; convert into local indices 1..Nsel / 1..Tsel
[~, trLocal] = ismember(S.rowInfo.trialIdx, S.select.trialSel);
[~, tLocal]  = ismember(S.rowInfo.timeIdx,  S.select.timeSel);

ok = (trLocal>0) & (tLocal>0);

% Defensive size check
if size(scoreRows,1) ~= height(S.rowInfo)
    error('local_scoresToTrialTimeCube:SizeMismatch', ...
        'scoreRows rows (%d) != rowInfo height (%d).', size(scoreRows,1), height(S.rowInfo));
end

for r = find(ok(:))'
    score3(trLocal(r), tLocal(r), :) = scoreRows(r,:);
end
end

function Xtrain = local_getTrainRows(S)
% Training rows exclude originally-missing trial_time samples (if available).
Xtrain = S.X;

if isfield(S,'rowMask') && isfield(S.rowMask,'wasMissing') && ~isempty(S.rowMask.wasMissing)
    keep = ~logical(S.rowMask.wasMissing(:));
    if numel(keep) == size(Xtrain,1)
        Xtrain = Xtrain(keep,:);
    end
end
end
