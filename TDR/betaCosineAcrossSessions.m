function [S_motif, B_aligned, predNames_master, sessLabels] = ...
    betaCosineAcrossSessions(betaC, XnameC, sessLabels)
% betaCosineAcrossSessions
% Align β across sessions (padding missing predictors with 0) and compute
% cosine similarity between sessions, **per motif**.
%
% Inputs
%   betaC     : 1xS cell, each {P_s x K} β (predictors x motifs)
%   XnameC    : 1xS cell, each {1 x P_s} predictor names (row order of β)
%   sessLabels: 1xS cellstr (optional)
%
% Outputs
%   S_motif          : 1xK cell, each {S x S} cosine similarity across sessions
%   B_aligned        : [P_master x K x S] β tensor (missing predictors padded with 0)
%   predNames_master : {1 x P_master} unified predictor list
%   sessLabels       : passthrough labels
%
% Notes
%   • If some sessions have different #motifs, columns are truncated to the
%     minimum K across sessions (warning issued).
%   • Sessions with all-zero β for a motif yield 0 similarity with others; the
%     diagonal is set to 1.

    % --- 1) Align & pad (reuses our previous function) ---
    [B_aligned, predNames_master, sessLabels] = ...
        betaCosineSimilarityAcrossSessions(betaC, XnameC, sessLabels); %#ok<NASGU>

    [P, K, S] = size(B_aligned);
    S_motif = cell(1, K);

    % --- 2) Cosine per motif ---
    for k = 1:K
        % P x S matrix: columns are sessions’ β for motif k
        Bk = squeeze(B_aligned(:, k, :));      % (P x S)
        % L2-normalize columns (guard against all-zeros)
        colNorm = sqrt(sum(Bk.^2, 1));
        unit = Bk;
        bad = colNorm < eps;                   % sessions with all-zeros (for this motif)
        colNorm(~bad) = 1 ./ colNorm(~bad);
        unit(:, ~bad) = Bk(:, ~bad) .* colNorm(~bad);   % normalize good columns
        unit(:, bad)  = 0;                               % keep zero vector if all-zeros

        Sk = unit.' * unit;                    % (S x S) cosine similarity
        % put 1s on the diagonal (even for all-zero columns)
        Sk(1:S+1:end) = 1;

        % clip numerical strays
        Sk = max(min(Sk, 1), -1);

        S_motif{k} = Sk;
    end
end

%% %%%%%%%%% Helper Function %%%%%%%%%%
function [B_aligned, predNames_master, sessLabels] = betaCosineSimilarityAcrossSessions(betaC, XnameC, sessLabels)
% BETA COSINE PREP — Align & pad β across sessions (no plotting yet)
%
% [B_aligned, predNames_master, sessLabels] = betaCosineSimilarityAcrossSessions(betaC, XnameC, sessLabels)
%
% Inputs
%   betaC     : 1xS cell, each {P_s x K} β-matrix (predictors x motifs)
%   XnameC    : 1xS cell, each {1 x P_s} cellstr of predictor names (order matches rows of β)
%   sessLabels: 1xS cellstr (labels for sessions; kept as-is)
%
% Outputs
%   B_aligned       : [P_master x K x S] β tensor, missing predictors padded with zeros
%   predNames_master: {1 x P_master} unified predictor name list (row order of B_aligned)
%   sessLabels      : passed through (for convenience)
%
% Notes
%   • If some sessions have fewer predictors, those rows are filled with 0.
%   • All sessions must share the same # of motifs K; if not, the function
%     truncates to the minimum K across sessions (with a warning).
%   • This function only prepares the data; downstream cosine-similarity
%     computation/visualization can be done next.

    % ---- basic checks ----
    assert(iscell(betaC)  && isvector(betaC),  'betaC must be a 1xS cell array.');
    assert(iscell(XnameC) && isvector(XnameC), 'XnameC must be a 1xS cell array.');
    S = numel(betaC);
    assert(numel(XnameC)==S, 'betaC and XnameC must have the same length.');
    if nargin < 3 || isempty(sessLabels)
        sessLabels = arrayfun(@(s) sprintf('sess%02d', s), 1:S, 'uni', 0);
    end

    % ---- check motif counts (columns) and harmonize K if needed ----
    K_each = cellfun(@(B) size(B,2), betaC);
    if numel(unique(K_each)) ~= 1
        K_common = min(K_each);
        warning('Motif counts differ across sessions. Truncating all to K = %d (first K columns).', K_common);
        for s = 1:S
            betaC{s} = betaC{s}(:, 1:K_common);
        end
    end
    K = size(betaC{1}, 2);

    % ---- build master predictor list (stable order) ----
    % Concatenate all names and take stable unique.
    allNames = [XnameC{:}];
    assert(iscellstr(allNames), 'Predictor names must be cellstr.');
    predNames_master = unique(allNames, 'stable');
    Pm = numel(predNames_master);

    % ---- allocate aligned tensor and fill ----
    B_aligned = zeros(Pm, K, S);   % zeros act as "missing predictor" padding

    for s = 1:S
        B  = betaC{s};
        nm = XnameC{s};
        % map master -> session indices
        [tf, loc] = ismember(predNames_master, nm);
        if any(tf)
            B_aligned(tf, :, s) = B(loc(tf), :);
        end
        % missing predictors remain zeros
        % quick sanity: protect against row-count mismatch
        if size(B,1) ~= numel(nm)
            error('Session %d: size(beta,1) != numel(X_names).', s);
        end
    end
end
