function [Y, stress, sessInfo, featMat, D, ell] = runSessionMDS_fromXcorr( ...
    xcorrPosLagMatC, mIdC, headerC, dim, fastLearnerIdC, varargin)
%RUNSESSIONMDS_FROMXCORR
% Build session-level feature vectors from directional motif xcorr matrices
% and compute a low-dimensional MDS embedding.
%
% INPUTS
%   xcorrPosLagMatC : {nMice x maxSessions} cell array of KxK positive-lag
%                     motif cross-correlograms (directional).
%   mIdC            : 1 x nMice cell array of mouse folder names or full paths.
%   headerC         : {nMice x maxSessions} cell array of session headers.
%   dim             : target MDS dimension (default = 3).
%
% OUTPUTS
%   Y        : [S x dim] MDS embedding (S = #sessions).
%   stress   : MDS stress value.
%   sessInfo : table with session metadata.
%   featMat  : [S x D] directional features per session (Fisher-z xcorr).
%   D        : [S x S] session distance matrix (1 - corr).
%   ell: struct with fields mu, Sigma, idxPts, confLevel, nFinalSess
% NOTES
%   • Each session's feature vector is constructed by:
%         - Fisher z-transform of 30x30 directional matrix
%         - dropping the diagonal
%         - flattening all K*(K-1) directional entries
%   • Distances = 1 - correlation between session feature vectors.
%
% EXAMPLE
%   [Y, stress, sessInfo] = runSessionMDS_fromXcorr(xcorrPosLagMatC, mIdC, headerC, 3);

if nargin < 4 || isempty(dim), dim = 3; end
if nargin < 5, fastLearnerIdC = {}; end

p = inputParser;
p.addParameter('nFinalSess', 3, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('confLevel', 0.95, @(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);
p.addParameter('robustCov', false, @(x)islogical(x)&&isscalar(x)); % optional
p.parse(varargin{:});
opt = p.Results;

% --- features + MDS (your existing code) ---
[featMat, sessInfo] = buildSessionFeaturesFromXcorr_dir(xcorrPosLagMatC, mIdC, headerC);
featMat_sm = smoothFeatMatByMouse(featMat, sessInfo, 1.0);

D_vec = pdist(featMat_sm, 'correlation');
D     = squareform(D_vec);

[Y, stress] = mdscale(D, dim, 'criterion','metricstress');
fprintf('MDS stress (dim=%d): %.4f\n', dim, stress);

% --- ellipsoid over last-N sessions for fast learners ---
ell = [];
if ~isempty(fastLearnerIdC)
    ell = computeFastLearnerEllipsoid(Y, sessInfo, fastLearnerIdC, ...
        opt.nFinalSess, opt.confLevel, opt.robustCov);
end
end


%% %%%%%%%%%%%%%%%%%%%%% Helper function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [featMat, sessInfo] = buildSessionFeaturesFromXcorr_dir(xcorrPosLagMatC, mListC, headerC)
% Build session-wise feature matrix from directional motif xcorr matrices.
%
% Each session's 30x30 matrix A is interpreted as A(i,j) = j after i
% (positive-lag correlation). We:
%   1) Fisher z-transform A
%   2) Drop the diagonal (i == j)
%   3) Flatten the remaining K*(K-1) entries into a feature vector
%
% OUTPUTS
%   featMat : [S x D] double, S = #sessions, D = K*(K-1)
%   sessInfo: table with session metadata

    nMice   = size(xcorrPosLagMatC,1);
    nSessMx = size(xcorrPosLagMatC,2);

    featMat_list    = {};
    mouseIdx_list   = [];
    mouseId_list    = {};
    sessWithin_list = [];
    header_list     = {};

    for i = 1:nMice
        thisPath = mListC{i};

        % Try to extract something like 'm1045' or 'm1893'
        tok = regexp(thisPath, '(m\d{3,5})', 'tokens', 'once');

        if isempty(tok)
            % Fallback: label by index if regex fails
            mId = sprintf('mouse_%02d', i);
            warning('buildSessionFeaturesFromXcorr_dir:NoMouseId', ...
                'Could not extract mouseId from "%s"; using "%s" instead.', ...
                thisPath, mId);
        else
            mId = tok{1};   % tok is a 1x1 cell containing the match
        end

        sessCount = 0;

        for j = 1:nSessMx
            A = xcorrPosLagMatC{i,j};
            if isempty(A)
                continue;
            end
            sessCount = sessCount + 1;

            % --- Fisher z-transform (clip to avoid Inf) ---
            A = max(min(A, 0.9999), -0.9999);
            Az = atanh(A);                % same size as A

            % --- drop diagonal, keep ALL off-diagonals (directional) ---
            K = size(Az,1);
            mask = ~eye(K);               % false on diag, true elsewhere

            featVec = Az(mask)';          % row vector, length K*(K-1)

            featMat_list{end+1,1} = featVec;      %#ok<AGROW>
            mouseIdx_list         = [mouseIdx_list;   i];         %#ok<AGROW>
            mouseId_list          = [mouseId_list;    {mId}];     %#ok<AGROW>
            sessWithin_list       = [sessWithin_list; sessCount]; %#ok<AGROW>
            header_list           = [header_list;     headerC{i,j}]; %#ok<AGROW>
        end
    end

    featMat = cell2mat(featMat_list);  % S x D

    sessInfo = table( ...
        (1:size(featMat,1))', ...
        mouseIdx_list, ...
        mouseId_list, ...
        sessWithin_list, ...
        header_list, ...
        'VariableNames', {'sessIdx','mouseIdx','mouseId','sessWithin','header'});
end


function featMat_sm = smoothFeatMatByMouse(featMat, sessInfo, sigmaSess)
% sigmaSess: smoothing width in "sessions" (e.g., 1.0)

if nargin < 3, sigmaSess = 1.0; end

featMat_sm = featMat;
mouseU = unique(sessInfo.mouseIdx);

for im = 1:numel(mouseU)
    idx = find(sessInfo.mouseIdx == mouseU(im));
    [~,ord] = sort(sessInfo.sessWithin(idx));
    idx = idx(ord);

    X = featMat(idx,:);                     % nSess x D
    Xs = smoothdata(X, 1, 'gaussian', max(3, 2*ceil(2*sigmaSess)+1));
    featMat_sm(idx,:) = Xs;
end
end

function ell = computeFastLearnerEllipsoid(Y, sessInfo, fastLearnerIdC, nFinalSess, confLevel, robustCov)

dim = size(Y,2);
idxPts = [];

for k = 1:numel(fastLearnerIdC)
    mid = fastLearnerIdC{k};
    iMouse = find(strcmp(sessInfo.mouseId, mid));
    if isempty(iMouse), continue; end

    % sort by within-mouse session index, take last nFinalSess
    [~, ord] = sort(sessInfo.sessWithin(iMouse), 'ascend');
    iMouse = iMouse(ord);

    take = max(1, numel(iMouse)-nFinalSess+1) : numel(iMouse);
    idxPts = [idxPts; iMouse(take)]; %#ok<AGROW>
end

idxPts = unique(idxPts, 'stable');
X = Y(idxPts, :);

mu = mean(X, 1);

if robustCov
    % requires Statistics Toolbox; otherwise set robustCov=false
    Sigma = robustcov(X);
else
    Sigma = cov(X);
end

% guard against singular covariance (common if few points)
epsReg = 1e-8;
Sigma = Sigma + epsReg * eye(dim);

ell = struct();
ell.idxPts     = idxPts;
ell.mu         = mu;
ell.Sigma      = Sigma;
ell.confLevel  = confLevel;
ell.nFinalSess = nFinalSess;
ell.dim        = dim;
ell.X          = X;
end

