function ell = computeGroupCentroidEllipsoid(Y, sessInfo, memberIdC, varargin)
%COMPUTEGROUPCENTROIDELLIPSOID
%   Builds a reference centroid (mu) and covariance (Sigma) from a given
%   group's LAST N sessions, pooled across all animals in that group.
%   Identical math to computeFastLearnerEllipsoid (the local function
%   embedded in the six-stream script) -- just extracted as its own
%   standalone, generically-named file so it can be called with ANY
%   reference group (e.g. slow learners, for a control analysis showing
%   that convergence is specific to the fast-learner reference and not a
%   generic drift toward any fixed point), not only from within that one
%   script.
%
%   ell = computeGroupCentroidEllipsoid(Y, sessInfo, memberIdC, ...)
%
% INPUTS
%   Y         : [S x dim] MDS coordinates for ONE stream.
%   sessInfo  : matching sessInfo table (needs mouseId, sessWithin).
%   memberIdC : cellstr of animal IDs defining the reference group (e.g.
%               groupDefs.slow, to build a SLOW-learner reference).
%
% NAME-VALUE ARGS
%   'nFinalSess' : how many of each member's final sessions to pool.
%                  Default: 3 (matches this project's convention).
%   'confLevel'  : stored for bookkeeping/downstream ellipsoid-drawing
%                  use (e.g. plotMDSWithEllipsoid's chi2inv scaling).
%                  Default: 0.95.
%   'robustCov'  : use robustcov() instead of plain cov(). Default: false.
%
% OUTPUT (ell)
%   Same structure as computeFastLearnerEllipsoid's output: .mu, .Sigma,
%   .idxPts, .confLevel, .nFinalSess, .dim, .X -- directly compatible
%   with plotMDSWithEllipsoid, plotMDSConvergenceAndAmongDistance,
%   runConvergenceRMANOVA_fastSlow, etc. anywhere else in this project
%   that expects an "ell" struct.

p = inputParser;
p.addParameter('nFinalSess', 3, @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('confLevel', 0.95, @(x) isnumeric(x) && isscalar(x) && x>0 && x<1);
p.addParameter('robustCov', false, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;

assert(all(ismember({'mouseId','sessWithin'}, sessInfo.Properties.VariableNames)), ...
    'sessInfo must have mouseId and sessWithin columns.');

mouseId    = string(sessInfo.mouseId);
sessWithin = sessInfo.sessWithin;

dim = size(Y, 2);
idxPts = [];

for k = 1:numel(memberIdC)
    iMouse = find(mouseId == string(memberIdC{k}));
    if isempty(iMouse)
        warning('computeGroupCentroidEllipsoid:MouseNotFound', ...
            'Animal "%s" not found in sessInfo.mouseId -- skipping.', memberIdC{k});
        continue;
    end

    [~, ord] = sort(sessWithin(iMouse), 'ascend');
    iMouse = iMouse(ord);

    take = max(1, numel(iMouse)-opt.nFinalSess+1) : numel(iMouse);
    idxPts = [idxPts; iMouse(take)]; %#ok<AGROW>
end

idxPts = unique(idxPts, 'stable');
assert(~isempty(idxPts), 'No valid sessions found for any animal in memberIdC.');

X = Y(idxPts, :);
mu = mean(X, 1);

if opt.robustCov
    Sigma = robustcov(X);
else
    Sigma = cov(X);
end

epsReg = 1e-8;
Sigma = Sigma + epsReg * eye(dim);

ell = struct();
ell.idxPts     = idxPts;
ell.mu         = mu;
ell.Sigma      = Sigma;
ell.confLevel  = opt.confLevel;
ell.nFinalSess = opt.nFinalSess;
ell.dim        = dim;
ell.X          = X;
ell.memberIdC  = memberIdC;   % which group this reference was actually built from -- kept for provenance/sanity checks

end