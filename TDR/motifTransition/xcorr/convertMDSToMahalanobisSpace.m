function [Yw, ellW] = convertMDSToMahalanobisSpace(Y, ell, dims)
%CONVERTMDSTOMAHALANOBISSPACE
%   Whitening transform that re-expresses MDS coordinates in units of
%   Mahalanobis distance from a reference ellipsoid (e.g. the
%   fast-learner late-session ellipsoid, ellipsoid_combinedCorrect).
%
%   [Yw, ellW] = convertMDSToMahalanobisSpace(Y, ell, dims)
%
%   AFTER THIS TRANSFORM:
%     - Euclidean distance between any two points in Yw EQUALS their
%       Mahalanobis distance (w.r.t. ell.Sigma) in the original space.
%     - The reference ellipsoid becomes a PERFECT SPHERE of radius
%       sqrt(chi2inv(ell.confLevel, 3)), centered at the origin --
%       returned as ellW (mu=[0 0 0], Sigma=eye(3)).
%     - The three output axes are the PRINCIPAL DIRECTIONS of ell.Sigma
%       (sorted by descending eigenvalue -- axis 1 is the direction of
%       greatest fast-learner late-session variability), NOT the
%       original MDS1/2/3 axes. Expect a visually different orientation
%       from the raw-MDS plot -- this is expected, not a bug.
%
% INPUTS
%   Y    : [S x dimFull] MDS coordinates (dimFull may exceed 3, e.g. if
%          mds_dim was set to 4 -- this function slices down to the 3
%          requested dims internally).
%   ell  : struct with .mu (1 x dimFull), .Sigma (dimFull x dimFull),
%          .confLevel -- e.g. ellipsoid_combinedCorrect from
%          computeFastLearnerEllipsoid.
%   dims : which 3 columns of Y / ell to use, e.g. [1 2 3].
%
% OUTPUTS
%   Yw   : [S x 3] whitened coordinates, in Mahalanobis-distance units.
%   ellW : struct with .mu=[0 0 0], .Sigma=eye(3), .confLevel (unchanged
%          from input), .dim=3 -- feed this directly into
%          plotMDSWithEllipsoid (or plotMDSConvergenceAndAmongDistance,
%          if you want convergence distances reported in Mahalanobis
%          units too -- note distToRef there would then just equal
%          ordinary Euclidean norm in Yw, since ellW.Sigma=eye already
%          IS the Mahalanobis metric) in place of the original Y/ell.
%
% EXAMPLE
%   [Yw, ellW] = convertMDSToMahalanobisSpace(Y_combinedCorrect, ellipsoid_combinedCorrect, [1 2 3]);
%   h = plotMDSWithEllipsoid(Yw, sessInfo_combinedCorrect, cmap9, ellW, ...
%           'ellAlpha', 0.10, 'ellN', 50, ...
%           'minAlpha', 0.25, 'maxAlpha', 1.0, ...
%           'minDotSize', 24, 'dotScaleMax', 3.0, 'dims', [1 2 3], ...
%           'trialTag', "combinedCorrectTrials_mahalUnits", ...
%           'printFig', false, 'figSaveDir', figSaveDir);
%   % axis labels / view angle will likely need resetting -- see note below.

dims = dims(:)';
if numel(dims) ~= 3
    error('convertMDSToMahalanobisSpace:dimsMustBeThree', ...
        'dims must select exactly 3 columns for a 3D whitened space (got %d).', numel(dims));
end

mu_full  = ell.mu(:)';
Sig_full = ell.Sigma;

if numel(mu_full) < max(dims) || any(size(Sig_full) < max(dims))
    error('ell.mu / ell.Sigma do not have enough dimensions for the requested dims.');
end

mu    = mu_full(dims);
Sigma = Sig_full(dims, dims);
Sigma = (Sigma + Sigma')/2;
Sigma = Sigma + 1e-8*eye(numel(dims));   % same defensive regularization used elsewhere in this pipeline

[V, L] = eig(Sigma);
eigVals = diag(L);

% Sort by descending eigenvalue so axis 1 = direction of greatest
% variability (purely for interpretability; does not affect distances).
[eigVals, ord] = sort(eigVals, 'descend');
V = V(:, ord);
eigVals = max(eigVals, 1e-12);   % guard against ~0/negative eigenvalues from numerical noise

A = V * diag(1./sqrt(eigVals));   % whitening matrix: A*A' = inv(Sigma)

Yw = (Y(:, dims) - mu) * A;

ellW = struct();
ellW.mu        = zeros(1, numel(dims));
ellW.Sigma     = eye(numel(dims));
ellW.confLevel = ell.confLevel;
ellW.dim       = numel(dims);

end