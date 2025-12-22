function [XcorrMat, lags] = computeMotifXcorr(HsY3, maxLag)
% COMPUTEMOTIFXCORR  Motif–motif cross-correlograms for one session.
%
%   [XcorrMat, lags] = computeMotifXcorr(HsY3, maxLag)
%
% INPUTS
%   HsY3   : [N x K x T] array
%            N = # trials, K = # motifs, T = # time bins.
%            Typically this is Hs.Y3 (often already z-scored).
%
%   maxLag : (optional) nonnegative integer, number of time bins
%            to include on each side of zero lag.
%            Default = 10 (i.e., lags = -10:10).
%
% OUTPUTS
%   XcorrMat : [K x K x L] array of cross-correlations,
%              where L = 2*maxLag + 1.
%              XcorrMat(i,j,:) is the cross-correlation between
%              motif i and motif j across all trials/time, using
%              MATLAB's xcorr(...,'coeff') normalization.
%
%   lags     : [1 x L] vector of integer lags in bins:
%              lags = -maxLag : maxLag
%
% NOTES
%   • Positive lags (lags > 0): motif j tends to FOLLOW motif i.
%   • Negative lags (lags < 0): motif j tends to PRECEDE motif i.
%   • Non-finite values in HsY3 (NaN/Inf) are set to 0 before
%     computing cross-correlations.
%
% EXAMPLE
%   [Xc, lags] = computeMotifXcorr(Hs.Y3, 15);
%   imagesc(lags, 1:K, squeeze(Xc(4,:,:))); axis xy;
%   xlabel('Lag (bins)'); ylabel('Target motif');
%   title('Cross-corr of motif 4 with all others');
%

    % ---- defaults ----
    if nargin < 2 || isempty(maxLag)
        maxLag = 10;
    end
    if maxLag < 0 || maxLag ~= round(maxLag)
        error('maxLag must be a nonnegative integer.');
    end

    % ---- basic sizes ----
    [N, K, T] = size(HsY3); %#ok<NASGU>  % N unused but might be useful to inspect
    L = 2 * maxLag + 1;
    lags = -maxLag:maxLag;

    % ---- preallocate ----
    XcorrMat = zeros(K, K, L);

    % ---- compute cross-correlations ----
    for i = 1:K
        Xi = reshape(HsY3(:, i, :), [], 1);  % (N*T) x 1
        Xi(~isfinite(Xi)) = 0;

        for j = 1:K
            Xj = reshape(HsY3(:, j, :), [], 1);  % (N*T) x 1
            Xj(~isfinite(Xj)) = 0;

            % normalized cross-correlation (coeff = correlation coef)
            XcorrMat(i, j, :) = xcorr(Xi, Xj, maxLag, 'coeff');
        end
    end
end
