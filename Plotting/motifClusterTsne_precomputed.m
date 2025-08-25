function Y = motifClusterTsne_precomputed(tcorr_mat, perplexity)
% motifClusterTsne_precomputed
% tcorr_mat : [nMotifs x nMotifs] OR [nMotifs x nMotifs x nLags]
%             pairwise (normalized) temporal correlations across lags.
% perplexity: scalar (optional), default 30

if nargin < 2 || isempty(perplexity), perplexity = 30; end

% Derive max-over-lag correlation matrix C in [-1, 1]
if ndims(tcorr_mat) == 3
    C = max(tcorr_mat, [], 3, 'omitnan');
else
    C = tcorr_mat;
end

% Clip to valid range and handle NaNs (treat unknown as zero corr)
C = max(min(C, 1), -1);
C(isnan(C)) = 0;

% Convert similarity to distance.
% Simple and effective for t-SNE:
%   D = 1 - C  (so perfectly correlated -> 0, anti-correlated -> 2)
D = 1 - C;

% Enforce symmetry and zero diagonal
D = (D + D.')/2;
n = size(D,1);
D(1:n+1:end) = 0;

% Embed distances into an intermediate Euclidean space
[Xmds,~] = cmdscale(D, 50);   % 50 dims is a good start (use mdscale if you prefer)

% Then run MATLAB t‑SNE on those features
perp = min(perplexity, floor((n-1)/3));    % guardrail for perplexity
rng(100);
Y = tsne(Xmds, 'Perplexity', perp, 'NumDimensions', 2);

end
