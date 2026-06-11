function hStats = computeHOverlapStats(h_train, varargin)
% computeHOverlapStats
%
% Quantifies temporal sparsity and cross-motif activation overlap in H.
%
% Input:
%   h_train : K x T matrix
%       K motifs by T frames.
%
% Optional name-value inputs:
%   'smoothWin' : smoothing window length in frames, default = 19
%   'threshold' : threshold for defining active H values, default = []
%                 If empty, uses 1e-6 * max(h_train(:)).
%   'normalizeRows' : whether to normalize each motif's H by its max before
%                     computing overlap, default = false.
%
% Output:
%   hStats : structure with H sparsity and overlap metrics.

% -------------------------------------------------------------------------
% Parse inputs
% -------------------------------------------------------------------------
p = inputParser;
p.addParameter('smoothWin', 19, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('threshold', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('normalizeRows', false, @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});

smoothWin = p.Results.smoothWin;
threshold = p.Results.threshold;
normalizeRows = logical(p.Results.normalizeRows);

% -------------------------------------------------------------------------
% Basic checks
% -------------------------------------------------------------------------
assert(isnumeric(h_train) && ismatrix(h_train), ...
    'h_train must be a numeric K x T matrix.');

H = h_train;
[K, T] = size(H);

% Avoid negative values if tiny numerical negatives exist
H(H < 0) = 0;

% -------------------------------------------------------------------------
% Optional row normalization
% -------------------------------------------------------------------------
if normalizeRows
    rowMax = max(H, [], 2);
    rowMax(rowMax == 0) = 1;
    H = H ./ rowMax;
end

% -------------------------------------------------------------------------
% Define active motifs and active time points
% -------------------------------------------------------------------------
if isempty(threshold)
    threshold = 1e-6 * max(H(:));
end

H_binary = H > threshold;

activeMotif = sum(H, 2) > 0;
K_active = sum(activeMotif);

% -------------------------------------------------------------------------
% Smooth H with a common post hoc kernel
% -------------------------------------------------------------------------
kernel = ones(1, smoothWin);
kernel = kernel ./ sum(kernel);

H_smooth = conv2(H, kernel, 'same');

% -------------------------------------------------------------------------
% Cross-motif H overlap
% -------------------------------------------------------------------------
offDiag = ~eye(K);

offActivation = offDiag * H_smooth;

H_overlap_raw = sum(H .* offActivation, 'all');

H_self = sum(H .* H_smooth, 'all');

H_overlap_total = H_self + H_overlap_raw;

H_overlap_norm = H_overlap_raw ./ (H_overlap_total + eps);

if K > 1
    H_overlap_raw_perOther = H_overlap_raw ./ (K - 1);
else
    H_overlap_raw_perOther = 0;
end

% -------------------------------------------------------------------------
% Sparsity metrics
% -------------------------------------------------------------------------

% Fraction of all H entries that are above threshold
H_density = nnz(H_binary) ./ numel(H_binary);
H_sparsity_fraction_zero = 1 - H_density;

% Hoyer sparsity, computed on vectorized H.
% 0 = dense/equal values, 1 = maximally sparse.
hVec = H(:);
n = numel(hVec);

if norm(hVec, 2) > 0
    H_hoyer_sparsity = (sqrt(n) - norm(hVec, 1) / norm(hVec, 2)) / ...
                       (sqrt(n) - 1 + eps);
else
    H_hoyer_sparsity = NaN;
end

% -------------------------------------------------------------------------
% Active motif count per frame
% -------------------------------------------------------------------------
activeCount = sum(H_binary, 1);

H_active_count_mean = mean(activeCount);
H_active_count_median = median(activeCount);
H_active_count_max = max(activeCount);

H_multi_active_fraction = mean(activeCount > 1);
H_any_active_fraction = mean(activeCount > 0);

% -------------------------------------------------------------------------
% Optional pairwise correlation / similarity of H rows
% -------------------------------------------------------------------------
if K > 1
    H_corr = corr(H');
    H_corr(1:K+1:end) = NaN;
    H_pairwise_corr_mean = mean(H_corr(:), 'omitnan');
    H_pairwise_corr_median = median(H_corr(:), 'omitnan');
else
    H_corr = NaN;
    H_pairwise_corr_mean = NaN;
    H_pairwise_corr_median = NaN;
end

% -------------------------------------------------------------------------
% Store outputs
% -------------------------------------------------------------------------
hStats = struct();

hStats.K = K;
hStats.K_active = K_active;
hStats.T = T;

hStats.smoothWin = smoothWin;
hStats.threshold = threshold;
hStats.normalizeRows = normalizeRows;

hStats.H_overlap_raw = H_overlap_raw;
hStats.H_overlap_norm = H_overlap_norm;
hStats.H_overlap_self = H_self;
hStats.H_overlap_total = H_overlap_total;
hStats.H_overlap_raw_perOther = H_overlap_raw_perOther;

hStats.H_density = H_density;
hStats.H_sparsity_fraction_zero = H_sparsity_fraction_zero;
hStats.H_hoyer_sparsity = H_hoyer_sparsity;

hStats.H_active_count_mean = H_active_count_mean;
hStats.H_active_count_median = H_active_count_median;
hStats.H_active_count_max = H_active_count_max;
hStats.H_multi_active_fraction = H_multi_active_fraction;
hStats.H_any_active_fraction = H_any_active_fraction;

hStats.H_pairwise_corr_mean = H_pairwise_corr_mean;
hStats.H_pairwise_corr_median = H_pairwise_corr_median;
hStats.H_pairwise_corr = H_corr;

hStats.H_active_count_trace = activeCount;

end