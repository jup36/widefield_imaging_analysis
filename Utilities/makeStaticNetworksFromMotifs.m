function W_basis_static = makeStaticNetworksFromMotifs(W_basis, varThresh)
% makeStaticNetworksFromMotifs
%
% Converts dynamic spatiotemporal motifs into static networks.
%
% For each motif:
%   1. Identify active lags as lags with variance across pixels > varThresh.
%   2. Compute the mean spatial activation pattern across active lags.
%   3. Replace each active lag with that mean spatial pattern.
%   4. Leave inactive lags as zero.
%
% INPUT
%   W_basis   : pixels × motifs × lags
%   varThresh : variance threshold for defining active lags
%               default = 0
%
% OUTPUT
%   W_basis_static : pixels × motifs × lags

if nargin < 2 || isempty(varThresh)
    varThresh = 0;
end

[P, K, L] = size(W_basis); %#ok<ASGLU>

W_basis_static = zeros(size(W_basis), 'like', W_basis);

for k = 1:K

    Wk = squeeze(W_basis(:, k, :));  % pixels × lags

    % Active time points are lags with variance across pixels > varThresh
    lagVar = var(Wk, 0, 1, 'omitnan');
    activeLags = lagVar > varThresh;

    if ~any(activeLags)
        continue
    end

    % Mean activation pattern across active lags
    meanMap = mean(Wk(:, activeLags), 2, 'omitnan');  % pixels × 1

    % Replace each active lag with the same static map
    for l = find(activeLags)
        W_basis_static(:, k, l) = meanMap;
    end

end

end