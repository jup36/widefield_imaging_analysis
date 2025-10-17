function [X, names, F] = convolve_events_with_basis(E, basis, baseName)
%CONVOLVE_EVENTS_WITH_BASIS Convolve event trains with a raised-cosine bank.
%
% INPUTS
%   E      : N x nW binary/count events (trials × time bins)
%   basis  : struct with field .B_evt (nLag x nB). Each column is a kernel,
%            centered at 0 lag (negative→positive lags).
%   baseName (optional) : string for column names, e.g., 'evtOn'
%
% OUTPUTS
%   X    : (N*nW) x nB in **bin-major** row order:
%          rows 1..N   = bin 1 (all trials)
%          rows N+1..2N= bin 2 (all trials)
%          ...
%          Index mapping: r = (j-1)*N + n
%   names: 1 x nB cell array of column names
%   F    : N x nW x nB (per-trial, per-bin, per-basis) for inspection
%
% NOTE
%   If your response Y is also flattened as Y(:), this bin-major X aligns
%   directly with Y(:) in MATLAB (column-major) memory layout.

% ---- basics ----
[N, nW] = size(E);
K = basis.B_evt;                       % nLag x nB
[nLag, nB] = size(K); %#ok<ASGLU>      % nLag unused explicitly, but validates dims
E = double(E);

% ---- convolve per trial, per basis ----
F = zeros(N, nW, nB);
for b = 1:nB
    k = K(:,b);                        % nLag x 1, centered at 0
    for n = 1:N
        F(n,:,b) = conv(E(n,:), k, 'same');   % 1 x nW (handles edges)
    end
end

% ---- pack to bin-major (matches Y(:)) ----
% Bin j block is rows ((j-1)*N + 1) : (j*N).
X = reshape(F, N*nW, nB);              % no permute needed for bin-major

% ---- names ----
if nargin < 3 || isempty(baseName), baseName = "event"; end
names = arrayfun(@(b) sprintf('%s_rc%02d', baseName, b), 1:nB, 'uni', 0);
end
