function [X, names, F] = convolve_events_with_basis(E, basis, baseName, varargin)
%CONVOLVE_EVENTS_WITH_BASIS  Convolve event trains with a basis bank.
% Optionally re-center kernels so lag==0 is at the conv anchor (center).
%
% Inputs
%   E               : [N x nW] event train (trials x time bins)
%   basis.B         : [nLag x nB] basis bank (columns = kernels)
%   basis.lagsSec   : [nLag x 1] lag axis (MUST contain 0 exactly)
%   baseName        : (optional) string base for column names
%
% Name-Value (optional)
%   'centerZeroBasis' (false) : if true, zero-pad shift each kernel so that
%                               the sample at lag==0 is at the center index
%                               used by conv(...,'same').
%   'assertCentered'  (false) : if true AND centerZeroBasis=false, assert
%                               that 0-lag is already at the center index.
%
% Outputs
%   X     : [(N*nW) x nB] design matrix in bin-major order
%   names : 1 x nB cell array of regressor names
%   F     : [N x nW x nB] per-trial, per-basis convolved predictors
%
% Notes
%   - Zero-padding (no wrap) is used when re-centering to avoid artifacts.
%   - Prefer odd nLag so the center index is unambiguous.
%   - Bin-major packing: rows 1..N = bin 1 (all trials), etc.

% ---------------- options ----------------
p = inputParser;
p.addParameter('centerZeroBasis', false, @(v) islogical(v) || ismember(v,[0 1]));
p.addParameter('assertCentered',   false, @(v) islogical(v) || ismember(v,[0 1]));
p.parse(varargin{:});
opt = p.Results;

% ---------------- basics ----------------
[N, nW] = size(E);
K    = basis.B;                 % [nLag x nB]
lags = basis.lagsSec(:);        % [nLag x 1]
[nLag, nB] = size(K);

% sanity: ensure 0-lag exists on the provided grid
idx0 = find(lags == 0, 1, 'first');
if isempty(idx0)
    error('convolve_events_with_basis: basis.lagsSec must contain 0 exactly.');
end

centerIdx = ceil((nLag + 1)/2);

% ---------------- optional (re)centering ----------------
if opt.centerZeroBasis
    % zero-padded shift so lag==0 sits at conv anchor
    Kc = zeros(size(K));
    for b = 1:nB
        Kc(:,b) = center_kernel_linear(K(:,b), idx0, centerIdx);
    end
else
    % assume already centered; optionally assert
    if opt.assertCentered && idx0 ~= centerIdx
        error(['Basis not centered: lag==0 at idx %d, center at %d. ' ...
               'Either rebuild centered or call with ''centerZeroBasis'',true.'], idx0, centerIdx);
    end
    Kc = K;
end

% ---------------- convolve per trial, per basis ----------------
E = double(E);
F = zeros(N, nW, nB);
for b = 1:nB
    k = Kc(:,b);
    for n = 1:N
        F(n,:,b) = conv(E(n,:), k, 'same');  % conv anchor = kernel center
    end
end

% ---------------- pack to bin-major ----------------
X = reshape(F, N*nW, nB);

% ---------------- names ----------------
if nargin < 3 || isempty(baseName), baseName = "event"; end
names = arrayfun(@(b) sprintf('%s_rc%02d', baseName, b), 1:nB, 'uni', 0);

end % main


% ===== helper (zero-padded shift) =====
function kc = center_kernel_linear(k, idx0, centerIdx)
%CENTER_KERNEL_LINEAR  Shift kernel so sample at lag==0 lands at centerIdx.
% Zero-pad instead of wrapping to keep positive-lag energy on the right.
    L = numel(k);
    s = centerIdx - idx0;   % desired shift in samples (+ = right)
    if s > 0
        % shift right by s: prepend s zeros, drop last s samples
        kc = [zeros(s,1); k(1:L-s)];
    elseif s < 0
        % shift left by |s|: drop first |s| samples, append zeros
        s = -s;
        kc = [k(1+s:L); zeros(s,1)];
    else
        kc = k;
    end
end
