function [B, lagsSec, info] = make_rcos_basis(nBases, lagRangeSec, binSec, varargin)
%MAKE_RCOS_BASIS  Create raised-cosine temporal basis functions (can include negative lags).
%
% [B, lagsSec, info] = make_rcos_basis(nBases, [lagMin lagMax], binSec, ...)
%
% INPUTS
%   nBases       : number of basis functions
%   lagRangeSec  : [lagMin lagMax] in seconds (can include negative values)
%   binSec       : temporal bin size (e.g., your GLM Step)
%
% NAME-VALUE (optional)
%   'nonlin'     : 'linear' | 'log' (default 'linear') - spacing of centers
%   'c'          : width scaling factor (default 1)
%   'normCols'   : true/false (default true)
%
% OUTPUTS
%   B        : L x nBases matrix of basis functions (each column = one basis)
%   lagsSec  : L x 1 vector of lag times (s)
%   info     : struct with centers, width, and other metadata
%
% NOTES
%   - Works for both causal ([0, +T]) and acausal ([-T1, +T2]) kernels.
%   - Each basis function is raised cosine centered at evenly spaced points.
%
% Example:
%   [B, lags] = make_rcos_basis(8, [-1 2], 0.05);
%   plot(lags, B); xlabel('Lag (s)'); ylabel('Basis amplitude');

p = inputParser;
p.addParameter('nonlin', 'linear', @(s) any(strcmpi(s, {'linear','log'})));
p.addParameter('c', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('normCols', true, @(v) islogical(v) || ismember(v,[0 1]));
p.parse(varargin{:});
opt = p.Results;

lagMin = lagRangeSec(1);
lagMax = lagRangeSec(2);
if lagMax <= lagMin
    error('lagMax must be greater than lagMin.');
end

% lag axis matching your GLM binning
lagsSec = (lagMin:binSec:lagMax)';
L = numel(lagsSec);

% centers of raised cosines
switch lower(opt.nonlin)
    case 'linear'
        ctrSec = linspace(lagMin, lagMax, nBases);
    case 'log'
        % ensure strictly positive range for log warping
        eps0 = max(binSec, 1e-3);
        shift = abs(lagMin) + eps0;
        u = linspace(log(eps0), log(lagMax + shift), nBases);
        ctrSec = exp(u) - shift;
end

% width so adjacent bases overlap by ~50%
if nBases == 1
    w = (lagMax - lagMin) / 2;
else
    w = opt.c * (ctrSec(2) - ctrSec(1));
end

B = zeros(L, nBases);
for k = 1:nBases
    x = (lagsSec - ctrSec(k)) / (w/2);
    bump = 0.5 * (cos(max(-pi, min(pi, x))) + 1);
    bump(abs(x) > pi) = 0;
    B(:,k) = bump;
end

if opt.normCols
    B = B ./ max(sqrt(sum(B.^2,1)), eps);
end

info = struct('centersSec', ctrSec, 'width', w, ...
              'lagsSec', lagsSec, 'nonlin', opt.nonlin, 'c', opt.c);
end
