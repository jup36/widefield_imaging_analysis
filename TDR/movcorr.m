function r = movcorr(x, y, w, varargin)
%MOVCORR  Moving window correlation (compatible with pre-R2023b MATLAB)
%
%   r = movcorr(x, y, w)
%
%   Computes Pearson correlation between x and y in a sliding window
%   of width w (in samples).  Returns vector of same length as x.
%
%   Optional 'Endpoints' behavior:
%       'shrink'  – compute correlation using smaller windows near edges (default)
%       'discard' – return NaN near edges
%
%   (Lightweight stand-in for MATLAB's built-in movcorr.)

    if nargin < 3
        error('movcorr(x, y, w) requires three inputs');
    end
    if isempty(x) || isempty(y)
        r = NaN(size(x));
        return;
    end

    x = x(:); y = y(:);
    N = numel(x);
    r = NaN(N,1);

    % parse optional 'Endpoints'
    endpoints = 'shrink';
    if nargin > 3 && ischar(varargin{1})
        endpoints = lower(varargin{1});
    end

    halfw = floor(w/2);
    for i = 1:N
        i1 = max(1, i - halfw);
        i2 = min(N, i + halfw);
        xi = x(i1:i2);
        yi = y(i1:i2);
        if numel(xi) > 2 && all(isfinite(xi)) && all(isfinite(yi))
            r(i) = corr(xi, yi, 'rows', 'pairwise');
        elseif strcmp(endpoints,'discard')
            r(i) = NaN;
        end
    end
end
