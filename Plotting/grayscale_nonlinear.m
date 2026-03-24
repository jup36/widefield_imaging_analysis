function cmap = grayscale_nonlinear(n, method, param, varargin)
%GRAYSCALE_NONLINEAR Create a nonlinear grayscale colormap.
%
% cmap = grayscale_nonlinear(n, method, param)
% cmap = grayscale_nonlinear(..., 'Reverse', true/false)
%
% PURPOSE
%   Create a grayscale colormap where intensity is mapped nonlinearly,
%   so you can enhance contrast in selected value ranges without changing
%   the underlying data.
%
% INPUTS
%   n      : number of colors (default = 256)
%   method : 'linear' | 'gamma' | 'logistic' | 'tail'
%            default = 'gamma'
%   param  : parameter controlling curvature
%            For 'linear'   : ignored
%            For 'gamma'    : gamma exponent (>1 emphasizes larger values)
%            For 'logistic' : steepness k (e.g. 8 or 10)
%            For 'tail'     : alpha (>1 emphasizes larger values)
%            default depends on method
%
% NAME-VALUE
%   'Reverse'  : false (default)
%   'Midpoint' : 0.7 (used only for 'logistic')
%
% OUTPUT
%   cmap : [n x 3] grayscale colormap
%
% EXAMPLES
%   colormap(grayscale_nonlinear(256, 'gamma', 2.2));
%   colormap(grayscale_nonlinear(256, 'logistic', 10, 'Midpoint', 0.75));
%   colormap(grayscale_nonlinear(256, 'tail', 2.5));
%

if nargin < 1 || isempty(n), n = 256; end
if nargin < 2 || isempty(method), method = 'gamma'; end
if nargin < 3, param = []; end

% defaults
reverseFlag = false;
midpoint = 0.7;

% lightweight name-value parsing
if ~isempty(varargin)
    assert(mod(numel(varargin),2)==0, 'Name-value args must come in pairs.');
    for i = 1:2:numel(varargin)
        key = lower(string(varargin{i}));
        val = varargin{i+1};
        switch key
            case "reverse"
                reverseFlag = logical(val);
            case "midpoint"
                midpoint = val;
            otherwise
                error('Unknown name-value argument: %s', key);
        end
    end
end

method = lower(string(method));
x = linspace(0,1,n)';

switch method
    case "linear"
        y = x;

    case "gamma"
        if isempty(param), param = 2; end
        gamma = param;
        assert(gamma > 0, 'For gamma mode, param must be > 0.');
        y = x.^gamma;

    case "logistic"
        if isempty(param), param = 10; end
        k = param;
        assert(k > 0, 'For logistic mode, param must be > 0.');
        assert(midpoint > 0 && midpoint < 1, 'Midpoint must be between 0 and 1.');

        y = 1 ./ (1 + exp(-k * (x - midpoint)));

        % rescale to [0,1]
        y = (y - min(y)) ./ max(max(y) - min(y), eps);

    case "tail"
        if isempty(param), param = 2; end
        alpha = param;
        assert(alpha > 0, 'For tail mode, param must be > 0.');
        y = 1 - (1 - x).^alpha;

    otherwise
        error('Unknown method: %s. Use linear, gamma, logistic, or tail.', method);
end

if reverseFlag
    y = flipud(y);
end

cmap = [y y y];
end