function pev = i_compute_pev(y, yhat, varargin)
% i_compute_pev
%
% Computes global fraction of variance explained:
%
%   PEV = 1 - SSE / SST
%
% Inputs:
%   y     : observed data, e.g. dffTr, size X x Y x T
%   yhat  : reconstructed data, e.g. whTr, size X x Y x T
%
% Name-value:
%   'asPercent' : true/false, default false
%
% Output:
%   pev : fraction explained variance, or percent if asPercent = true

% --------------------
% Parse inputs
% --------------------
p = inputParser;
p.addParameter('asPercent', false, @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});

asPercent = logical(p.Results.asPercent);

% --------------------
% Basic checks
% --------------------
pev = nan;

if isempty(y) || isempty(yhat)
    warning('Input y or yhat is empty.');
    return
end

if ~isequal(size(y), size(yhat))
    warning('Input y and yhat must have the same size.');
    return
end

% --------------------
% Convert to double
% --------------------
y = double(y);
yhat = double(yhat);

% --------------------
% Use only finite values present in both y and yhat
% --------------------
valid = isfinite(y) & isfinite(yhat);

if ~any(valid(:))
    warning('No valid finite samples found.');
    return
end

yv = y(valid);
yhv = yhat(valid);

% --------------------
% Compute global PEV across all pixels and frames
% --------------------
ss_res = sum((yv - yhv).^2);
ss_tot = sum((yv - mean(yv)).^2);

if ss_tot <= 0
    warning('Total variance is zero or negative.');
    return
end

pev = 1 - ss_res / ss_tot;

if asPercent
    pev = pev * 100;
end

end