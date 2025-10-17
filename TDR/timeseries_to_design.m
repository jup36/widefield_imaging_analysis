function [X, meta] = timeseries_to_design(timeseriesC, timeC, varargin)
%TIMESERIES_TO_DESIGN  Bin per-trial continuous time series into GLM windows.
%
% [X, meta] = TIMESERIES_TO_DESIGN(timeseriesC, timeC, ...)
%
% (Docstring unchanged for brevity — see earlier version.)

% -------- parse args --------
p = inputParser;
p.addParameter('Epoch', [-0.9 5], @(v)isnumeric(v)&&numel(v)==2);
p.addParameter('Win',   0.100,     @(v)isnumeric(v)&&isscalar(v)&&v>0);
p.addParameter('Step',  0.050,     @(v)isnumeric(v)&&isscalar(v)&&v>0);
p.addParameter('Method','mean',    @(s)ischar(s)||isstring(s));
p.addParameter('Fill',  NaN,       @(v)isnumeric(v)&&isscalar(v));
p.addParameter('zscore', true,     @(v)islogical(v)||ismember(v,[0 1]));
p.addParameter('minSD', 1e-3,      @(v)isnumeric(v)&&isscalar(v)&&v>0);
p.addParameter('clipExtremes', true, @(v)islogical(v)||ismember(v,[0 1]));
p.addParameter('zThresh', 20,      @(v)isnumeric(v)&&isscalar(v)&&v>0);
p.parse(varargin{:});
prm = p.Results;

N = numel(timeseriesC);
assert(N>0, 'timeseries_to_design: Empty input.');

% ensure timeseries and time vectors are aligned in length
ts_len = unique(cellfun(@(a) size(a, 2), timeseriesC));
assert(isscalar(ts_len), 'All trials must have the same # of samples per trial.');
if iscell(timeC)
    t_len = unique(cellfun(@(a) size(a, 2), timeC));
    assert(isscalar(t_len) && t_len==ts_len, 'timeC lengths must match timeseriesC lengths.');
elseif isnumeric(timeC)
    assert(numel(timeC)==ts_len, 'Shared time vector length must match per-trial series length.');
end

% windows
t0 = prm.Epoch(1); t1 = prm.Epoch(2);
winCtrs   = (t0 + prm.Win/2) : prm.Step : (t1 - prm.Win/2);
nW        = numel(winCtrs);
winBounds = [winCtrs(:)-prm.Win/2, winCtrs(:)+prm.Win/2];

% allow timeC to be a single numeric vector for all trials
useSharedTime = isnumeric(timeC) && isvector(timeC);
if ~useSharedTime
    assert(iscell(timeC) && numel(timeC)==N, ...
        'timeC must be a 1xN cell (timestamps per trial) or a single numeric vector.');
end

X = nan(N, nW);
validMask = false(N, nW);

% aggregator
switch lower(string(prm.Method))
    case "mean",   aggFun = @(x) mean(x, 'omitnan');
    case "sum",    aggFun = @(x) sum(x, 'omitnan');
    case "median", aggFun = @(x) median(x, 'omitnan');
    otherwise, error('Unsupported Method: %s', prm.Method);
end

% -------- main loop --------
for n = 1:N
    v = double(timeseriesC{n}(:));
    if isempty(v), continue; end
    
    if useSharedTime
        t = double(timeC(:));
    else
        t = double(timeC{n}(:));
    end
    assert(numel(t)==numel(v), 'Trial %d: time and value vectors must match in length.', n);
   
    for j = 1:nW
        lb = winBounds(j,1); ub = winBounds(j,2);
        idx = (t >= lb) & (t < ub);
        if any(idx)
            X(n,j) = aggFun(v(idx));
            validMask(n,j) = true;
        else
            X(n,j) = prm.Fill;
        end
    end
end

% -------- optional global z-score --------
if prm.zscore
    mu = mean(X(:), 'omitnan');
    sd = std(X(:), 0, 'omitnan');
    if ~isfinite(sd) || sd < prm.minSD, sd = prm.minSD; end
    X = (X - mu) ./ sd;
end

% -------- optional extreme-value masking (by z-score) --------
if prm.clipExtremes
    if prm.zscore
        Z = X;  % already standardized
    else
        mu = mean(X(:), 'omitnan');
        sd = std(X(:), 0, 'omitnan');
        if ~isfinite(sd) || sd < prm.minSD, sd = prm.minSD; end
        Z = (X - mu) ./ sd;
    end
    bad = isfinite(Z) & abs(Z) > prm.zThresh;
    X(bad) = NaN;
end

% -------- pack meta --------
meta = struct( ...
    'winCtrs',   winCtrs, ...
    'winBounds', winBounds, ...
    'validMask', validMask, ...
    'params',    rmfield(prm, [{'Method','Fill'}]) ...
    );
end
