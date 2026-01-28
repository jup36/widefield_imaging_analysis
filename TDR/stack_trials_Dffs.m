function out = stack_trials_Dffs(dff_itp, varargin)
%STACK_TRIALS_DFFS  Bin & stack per-trial DFF time series across time.
%
% out = STACK_TRIALS_DFFS(dff_itp, 'Epoch', [t0 t1], 'Win', 0.1, 'Step', 0.05, ...)
%
% INPUT
%   dff_itp : 2xN cell array per trial {dff (1 x Tn OR Tn x 1), t (1 x Tn)}
%             Can contain empty entries.
%
% NAME-VALUE OPTIONS
%   'Epoch'       [t0 t1]   time seconds to analyze (default: [min_t max_t] across trials)
%   'Win'         scalar    bin width in seconds (default: 0.100)
%   'Step'        scalar    step in seconds (default: 0.050)
%   'zscore'      logical   z-score using global mean/std across all bins (default: false)
%   'minPerBin'   integer   min #valid trials required to keep a bin (masking only; default: 1)
%   'minSD'       scalar    minimum SD (default: 1e-3)
%
% OUTPUT (struct)
%   out.Yw        1 x nW cell, each N x 1 vector of window-averaged dff (NaNs where no samples)
%   out.Y3        N x 1 x nW array, same as Yw but stacked (NaN-padded)
%   out.validMask N x nW logical, true if trial contributed to that time bin
%   out.winCtrs   1 x nW vector of window centers (s)
%   out.winBounds nW x 2 window edges [lb ub] (s)
%   out.params    struct with Epoch/Win/Step/zscore/minPerBin
%
% NOTES
%   - A "valid" trial for a bin is one with at least 1 sample inside [lb, ub).
%   - Robust to empty dff/t entries (skips them).

% -------------------- Parse args --------------------
p = inputParser;
p.addParameter('Epoch', [], @(v)isnumeric(v)&&numel(v)==2);
p.addParameter('Win', 0.100, @(v)isnumeric(v)&&isscalar(v)&&v>0);
p.addParameter('Step', 0.050, @(v)isnumeric(v)&&isscalar(v)&&v>0);
p.addParameter('zscore', false, @(v)islogical(v)||ismember(v,[0 1]));
p.addParameter('minPerBin', 1, @(v)isnumeric(v)&&isscalar(v)&&v>=0);
p.addParameter('minSD', 1e-3, @(v)isnumeric(v)&&isscalar(v)&&v>0);
p.parse(varargin{:});
prm = p.Results;

% -------------------- Basics --------------------
if size(dff_itp,1) ~= 2
    error('stack_trials_Dffs: input must be 2xN cell array {dff; t}.');
end
N = size(dff_itp, 2);

% -------------------- Epoch default from available timestamps --------------------
if isempty(prm.Epoch)
    tmin = +inf; tmax = -inf;
    for n = 1:N
        t = dff_itp{2,n};
        if isempty(t), continue; end
        t = t(:);
        if all(~isfinite(t)), continue; end
        tmin = min(tmin, min(t,[],'omitnan'));
        tmax = max(tmax, max(t,[],'omitnan'));
    end
    if ~isfinite(tmin) || ~isfinite(tmax)
        error('stack_trials_Dffs: could not infer Epoch (all timestamps empty/invalid).');
    end
    prm.Epoch = [tmin tmax];
end

t0 = prm.Epoch(1); t1 = prm.Epoch(2);
ctr = (t0 + prm.Win/2):prm.Step:(t1 - prm.Win/2);
nW  = numel(ctr);
winBounds = [ctr(:)-prm.Win/2, ctr(:)+prm.Win/2];

% -------------------- Window & stack --------------------
Yw = cell(1, nW);              % each: N x 1 (NaN where empty)
validMask = false(N, nW);

for j = 1:nW
    lb = winBounds(j,1); ub = winBounds(j,2);
    Y = nan(N, 1);

    for n = 1:N
        t = dff_itp{2,n};
        x = dff_itp{1,n};

        if isempty(t) || isempty(x)
            continue;
        end

        t = t(:); % column
        x = x(:); % column

        % require matching lengths
        if numel(t) ~= numel(x)
            % try to be forgiving if one is row and the other is column but same numel (handled)
            warning('stack_trials_Dffs: trial %d has mismatched lengths (t=%d, x=%d); skipping.', ...
                n, numel(t), numel(x));
            continue;
        end

        idx = (t >= lb) & (t < ub) & isfinite(x) & isfinite(t);
        if any(idx)
            Y(n,1) = mean(x(idx), 'omitnan');
            validMask(n,j) = true;
        end
    end

    Yw{j} = Y;
end

% enforce minPerBin (mask only; data remain as-is)
if prm.minPerBin > 1
    keep = sum(validMask, 1) >= prm.minPerBin;
    validMask(:, ~keep) = false;
end

% -------------------- Global z-scoring (optional) --------------------
if prm.zscore
    Ystack = cat(1, Yw{:}); % (N*nW) x 1
    muG = mean(Ystack, 1, 'omitnan');
    sdG = std( Ystack, 0, 1, 'omitnan');
    if ~isfinite(sdG) || sdG < prm.minSD, sdG = prm.minSD; end

    for j = 1:nW
        Yw{j} = (Yw{j} - muG) ./ sdG;
    end
end

% -------------------- Pack 3D array --------------------
Y3 = nan(N, 1, nW);
for j = 1:nW
    Y3(:,:,j) = Yw{j};
end

out = struct( ...
    'Yw',        {Yw}, ...
    'Y3',        Y3, ...
    'validMask', validMask, ...
    'winCtrs',   ctr, ...
    'winBounds', winBounds, ...
    'params',    struct('Epoch',prm.Epoch,'Win',prm.Win,'Step',prm.Step, ...
                        'zscore',logical(prm.zscore),'minPerBin',prm.minPerBin) );

end
