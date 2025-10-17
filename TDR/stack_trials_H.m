function out = stack_trials_H(tbytDat_hAligned, varargin)
%STACK_TRIALS_H  Bin & stack per-trial motif activity (H) across time.
%
% out = STACK_TRIALS_H(tbytDat_hAligned, 'Epoch', [t0 t1], 'Win', 0.1, 'Step', 0.05, ...)
%
% INPUT
%   tbytDat_hAligned : 2xN cell array per trial {H (K x Tn), t (1 x Tn)}.
%
% NAME-VALUE OPTIONS
%   'Epoch'       [t0 t1]   time seconds to analyze (default: [min_t max_t] across trials)
%   'Win'         scalar    bin width in seconds (default: 0.100)
%   'Step'        scalar    step in seconds (default: 0.050)
%   'motifIdx'    vector    subset of motifs (rows of H) to keep (default: all)
%   'zscore'      logical   z-score each motif using global mean/std across all bins (default: false)
%   'minPerBin'   integer   min #valid trials required to keep a bin (masking only; default: 1)
%   'minSD'       scalar    minimum SD (default: 1e-3)
%
% OUTPUT (struct)
%   out.Yw        1 x nW cell, each N x K matrix of window-averaged H (NaNs where no samples)
%   out.Y3        N x K x nW array, same as Yw but stacked (NaN-padded)
%   out.validMask N x nW logical, true if trial contributed to that time bin
%   out.winCtrs   1 x nW vector of window centers (s)
%   out.winBounds nW x 2 window edges [lb ub] (s)
%   out.params    struct with Epoch/Win/Step/motifIdx/zscore
%
% NOTES
%   - A "valid" trial for a bin is one with at least 1 sample inside [lb, ub).
%   - If 'zscore' is true, z-scoring is done after windowing using global stats
%     across all trials/bins (per motif).
%
% Junchol Park / Buschman Lab — 2025

% -------------------- Parse args --------------------
p = inputParser;
p.addParameter('Epoch', [], @(v)isnumeric(v)&&numel(v)==2);
p.addParameter('Win', 0.100, @(v)isnumeric(v)&&isscalar(v)&&v>0);
p.addParameter('Step', 0.050, @(v)isnumeric(v)&&isscalar(v)&&v>0);
p.addParameter('motifIdx', [], @(v)isnumeric(v)&&isvector(v));
p.addParameter('zscore', true, @(v)islogical(v)||ismember(v,[0 1]));
p.addParameter('minPerBin', 1, @(v)isnumeric(v)&&isscalar(v)&&v>=0);
p.addParameter('minSD', 1e-3, @(v)isnumeric(v)&&isscalar(v)&&v>0);     % variance floor for unscaling
p.parse(varargin{:});
prm = p.Results;

% -------------------- Basics --------------------
N = size(tbytDat_hAligned, 2);
K_all = size(tbytDat_hAligned{1,1}, 1);
motifIdx = ternary(isempty(prm.motifIdx), (1:K_all), prm.motifIdx(:)');
K = numel(motifIdx);

% Epoch default from available timestamps
if isempty(prm.Epoch)
    tmin = +inf; tmax = -inf;
    for n = 1:N
        t = tbytDat_hAligned{2,n};
        if isempty(t), continue; end
        tmin = min(tmin, min(t));
        tmax = max(tmax, max(t));
    end
    prm.Epoch = [tmin tmax];
end

t0 = prm.Epoch(1); t1 = prm.Epoch(2);
ctr = (t0 + prm.Win/2):prm.Step:(t1 - prm.Win/2);
nW  = numel(ctr);
winBounds = [ctr(:)-prm.Win/2, ctr(:)+prm.Win/2];

% -------------------- Window & stack --------------------
Yw = cell(1, nW);              % each: N x K (NaN where empty)
validMask = false(N, nW);

for j = 1:nW
    lb = winBounds(j,1); ub = winBounds(j,2);
    Y = nan(N, K);

    for n = 1:N
        t = tbytDat_hAligned{2,n};
        if isempty(t), continue; end

        idx = (t >= lb) & (t <  ub);
        if any(idx)
            H = tbytDat_hAligned{1,n};
            H = H(motifIdx, :);
            Y(n,:) = mean(H(:, idx), 2, 'omitnan')';
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
    % stack all Y into (N*nW) x K with NaNs
    Ystack = cat(1, Yw{:});
    muG = mean(Ystack, 1, 'omitnan');
    sdG = std( Ystack, 0, 1, 'omitnan');
    sdG(sdG < prm.minSD | ~isfinite(sdG)) = prm.minSD;        % variance floor

    for j = 1:nW
        Y = Yw{j};
        Yw{j} = (Y - muG) ./ sdG;
    end
end

% -------------------- Pack 3D array as convenience --------------------
Y3 = nan(N, K, nW);
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
                        'motifIdx',motifIdx,'zscore',prm.zscore,'minPerBin',prm.minPerBin) );

end

% --- tiny helper (local) ---
function y = ternary(cond, a, b)
if cond, y = a; else, y = b; end
end
