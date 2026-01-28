function rezInterp = interp_rez_to_hAligned(tbytDat_hAligned, rezStim)
%INTERP_REZ_TO_HALIGNED
% Interpolate sparse rez time series onto the time grid used in
% tbytDat_hAligned.
%
% INPUTS
%   tbytDat_hAligned : 2 x N cell
%       row 1: K x T matrices (motifs x time) [not used except for T]
%       row 2: 1 x T or T x 1 time vector (e.g., -0.9:0.01:5)
%
%   rezStim : N x 2 cell
%       col 1: signal values (vector or [D x T_sparse])
%       col 2: corresponding time vector (T_sparse)
%
% OUTPUT
%   rezInterp : 2 x N cell
%       row 1: interpolated signal (same rows as rezStim{:,1}, length = T)
%       row 2: target time vector (copied from tbytDat_hAligned)
%
% NOTES
%   - Empty entries are skipped and returned as empty.
%   - Uses linear interpolation with NaN extrapolation.
%   - Assumes rezStim{n,1} varies along time dimension 2 if matrix.

% ---- basic checks ----
assert(size(tbytDat_hAligned,1) == 2, 'tbytDat_hAligned must be 2 x N');
assert(size(rezStim,2) == 2, 'rezStim must be N x 2');

N = size(tbytDat_hAligned, 2);

rezInterp = cell(2, N);

for n = 1:N

    % target time grid
    tTarget = tbytDat_hAligned{2, n};

    % propagate time axis regardless
    rezInterp{2, n} = tTarget;

    % skip empty trials
    if isempty(rezStim{n,1}) || isempty(rezStim{n,2}) || isempty(tTarget)
        rezInterp{1, n} = [];
        continue
    end

    y = rezStim{n,1};
    t = rezStim{n,2};

    % enforce column time vectors
    t = t(:);
    tTarget = tTarget(:);

    % ensure y is [D x T]
    if isvector(y)
        y = y(:).';   % 1 x T
    end

    % interpolate along time dimension
    try
        yInterp = interp1( ...
            t, ...
            y.', ...                 % T x D
            tTarget, ...
            'linear', ...
            NaN ...
        ).';                          % D x T

    catch ME
        warning('Interpolation failed at trial %d: %s', n, ME.message);
        yInterp = [];
    end

    rezInterp{1, n} = yInterp;
end
end
