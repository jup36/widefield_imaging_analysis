function rezInterp = interp_rez_to_target_time(rezStim, tTarget)
%INTERP_REZ_TO_TARGET_TIME
% Interpolate sparse rez time series onto a fixed target time grid.
%
% INPUTS
%   rezStim : N x 2 cell
%       col 1: signal values (vector or [D x T_sparse])
%       col 2: corresponding time vector (T_sparse)
%
%   tTarget : 1 x T or T x 1 numeric vector
%       Target timestamps (e.g., -0.9:0.01:5)
%
% OUTPUT
%   rezInterp : 2 x N cell
%       row 1: interpolated signal (same rows as rezStim{:,1}, length = T)
%       row 2: target timestamps (copied from tTarget)
%
% NOTES
%   - Empty entries are skipped and returned as empty.
%   - Uses linear interpolation with NaN extrapolation.
%   - Interpolates along time dimension (dim 2).

assert(size(rezStim,2) == 2, 'rezStim must be N x 2 cell array');

N = size(rezStim, 1);
tTarget = tTarget(:);    % force column

rezInterp = cell(2, N);

for n = 1:N

    % skip empty trials entirely
    if isempty(rezStim{n,1}) || isempty(rezStim{n,2})
        rezInterp{1, n} = [];
        rezInterp{2, n} = [];
        continue
    end

    % fill timestamps ONLY when data exist
    rezInterp{2, n} = tTarget(:).';   % row vector

    y = rezStim{n,1};
    t = rezStim{n,2};

    t = t(:);   % column

    if isvector(y)
        y = y(:).';   % enforce row
    end

    % interpolate (NaN outside support)
    yInterp = interp1(t, y.', tTarget(:), 'linear', NaN).';
    rezInterp{1, n} = yInterp;
end

end
