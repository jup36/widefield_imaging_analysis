function [Wz_aligned, flipSign] = align_decoder_signs(Wz, ctr, refEpoch)
% Aligns decoder weights Wz across time bins based on a reference epoch
%   Wz       - K x T matrix of decoder weights (K: features, T: time bins)
%   ctr      - 1 x T vector of window centers (in seconds)
%   refEpoch - [start, end] reference period for sign alignment (e.g., [2, 4])

T = size(Wz,2);
Wz_aligned = Wz;
flipSign = false(1, T);

% Identify reference bins
refIdx = find(ctr >= refEpoch(1) & ctr <= refEpoch(2));
if isempty(refIdx)
    warning('No reference bins found in specified refEpoch.');
    return
end

% Compute reference vector (mean over ref bins)
refVec = mean(Wz(:,refIdx), 2, 'omitnan');
refVecNorm = norm(refVec(~isnan(refVec)));
if refVecNorm == 0
    warning('Reference vector has zero norm; skipping alignment.');
    return
end
refUnit = refVec / refVecNorm;

% Align each time bin
for tt = 1:T
    w = Wz(:,tt);
    if all(isnan(w)), continue; end
    val = isfinite(w) & isfinite(refUnit);
    if ~any(val), continue; end
    cosSim = dot(w(val), refUnit(val)) / (norm(w(val)) * norm(refUnit(val)));
    if cosSim < 0
        Wz_aligned(:,tt) = -w;
        flipSign(tt) = true;
    end
end
end