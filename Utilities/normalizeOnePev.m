% ---------- helper ----------
function out = normalizeOnePev(pev, WnormCoef)
    % ensure row
    pev = pev(:).';
    n   = numel(pev);
    w   = WnormCoef(1:n);

    % guard against zeros in w
    w(w==0) = NaN;

    scaled = pev ./ w;                % reweight by W
    s = nansum(scaled);
    if s > 0
        out = scaled ./ s;            % renormalize to sum = 1
    else
        out = nan(size(pev));         % fallback if all zero/NaN
    end
end