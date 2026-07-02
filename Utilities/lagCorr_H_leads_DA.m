
function [rLag, lagBins] = lagCorr_H_leads_DA(H, DA, maxLagBins)
% lagCorr_H_leads_DA
%
% Positive lag means H leads DA:
%
%   lag > 0:
%       corr(H(t), DA(t + lag))
%
%   lag < 0:
%       corr(H(t), DA(t + lag))
%       where DA occurs earlier than H.
%
% Inputs should already be on a common time grid.

lagBins = -maxLagBins:maxLagBins;
rLag = NaN(size(lagBins));

H = double(H(:)');
DA = double(DA(:)');

n = min(numel(H), numel(DA));
H = H(1:n);
DA = DA(1:n);

for i = 1:numel(lagBins)

    lag = lagBins(i);

    if lag > 0

        hUse = H(1:n-lag);
        daUse = DA(1+lag:n);

    elseif lag < 0

        lagAbs = abs(lag);
        hUse = H(1+lagAbs:n);
        daUse = DA(1:n-lagAbs);

    else

        hUse = H;
        daUse = DA;
    end

    validI = isfinite(hUse) & isfinite(daUse);

    if sum(validI) < 10
        continue;
    end

    rLag(i) = corr(hUse(validI)', daUse(validI)');
end

end