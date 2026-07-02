function rez = computeTrialLevelSameMotifXcorr(hAligned, daAligned, trialMask, varargin)
% computeTrialLevelSameMotifXcorr
%
% Computes same-motif trial-level H-DA xcorr.
%
% Inputs:
%   hAligned{1, tr}  = K x time H
%   daAligned{1, tr} = K x time DA
%
% Positive lag means H leads DA:
%   corr(H(t), DA(t + lag))

p = inputParser;
addParameter(p, 'maxLagBins', 200, @isnumeric);
addParameter(p, 'dt', 0.01, @isnumeric);
addParameter(p, 'subtractConditionMean', false, @islogical);
parse(p, varargin{:});

maxLagBins = p.Results.maxLagBins;
dt = p.Results.dt;
subtractConditionMean = p.Results.subtractConditionMean;

trialIdx = find(trialMask(:)');

% Identify valid trials with non-empty H and DA.
validTrials = [];

for ii = 1:numel(trialIdx)

    tr = trialIdx(ii);

    if tr > size(hAligned, 2) || tr > size(daAligned, 2)
        continue;
    end

    if isempty(hAligned{1, tr}) || isempty(daAligned{1, tr})
        continue;
    end

    H = hAligned{1, tr};
    DA = daAligned{1, tr};

    if ndims(H) ~= 2 || ndims(DA) ~= 2
        continue;
    end

    if size(H, 1) ~= size(DA, 1)
        continue;
    end

    if size(H, 2) ~= size(DA, 2)
        continue;
    end

    validTrials(end+1) = tr; %#ok<AGROW>
end

if isempty(validTrials)
    warning('No valid trials found for this condition.');
    rez = emptyTrialXcorrRez(maxLagBins, dt);
    return;
end

% Use first valid trial to define K and T.
H0 = hAligned{1, validTrials(1)};
[K, T] = size(H0);

nTrials = numel(validTrials);

H3 = NaN(K, T, nTrials);
DA3 = NaN(K, T, nTrials);

for ii = 1:nTrials

    tr = validTrials(ii);

    H = double(hAligned{1, tr});
    DA = double(daAligned{1, tr});

    H3(:, :, ii) = H;
    DA3(:, :, ii) = DA;
end

% Optional: subtract condition-average PETH.
% This asks about trial-by-trial residual coupling, not just shared event-locking.
if subtractConditionMean
    Hmean = mean(H3, 3, 'omitnan');
    DAmean = mean(DA3, 3, 'omitnan');

    for ii = 1:nTrials
        H3(:, :, ii) = H3(:, :, ii) - Hmean;
        DA3(:, :, ii) = DA3(:, :, ii) - DAmean;
    end
end

lagBins = -maxLagBins:maxLagBins;
lagSec = lagBins * dt;

rTrial = NaN(K, nTrials, numel(lagBins));
peakR = NaN(K, nTrials);
peakLagSec = NaN(K, nTrials);

for ii = 1:nTrials

    for k = 1:K

        h = squeeze(H3(k, :, ii));
        da = squeeze(DA3(k, :, ii));

        h = zscore_omitnan(h);
        da = zscore_omitnan(da);

        if all(~isfinite(h)) || all(~isfinite(da))
            continue;
        end

        [rLag, ~] = lagCorr_H_leads_DA(h, da, maxLagBins);

        rTrial(k, ii, :) = rLag;

        if all(~isfinite(rLag))
            continue;
        end

        [~, maxI] = max(abs(rLag));

        peakR(k, ii) = rLag(maxI);
        peakLagSec(k, ii) = lagSec(maxI);
    end
end

meanR = squeeze(mean(rTrial, 2, 'omitnan'));
semR = squeeze(std(rTrial, 0, 2, 'omitnan')) ./ sqrt(squeeze(sum(isfinite(rTrial), 2)));

nTrialsPerMotif = squeeze(sum(any(isfinite(rTrial), 3), 2));

rez = struct();
rez.rTrial = rTrial;
rez.meanR = meanR;
rez.semR = semR;
rez.peakR = peakR;
rez.peakLagSec = peakLagSec;
rez.lagBins = lagBins;
rez.lagSec = lagSec;
rez.validTrials = validTrials;
rez.nTrials = nTrials;
rez.nTrialsPerMotif = nTrialsPerMotif;
rez.subtractConditionMean = subtractConditionMean;

end