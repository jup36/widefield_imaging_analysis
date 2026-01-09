function hFig = imageHmmStateOccupancy(gammaC, trialId, time)
% imageHmmStateOccupancy
%
% Visualize HMM state occupancy (gamma) for a single trial using imagesc.
%
% Inputs:
%   gammaC  : 1×N cell array, each cell is S×T (state occupancy)
%   trialId : scalar index of the trial to visualize
%   time    : 1×T time vector (used for x-axis labeling)
%
% Output:
%   hFig    : figure handle
%
% Example:
%   hFig = imageHmmStateOccupancy(gammaC, 3, hmmInput.time);

% -------------------- sanity checks --------------------
assert(trialId >= 1 && trialId <= numel(gammaC), ...
    'trialId out of range');

gamma = gammaC{trialId};   % S×T
[S, T] = size(gamma);

assert(S == 4, 'Expected 4 states (S=4), got S=%d', S);
assert(numel(time) == T, 'Time vector length must match T');

% -------------------- plot --------------------
hFig = figure('Color','w');
imagesc(time, 1:S, gamma);

axis tight;
set(gca, 'YDir', 'normal');

xlabel('Time (s)', 'FontSize', 14, 'FontWeight','bold');
ylabel('States',   'FontSize', 14, 'FontWeight','bold');

set(gca, ...
    'YTick', 1:S, ...
    'YTickLabel', arrayfun(@num2str, 1:S, 'UniformOutput', false), ...
    'FontSize', 12, ...
    'FontWeight','bold');

colormap(parula);
colorbar;

title(sprintf('HMM state occupancy (trial %d)', trialId), ...
    'FontSize', 16, 'FontWeight','bold');

end
