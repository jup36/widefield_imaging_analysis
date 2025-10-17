function plotIndividualTimeSeries(datCell, xTime, animalID)
% plotIndividualTimeSeries(datCell, xTime, animalID)
%
% Inputs:
%   datCell  : {nAnimals x 1} cell array, each containing 1xN time series (or empty)
%   xTime    : 1xN numeric vector (e.g., timepoints)
%   animalID : 1xn cell array of strings for legend labels (optional)
%
% Example:
%   plotIndividualTimeSeries(myDat, xTime, mIdC)

% --------- Input checks
nAnimals = numel(datCell);
if nargin < 3 || isempty(animalID)
    animalID = arrayfun(@(i) sprintf('Animal %d', i), 1:nAnimals, 'UniformOutput', false);
end
if length(animalID) ~= nAnimals
    error('Length of animalID must match number of animals in datCell.');
end

% --------- Color assignment
animalCols = lines(max(nAnimals,7));

% --------- Setup figure
% Get the current position
fig = gcf;
pos = get(fig, 'Position');

% Double the width (the 3rd element)
pos(3) = pos(3) * 1.5;

% Set the new position
set(fig, 'Position', pos);
hold on;
box off;
xlabel('Time');
ylabel('PEV (%)');
title('Individual Animal Time Series');

% --------- Plotting
hasData = false(1, nAnimals);
hLegend = gobjects(1, nAnimals);

for a = 1:nAnimals
    d = datCell{a};
    if isempty(d) || all(isnan(d)), continue; end

    plot(xTime, d, '-', ...
        'Color', animalCols(a,:), ...
        'LineWidth', 1.4);

    % Create dummy for legend
    hLegend(a) = plot(nan, nan, '-', ...
        'Color', animalCols(a,:), ...
        'LineWidth', 2);
    
    hasData(a) = true;
end

% --------- Legend (only non-empty animals)
valid = hasData;
legend(hLegend(valid), animalID(valid), ...
    'Location', 'northeastoutside', ...
    'Box', 'off');

hold off; 
end
