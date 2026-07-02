function [animalLK_T, groupLK_T] = summarize_LK_acrossAnimals(pevSummaryT, varargin)
% summarize_LK_acrossAnimals
%
% Summarize Lag-K combo metrics across animals in an animal-balanced way.
%
% First averages sessions within each mouse x L x K_input.
% Then computes group mean/std/SEM across animals.
%
% Input:
%   pevSummaryT : session-level summary table
%
% Required columns:
%   mouseID, L, K_input
%
% Optional name-value pairs:
%   'metricNames' : cell array of metric columns to summarize.
%                   Default uses common LK-combo metrics if present.
%   'lags'        : lags to include. Default = unique(pevSummaryT.L)
%   'Ks'          : K_input values to include. Default = unique(pevSummaryT.K_input)
%
% Outputs:
%   animalLK_T : one row per mouse x L x K_input
%   groupLK_T  : one row per L x K_input, summarized across animals

%% Parse inputs
p = inputParser;

defaultMetricNames = { ...
    'meanPEV', ...
    'meanK_discovered', ...
    'mean_H_overlap_raw', ...
    'mean_H_overlap_norm', ...
    'mean_H_multi_active_fraction', ...
    'mean_H_active_count_mean', ...
    'mean_H_hoyer_sparsity' ...
    };

addParameter(p, 'metricNames', defaultMetricNames, @(x) iscell(x) || isstring(x));
addParameter(p, 'lags', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'Ks', [], @(x) isempty(x) || isnumeric(x));

parse(p, varargin{:});

metricNames = cellstr(p.Results.metricNames);
lags = p.Results.lags;
Ks = p.Results.Ks;

%% Checks
requiredVars = {'mouseID', 'L', 'K_input'};
missingRequired = setdiff(requiredVars, pevSummaryT.Properties.VariableNames);

if ~isempty(missingRequired)
    error('pevSummaryT is missing required variable(s): %s', strjoin(missingRequired, ', '));
end

% Keep only metrics that exist
metricNames = metricNames(ismember(metricNames, pevSummaryT.Properties.VariableNames));

if isempty(metricNames)
    error('None of the requested metricNames exist in pevSummaryT.');
end

if isempty(lags)
    lags = sort(unique(pevSummaryT.L(:)))';
end

if isempty(Ks)
    Ks = sort(unique(pevSummaryT.K_input(:)))';
end

mouseList = unique(pevSummaryT.mouseID, 'stable');

%% ------------------------------------------------------------------------
%  Animal-level summary: average sessions within each animal x L x K
% -------------------------------------------------------------------------

animalRows = {};

for m = 1:numel(mouseList)

    thisMouse = mouseList(m);

    if iscategorical(pevSummaryT.mouseID)
        mouseRows = pevSummaryT.mouseID == thisMouse;
        mouseLabel = char(thisMouse);
    else
        mouseRows = strcmp(string(pevSummaryT.mouseID), string(thisMouse));
        mouseLabel = char(string(thisMouse));
    end

    thisMouseT = pevSummaryT(mouseRows, :);

    for iL = 1:numel(lags)
        for iK = 1:numel(Ks)

            L = lags(iL);
            K_input = Ks(iK);

            rowI = thisMouseT.L == L & thisMouseT.K_input == K_input;

            nSessions = sum(rowI);

            metricVals = nan(1, numel(metricNames));

            if nSessions > 0
                for mm = 1:numel(metricNames)
                    vals = thisMouseT.(metricNames{mm})(rowI);
                    metricVals(mm) = mean(vals, 'omitnan');
                end
            end

            animalRows(end+1, :) = [ ...
                {mouseLabel, L, K_input, nSessions}, ...
                num2cell(metricVals) ...
                ];
        end
    end
end

animalLK_T = cell2table(animalRows, ...
    'VariableNames', [{'mouseID', 'L', 'K_input', 'nSessions'}, metricNames]);

animalLK_T.mouseID = categorical(animalLK_T.mouseID);

%% ------------------------------------------------------------------------
%  Group-level summary: summarize animal means across animals
% -------------------------------------------------------------------------

groupRows = {};

for iL = 1:numel(lags)
    for iK = 1:numel(Ks)

        L = lags(iL);
        K_input = Ks(iK);

        rowI = animalLK_T.L == L & animalLK_T.K_input == K_input & animalLK_T.nSessions > 0;

        nAnimals = sum(rowI);

        rowOut = {L, K_input, nAnimals};

        for mm = 1:numel(metricNames)

            vals = animalLK_T.(metricNames{mm})(rowI);

            groupMean = mean(vals, 'omitnan');
            groupStd  = std(vals, 'omitnan');
            groupN    = sum(~isnan(vals));

            if groupN > 0
                groupSem = groupStd ./ sqrt(groupN);
            else
                groupSem = NaN;
            end

            rowOut = [rowOut, {groupMean, groupStd, groupSem}]; %#ok<AGROW>
        end

        groupRows(end+1, :) = rowOut;
    end
end

groupVarNames = {'L', 'K_input', 'nAnimals'};

for mm = 1:numel(metricNames)
    groupVarNames = [groupVarNames, ...
        {['groupMean_' metricNames{mm}], ...
         ['groupStd_'  metricNames{mm}], ...
         ['groupSem_'  metricNames{mm}]}]; %#ok<AGROW>
end

groupLK_T = cell2table(groupRows, 'VariableNames', groupVarNames);

end