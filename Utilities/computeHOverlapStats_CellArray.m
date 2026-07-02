function [hSummaryT, hChunkT] = computeHOverlapStats_CellArray(hCell, varargin)
% computeHOverlapStats_CellArray
%
% Computes H-overlap / H-sparsity metrics for a cell array of H matrices.
%
% Each cell should contain one chunk's H:
%   H = K x T
%
% Example:
%   hM = hC(2, :);   % dynamic / motif H chunks
%   hS = hC(1, :);   % static H chunks
%
%   [hSummaryM, hChunkM] = computeHOverlapStats_CellArray(hM);
%   [hSummaryS, hChunkS] = computeHOverlapStats_CellArray(hS);
%
% Required helper:
%   computeHOverlapStats.m
%
% Outputs:
%   hSummaryT : 1-row table containing mean/std/SEM across chunks
%   hChunkT   : chunk-level table

%% Parse inputs
p = inputParser;

addParameter(p, 'smoothWin', 9, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'threshold', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'normalizeRows', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'label', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'verbose', true, @(x) islogical(x) || isnumeric(x));

parse(p, varargin{:});

smoothWin = p.Results.smoothWin;
threshold = p.Results.threshold;
normalizeRows = logical(p.Results.normalizeRows);
label = char(p.Results.label);
verbose = logical(p.Results.verbose);

%% Metrics to collect
hMetricNames = { ...
    'K', ...
    'K_active', ...
    'T', ...
    'H_overlap_raw', ...
    'H_overlap_norm', ...
    'H_overlap_self', ...
    'H_overlap_total', ...
    'H_overlap_raw_perOther', ...
    'H_density', ...
    'H_sparsity_fraction_zero', ...
    'H_hoyer_sparsity', ...
    'H_active_count_mean', ...
    'H_active_count_median', ...
    'H_active_count_max', ...
    'H_multi_active_fraction', ...
    'H_any_active_fraction', ...
    'H_pairwise_corr_mean', ...
    'H_pairwise_corr_median' ...
    };

hMetricMeanNames = strcat('mean_', hMetricNames);
hMetricStdNames  = strcat('std_',  hMetricNames);
hMetricSemNames  = strcat('sem_',  hMetricNames);

%% Basic checks
if ~iscell(hCell)
    error('hCell must be a cell array where each cell contains a K x T H matrix.');
end

hCell = hCell(:);  % force column vector
nChunksTotal = numel(hCell);

hMetricMat = nan(nChunksTotal, numel(hMetricNames));
validChunk = false(nChunksTotal, 1);

%% Loop over chunks
for c = 1:nChunksTotal

    H = hCell{c};

    if isempty(H)
        if verbose
            fprintf('  Chunk %d/%d: empty H. Skipping.\n', c, nChunksTotal);
        end
        continue
    end

    if ~isnumeric(H) || ~ismatrix(H)
        warning('Chunk %d/%d: H is not a numeric matrix. Skipping.', c, nChunksTotal);
        continue
    end

    % Defensive orientation check:
    % H should be K x T. Usually K << T.
    % If rows look much larger than columns, likely transposed.
    if size(H, 1) > size(H, 2)
        H = H';
        if verbose
            fprintf('  Chunk %d/%d: transposed H to K x T orientation.\n', c, nChunksTotal);
        end
    end

    try
        hStats = computeHOverlapStats(H, ...
            'smoothWin', smoothWin, ...
            'threshold', threshold, ...
            'normalizeRows', normalizeRows);

        for hm = 1:numel(hMetricNames)
            if isfield(hStats, hMetricNames{hm})
                hMetricMat(c, hm) = hStats.(hMetricNames{hm});
            end
        end

        validChunk(c) = true;

        if verbose
            fprintf(['  Chunk %d/%d: K=%d, T=%d, raw overlap=%.4g, ' ...
                     'norm overlap=%.4f, multi-active=%.4f\n'], ...
                c, nChunksTotal, ...
                hStats.K, ...
                hStats.T, ...
                hStats.H_overlap_raw, ...
                hStats.H_overlap_norm, ...
                hStats.H_multi_active_fraction);
        end

    catch ME
        warning('computeHOverlapStats failed for chunk %d/%d. Error: %s', ...
            c, nChunksTotal, ME.message);
    end
end

%% Chunk-level table
chunkID = (1:nChunksTotal)';
conditionLabel = repmat({label}, nChunksTotal, 1);

hChunkT = array2table(hMetricMat, 'VariableNames', hMetricNames);
hChunkT = [table(conditionLabel, chunkID, validChunk, ...
    'VariableNames', {'label', 'chunkID', 'validChunk'}), hChunkT];

%% Summary table across chunks
hMetricMean = mean(hMetricMat, 1, 'omitnan');
hMetricStd  = std(hMetricMat, 0, 1, 'omitnan');
hMetricN    = sum(~isnan(hMetricMat), 1);
hMetricSem  = hMetricStd ./ sqrt(hMetricN);
hMetricSem(hMetricN == 0) = NaN;

summaryCell = [ ...
    {label, nChunksTotal, sum(validChunk)}, ...
    num2cell(hMetricMean), ...
    num2cell(hMetricStd), ...
    num2cell(hMetricSem) ...
    ];

hSummaryT = cell2table(summaryCell, ...
    'VariableNames', [ ...
    {'label', 'nChunksTotal', 'nChunksValid'}, ...
    hMetricMeanNames, ...
    hMetricStdNames, ...
    hMetricSemNames ...
    ]);

if verbose
    fprintf('\nSummary for %s:\n', label);
    fprintf('  Valid chunks: %d/%d\n', sum(validChunk), nChunksTotal);
    fprintf('  Mean raw H overlap: %.4g\n', hSummaryT.mean_H_overlap_raw);
    fprintf('  Mean norm H overlap: %.4f\n', hSummaryT.mean_H_overlap_norm);
    fprintf('  Mean multi-active fraction: %.4f\n', hSummaryT.mean_H_multi_active_fraction);
    fprintf('  Mean Hoyer sparsity: %.4f\n', hSummaryT.mean_H_hoyer_sparsity);
end

end