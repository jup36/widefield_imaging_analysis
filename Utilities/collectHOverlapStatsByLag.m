function hOverlapOut = collectHOverlapStatsByLag(filePath_base, mlist, lagNum, chunkSearchStr, varargin)
% collectHOverlapStatsByLag
%
% Collects H-only temporal overlap/sparsity metrics across mice/sessions/chunks
% for a given motif lag condition.
%
% This mirrors collectMotifStatsByLag, but instead of loading stats_train/test,
% it loads h_train from each chunk file and runs computeHOverlapStats.
%
% Inputs:
%   filePath_base
%       Base directory containing preprocessed motif folders.
%
%   mlist
%       Cell array of mouse IDs, e.g. {'m1044', 'm1045'}.
%
%   lagNum
%       Motif lag condition.
%
%   chunkSearchStr
%       Search pattern for chunk files inside each session folder.
%       Examples:
%           '*green*chunk*.mat'
%           '*red*chunk*.mat'
%
% Optional name-value inputs:
%   'hVarName'
%       Variable name to load from each chunk file. Default = 'h_train'.
%
%   'smoothWin'
%       Common smoothing window for computeHOverlapStats. Default = 19.
%
%   'threshold'
%       Threshold passed to computeHOverlapStats. Default = [].
%
%   'normalizeRows'
%       Whether to row-normalize H before computing overlap. Default = false.
%
%   'useOriginalLag10'
%       If true and lagNum == 10, search for original/default motif folders
%       ending in 'motif' and not containing 'lag'. Default = false.
%
% Output:
%   hOverlapOut
%       Structure containing mouse x session cell arrays for each metric.
%       Each cell contains the mean across chunks within that session.

%% Parse inputs
p = inputParser;
p.addParameter('hVarName', 'h_train', @(x) ischar(x) || isstring(x));
p.addParameter('smoothWin', 19, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('threshold', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('normalizeRows', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('useOriginalLag10', false, @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});

hVarName         = char(p.Results.hVarName);
smoothWin        = p.Results.smoothWin;
threshold        = p.Results.threshold;
normalizeRows    = logical(p.Results.normalizeRows);
useOriginalLag10 = logical(p.Results.useOriginalLag10);

%% Define metrics to collect
metricNames = { ...
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
    'H_pairwise_corr_median'};

nMice = numel(mlist);
maxSessions = 20;

%% Initialize output structure
hOverlapOut = struct();

hOverlapOut.mlist = mlist;
hOverlapOut.lagNum = lagNum;
hOverlapOut.chunkSearchStr = chunkSearchStr;
hOverlapOut.hVarName = hVarName;
hOverlapOut.smoothWin = smoothWin;
hOverlapOut.threshold = threshold;
hOverlapOut.normalizeRows = normalizeRows;
hOverlapOut.useOriginalLag10 = useOriginalLag10;

hOverlapOut.sessionPathC = cell(nMice, maxSessions);
hOverlapOut.chunkPathC   = cell(nMice, maxSessions);
hOverlapOut.chunkStatsC  = cell(nMice, maxSessions);

for m = 1:numel(metricNames)
    hOverlapOut.([metricNames{m}, 'C']) = cell(nMice, maxSessions);
end

%% Main loop
for f = 1:nMice

    %% Find motif session folders
    if lagNum == 10 && useOriginalLag10

        % Original/default motifs:
        % folder name should end with 'motif' and should not contain 'lag'
        sessions = GrabFiles_sort_trials([mlist{f} '*motif'], 0, {filePath_base});

        [~, sessionNames] = cellfun(@fileparts, sessions, 'UniformOutput', false);

        keepI = endsWith(sessionNames, 'motif') & ...
                ~contains(sessionNames, 'lag', 'IgnoreCase', true);

        sessions = sessions(keepI);

    else

        % Explicit lag folders
        sessions = GrabFiles_sort_trials( ...
            [mlist{f} '*motif*lag' num2str(lagNum)], ...
            0, {filePath_base});

    end

    fprintf('\nMouse %s | Lag %d | Search: %s | %d sessions found\n', ...
        mlist{f}, lagNum, chunkSearchStr, numel(sessions));

    %% Iterate through sessions
    for ff = 1:numel(sessions)

        % Initialize as NaN by default
        for m = 1:numel(metricNames)
            hOverlapOut.([metricNames{m}, 'C']){f, ff} = NaN;
        end

        hOverlapOut.sessionPathC{f, ff} = sessions{ff};

        header = extract_date_animalID_header(sessions{ff});

        % Example:
        %   header = 'm1045_122424'
        %   chunkSearchStr = '*green*chunk*.mat'
        %   searchStr = 'm1045_122424*green*chunk*.mat'
        searchStr = [header, chunkSearchStr];

        filePath_chunk = GrabFiles_sort_trials(searchStr, 0, sessions(ff));
        hOverlapOut.chunkPathC{f, ff} = filePath_chunk;

        if isempty(filePath_chunk)
            fprintf('%s | lag %d | session %d/%d: no chunks found with %s\n', ...
                mlist{f}, lagNum, ff, numel(sessions), searchStr);
            continue
        end

        %% Run computeHOverlapStats for each chunk
        metricMat = nan(numel(filePath_chunk), numel(metricNames));
        chunkStats = cell(numel(filePath_chunk), 1);

        for j = 1:numel(filePath_chunk)

            try
                temp = load(filePath_chunk{j}, hVarName);
            catch ME
                warning('Could not load %s from chunk file:\n%s\nError: %s', ...
                    hVarName, filePath_chunk{j}, ME.message);
                continue
            end

            if ~isfield(temp, hVarName)
                warning('Variable %s not found in chunk file:\n%s', ...
                    hVarName, filePath_chunk{j});
                continue
            end

            h_train = temp.(hVarName);

            if isempty(h_train) || ~isnumeric(h_train) || ~ismatrix(h_train)
                warning('Invalid %s in chunk file:\n%s', ...
                    hVarName, filePath_chunk{j});
                continue
            end

            try
                hStats = computeHOverlapStats(h_train, ...
                    'smoothWin', smoothWin, ...
                    'threshold', threshold, ...
                    'normalizeRows', normalizeRows);

                chunkStats{j} = hStats;

                for m = 1:numel(metricNames)
                    if isfield(hStats, metricNames{m})
                        metricMat(j, m) = hStats.(metricNames{m});
                    end
                end

            catch ME
                warning('computeHOverlapStats failed for chunk file:\n%s\nError: %s', ...
                    filePath_chunk{j}, ME.message);
                continue
            end
        end

        hOverlapOut.chunkStatsC{f, ff} = chunkStats;

        %% Average across chunks within this session
        for m = 1:numel(metricNames)
            hOverlapOut.([metricNames{m}, 'C']){f, ff} = ...
                mean(metricMat(:, m), 'omitnan');
        end

        fprintf(['%s | lag %d | session %d/%d | %d chunks | ', ...
            'H overlap norm %.4f | multi-active frac %.4f | active count mean %.3f | Hoyer %.4f\n'], ...
            mlist{f}, lagNum, ff, numel(sessions), numel(filePath_chunk), ...
            hOverlapOut.H_overlap_normC{f, ff}, ...
            hOverlapOut.H_multi_active_fractionC{f, ff}, ...
            hOverlapOut.H_active_count_meanC{f, ff}, ...
            hOverlapOut.H_hoyer_sparsityC{f, ff});

    end
end

end

function hStats = computeHOverlapStats(h_train, varargin)
% computeHOverlapStats
%
% Quantifies temporal sparsity and cross-motif activation overlap in H.
%
% Input:
%   h_train : K x T matrix
%       K motifs by T frames.
%
% Optional name-value inputs:
%   'smoothWin' : smoothing window length in frames, default = 19
%   'threshold' : threshold for defining active H values, default = []
%                 If empty, uses 1e-6 * max(h_train(:)).
%   'normalizeRows' : whether to normalize each motif's H by its max before
%                     computing overlap, default = false.
%
% Output:
%   hStats : structure with H sparsity and overlap metrics.

% -------------------------------------------------------------------------
% Parse inputs
% -------------------------------------------------------------------------
p = inputParser;
p.addParameter('smoothWin', 19, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('threshold', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('normalizeRows', false, @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});

smoothWin = p.Results.smoothWin;
threshold = p.Results.threshold;
normalizeRows = logical(p.Results.normalizeRows);

% -------------------------------------------------------------------------
% Basic checks
% -------------------------------------------------------------------------
assert(isnumeric(h_train) && ismatrix(h_train), ...
    'h_train must be a numeric K x T matrix.');

H = h_train;
[K, T] = size(H);

% Avoid negative values if tiny numerical negatives exist
H(H < 0) = 0;

% -------------------------------------------------------------------------
% Optional row normalization
% -------------------------------------------------------------------------
if normalizeRows
    rowMax = max(H, [], 2);
    rowMax(rowMax == 0) = 1;
    H = H ./ rowMax;
end

% -------------------------------------------------------------------------
% Define active motifs and active time points
% -------------------------------------------------------------------------
if isempty(threshold)
    threshold = 1e-6 * max(H(:));
end

H_binary = H > threshold;

activeMotif = sum(H, 2) > 0;
K_active = sum(activeMotif);

% -------------------------------------------------------------------------
% Smooth H with a common post hoc kernel
% -------------------------------------------------------------------------
kernel = ones(1, smoothWin);
kernel = kernel ./ sum(kernel);

H_smooth = conv2(H, kernel, 'same');

% -------------------------------------------------------------------------
% Cross-motif H overlap
% -------------------------------------------------------------------------
offDiag = ~eye(K);

offActivation = offDiag * H_smooth;

H_overlap_raw = sum(H .* offActivation, 'all');

H_self = sum(H .* H_smooth, 'all');

H_overlap_total = H_self + H_overlap_raw;

H_overlap_norm = H_overlap_raw ./ (H_overlap_total + eps);

if K > 1
    H_overlap_raw_perOther = H_overlap_raw ./ (K - 1);
else
    H_overlap_raw_perOther = 0;
end

% -------------------------------------------------------------------------
% Sparsity metrics
% -------------------------------------------------------------------------

% Fraction of all H entries that are above threshold
H_density = nnz(H_binary) ./ numel(H_binary);
H_sparsity_fraction_zero = 1 - H_density;

% Hoyer sparsity, computed on vectorized H.
% 0 = dense/equal values, 1 = maximally sparse.
hVec = H(:);
n = numel(hVec);

if norm(hVec, 2) > 0
    H_hoyer_sparsity = (sqrt(n) - norm(hVec, 1) / norm(hVec, 2)) / ...
                       (sqrt(n) - 1 + eps);
else
    H_hoyer_sparsity = NaN;
end

% -------------------------------------------------------------------------
% Active motif count per frame
% -------------------------------------------------------------------------
activeCount = sum(H_binary, 1);

H_active_count_mean = mean(activeCount);
H_active_count_median = median(activeCount);
H_active_count_max = max(activeCount);

H_multi_active_fraction = mean(activeCount > 1);
H_any_active_fraction = mean(activeCount > 0);

% -------------------------------------------------------------------------
% Optional pairwise correlation / similarity of H rows
% -------------------------------------------------------------------------
if K > 1
    H_corr = corr(H');
    H_corr(1:K+1:end) = NaN;
    H_pairwise_corr_mean = mean(H_corr(:), 'omitnan');
    H_pairwise_corr_median = median(H_corr(:), 'omitnan');
else
    H_corr = NaN;
    H_pairwise_corr_mean = NaN;
    H_pairwise_corr_median = NaN;
end

% -------------------------------------------------------------------------
% Store outputs
% -------------------------------------------------------------------------
hStats = struct();

hStats.K = K;
hStats.K_active = K_active;
hStats.T = T;

hStats.smoothWin = smoothWin;
hStats.threshold = threshold;
hStats.normalizeRows = normalizeRows;

hStats.H_overlap_raw = H_overlap_raw;
hStats.H_overlap_norm = H_overlap_norm;
hStats.H_overlap_self = H_self;
hStats.H_overlap_total = H_overlap_total;
hStats.H_overlap_raw_perOther = H_overlap_raw_perOther;

hStats.H_density = H_density;
hStats.H_sparsity_fraction_zero = H_sparsity_fraction_zero;
hStats.H_hoyer_sparsity = H_hoyer_sparsity;

hStats.H_active_count_mean = H_active_count_mean;
hStats.H_active_count_median = H_active_count_median;
hStats.H_active_count_max = H_active_count_max;
hStats.H_multi_active_fraction = H_multi_active_fraction;
hStats.H_any_active_fraction = H_any_active_fraction;

hStats.H_pairwise_corr_mean = H_pairwise_corr_mean;
hStats.H_pairwise_corr_median = H_pairwise_corr_median;
hStats.H_pairwise_corr = H_corr;

hStats.H_active_count_trace = activeCount;

end