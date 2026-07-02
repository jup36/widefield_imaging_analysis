function [summaryRows, chunkRows, nChunksValid] = collectChunkSetToRows( ...
    chunkFileC, ...
    summaryRows, ...
    chunkRows, ...
    mouseID, ...
    sessionID, ...
    sourceLabel, ...
    L, ...
    K_input, ...
    hMetricNames, ...
    commonSmoothWin, ...
    hThreshold, ...
    normalizeHRows)
% collectChunkSetToRows
%
% Collects PEV, discovered K, H metrics, and fitted lambda from one set of
% chunk files.
%
% This function is used for both:
%   1) LKcombo chunk folders
%   2) legacy/original motif folders outside LKcombo
%
% Important:
%   stats_train.lambda is collected only for legacy files, because those
%   runs used fitted lambda. For LKcombo rows, lambda_train is stored as NaN.
%
% It appends both chunk-level and summary-level rows.

nChunkFiles = numel(chunkFileC);

pevChunks = nan(nChunkFiles, 1);
KdiscoveredChunks = nan(nChunkFiles, 1);
lambdaTrainChunks = nan(nChunkFiles, 1);
hMetricMat = nan(nChunkFiles, numel(hMetricNames));

nChunksValid = 0;

% Only collect fitted lambda for legacy/original motif results.
collectLambda = contains(lower(string(sourceLabel)), "legacy");

for ii = 1:nChunkFiles
    
    chunkFile = chunkFileC{ii};
    
    try
        loadedData = load(chunkFile, 'stats_test', 'stats_train', 'w', 'h');
        
        %% ----------------------------------------------------------------
        %  PEV
        % -----------------------------------------------------------------
        
        if isfield(loadedData, 'stats_test') && ...
                isfield(loadedData.stats_test, 'pev')
            
            stats_test = loadedData.stats_test;
            pevVal = mean(stats_test.pev(:), 'omitnan');
            
        elseif isfield(loadedData, 'stats_train') && ...
                isfield(loadedData.stats_train, 'pev')
            
            % Fallback only. Prefer stats_test whenever available.
            stats_train = loadedData.stats_train;
            pevVal = mean(stats_train.pev(:), 'omitnan');
            
            warning('Using stats_train.pev fallback for file:\n%s', chunkFile);
            
        else
            warning('Missing stats_test.pev and stats_train.pev in file:\n%s', chunkFile);
            pevVal = NaN;
        end
        
        pevChunks(ii) = pevVal;
        
        %% ----------------------------------------------------------------
        %  Fitted lambda from stats_train.lambda
        % -----------------------------------------------------------------
        
        lambdaVal = NaN;
        
        if collectLambda
            
            if isfield(loadedData, 'stats_train') && ...
                    isfield(loadedData.stats_train, 'lambda') && ...
                    ~isempty(loadedData.stats_train.lambda)
                
                lambdaRaw = loadedData.stats_train.lambda;
                
                if isnumeric(lambdaRaw)
                    % Usually scalar. If vector/array, summarize to one value
                    % for this chunk.
                    lambdaVal = mean(lambdaRaw(:), 'omitnan');
                    
                    if numel(lambdaRaw) > 1
                        warning(['stats_train.lambda has %d values; using mean(lambda(:)) for file:\n%s'], ...
                            numel(lambdaRaw), chunkFile);
                    end
                else
                    warning('stats_train.lambda is non-numeric in file:\n%s', chunkFile);
                end
                
            else
                warning('Missing stats_train.lambda in legacy file:\n%s', chunkFile);
            end
        end
        
        lambdaTrainChunks(ii) = lambdaVal;
        
        %% ----------------------------------------------------------------
        %  Discovered K from w
        % -----------------------------------------------------------------
        
        if isfield(loadedData, 'w') && ~isempty(loadedData.w)
            % w is expected to be [pixels x K_discovered x lags]
            K_discovered = size(loadedData.w, 2);
        else
            warning('Missing w in file: %s', chunkFile);
            K_discovered = NaN;
        end
        
        KdiscoveredChunks(ii) = K_discovered;
        
        %% ----------------------------------------------------------------
        %  H-overlap / H-sparsity metrics
        % -----------------------------------------------------------------
        
        if isfield(loadedData, 'h') && ~isempty(loadedData.h)
            
            H = loadedData.h;
            
            % Expected orientation: K x T.
            % If H appears transposed relative to discovered K, transpose it.
            if ~isnan(K_discovered)
                
                if size(H, 1) ~= K_discovered && size(H, 2) == K_discovered
                    H = H';
                end
                
                if size(H, 1) ~= K_discovered
                    warning(['H dimension does not match K_discovered in file:\n%s\n' ...
                             'size(H) = [%d %d], K_discovered = %d'], ...
                             chunkFile, size(H, 1), size(H, 2), K_discovered);
                end
            end
            
            if isempty(H) || ~isnumeric(H) || ~ismatrix(H)
                
                warning('Invalid h in file: %s', chunkFile);
                
            else
                
                try
                    hStats = computeHOverlapStats(H, ...
                        'smoothWin', commonSmoothWin, ...
                        'threshold', hThreshold, ...
                        'normalizeRows', normalizeHRows);
                    
                    for hm = 1:numel(hMetricNames)
                        if isfield(hStats, hMetricNames{hm})
                            hMetricMat(ii, hm) = hStats.(hMetricNames{hm});
                        end
                    end
                    
                catch ME
                    warning('computeHOverlapStats failed for chunk file:\n%s\nError: %s', ...
                        chunkFile, ME.message);
                end
            end
            
        else
            warning('Missing h in file: %s', chunkFile);
        end
        
        %% ----------------------------------------------------------------
        %  Store chunk-level row
        % -----------------------------------------------------------------
        %
        % IMPORTANT:
        % chunkRows variable names must include lambda_train after
        % K_discovered.
        
        chunkRows(end+1, :) = [ ...
            {mouseID, ...
             sessionID, ...
             sourceLabel, ...
             L, ...
             K_input, ...
             ii, ...
             pevVal, ...
             K_discovered, ...
             lambdaVal}, ...
            num2cell(hMetricMat(ii, :)), ...
            {chunkFile} ...
            ];
        
        if ~isnan(pevVal)
            nChunksValid = nChunksValid + 1;
        end
        
    catch ME
        warning('Failed to load/process chunk file:\n%s\nError: %s', ...
            chunkFile, ME.message);
    end
end

%% ------------------------------------------------------------------------
%  Summarize PEV across chunks
% -------------------------------------------------------------------------

meanPEV = mean(pevChunks, 'omitnan');
stdPEV  = std(pevChunks, 'omitnan');
nChunks = sum(~isnan(pevChunks));

if nChunks > 0
    semPEV = stdPEV ./ sqrt(nChunks);
else
    semPEV = NaN;
end

%% ------------------------------------------------------------------------
%  Summarize discovered K across chunks
% -------------------------------------------------------------------------

meanK_discovered = mean(KdiscoveredChunks, 'omitnan');
stdK_discovered  = std(KdiscoveredChunks, 'omitnan');
nKvalid = sum(~isnan(KdiscoveredChunks));

if nKvalid > 0
    semK_discovered = stdK_discovered ./ sqrt(nKvalid);
else
    semK_discovered = NaN;
end

%% ------------------------------------------------------------------------
%  Summarize fitted lambda across chunks
% -------------------------------------------------------------------------

meanLambda_train = mean(lambdaTrainChunks, 'omitnan');
stdLambda_train  = std(lambdaTrainChunks, 'omitnan');
nLambda_train    = sum(~isnan(lambdaTrainChunks));

if nLambda_train > 0
    semLambda_train = stdLambda_train ./ sqrt(nLambda_train);
else
    semLambda_train = NaN;
end

%% ------------------------------------------------------------------------
%  Summarize H metrics across chunks
% -------------------------------------------------------------------------

hMetricMean = mean(hMetricMat, 1, 'omitnan');
hMetricStd  = std(hMetricMat, 0, 1, 'omitnan');
hMetricN    = sum(~isnan(hMetricMat), 1);
hMetricSem  = hMetricStd ./ sqrt(hMetricN);
hMetricSem(hMetricN == 0) = NaN;

%% ------------------------------------------------------------------------
%  Store summary-level row
% -------------------------------------------------------------------------

key = sprintf('%s_%s_L%d_K%d', sourceLabel, sessionID, L, K_input);

% IMPORTANT:
% summaryRows variable names must include:
%   meanLambda_train, stdLambda_train, semLambda_train, nLambda_train
% after semK_discovered.

summaryRows(end+1, :) = [ ...
    {mouseID, ...
     sessionID, ...
     sourceLabel, ...
     L, ...
     K_input, ...
     meanPEV, ...
     stdPEV, ...
     semPEV, ...
     nChunks, ...
     meanK_discovered, ...
     stdK_discovered, ...
     semK_discovered, ...
     meanLambda_train, ...
     stdLambda_train, ...
     semLambda_train, ...
     nLambda_train}, ...
    num2cell(hMetricMean), ...
    num2cell(hMetricStd), ...
    num2cell(hMetricSem), ...
    {key} ...
    ];

%% ------------------------------------------------------------------------
%  Console report
% -------------------------------------------------------------------------

idxRawOverlap   = strcmp(hMetricNames, 'H_overlap_raw');
idxNormOverlap  = strcmp(hMetricNames, 'H_overlap_norm');
idxMultiActive  = strcmp(hMetricNames, 'H_multi_active_fraction');
idxActiveCount  = strcmp(hMetricNames, 'H_active_count_mean');
idxHoyer        = strcmp(hMetricNames, 'H_hoyer_sparsity');

fprintf(['    [%s] L=%d, K_input=%d: mean PEV = %.4f, std = %.4f, SEM = %.4f, ' ...
    'nChunks = %d, mean K_discovered = %.2f +/- %.2f, ' ...
    'lambda_train = %.4g +/- %.4g, nLambda = %d, ' ...
    'raw H overlap = %.4g, norm H overlap = %.4f, ' ...
    'multi-active frac = %.4f, active count = %.3f, Hoyer = %.4f\n'], ...
    sourceLabel, ...
    L, K_input, ...
    meanPEV, stdPEV, semPEV, nChunks, ...
    meanK_discovered, semK_discovered, ...
    meanLambda_train, semLambda_train, nLambda_train, ...
    hMetricMean(idxRawOverlap), ...
    hMetricMean(idxNormOverlap), ...
    hMetricMean(idxMultiActive), ...
    hMetricMean(idxActiveCount), ...
    hMetricMean(idxHoyer));

end