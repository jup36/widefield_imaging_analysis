function rez = compute_refit_recon_pev_loo_fast(dffC, hC, wC, tbytDat, nanpxs, nanpxsC, varargin)
% compute_refit_recon_pev_loo_fast
%
% Fast computation of full and leave-one-out reconstruction PEV for chunked refit motifs.
% Uses linearity:
%   recon_drop_k = recon_full - recon_single_k
%
% INPUTS
%   dffC    : 1 x nChunks cell, each [X x Y x Tchunk]
%   hC      : 1 x nChunks cell, each [K x Tchunk]
%   wC      : 1 x nChunks cell, each [Pvalid x K x L]
%   tbytDat : trial struct array with field blueLEDPulsesOfTrain
%   nanpxs  : global nan mask (logical or indices)
%   nanpxsC : per-chunk nan mask cell array (logical or indices)
%
% NAME-VALUE
%   'TrialsPerChunk'         : default 10
%   'ImSize'                 : default [64 64]
%   'StoreRecon'             : default true
%   'StoreSingleMotifRecon'  : default false
%
% OUTPUT
%   rez.pevChunkFull   [1 x nChunks]
%   rez.pevChunkLOO    [K x nChunks]
%   rez.pevChunkLoss   [K x nChunks]
%   rez.pevTrialFull   [1 x nTrials]
%   rez.pevTrialLOO    [K x nTrials]
%   rez.pevTrialLoss   [K x nTrials]
%
% Optional stored outputs:
%   rez.reconC
%   rez.reconTrialC
%   rez.dffTrialC
%   rez.singleMotifChunkPEV [K x nChunks]
%
% Notes:
% - pevChunkLoss = pevChunkFull - pevChunkLOO
% - pevTrialLoss = pevTrialFull - pevTrialLOO
% - negative loss is possible and indicates the dropped motif slightly improved fit
%   (e.g. redundancy / overlap / mild overfit)
%
% Assumption:
% - trial-to-chunk mapping is sequential with TrialsPerChunk trials per chunk.
%   If your chunking is irregular, replace that logic with chunk metadata.

p = inputParser;
p.addParameter('TrialsPerChunk', 10, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('ImSize', [64 64], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('StoreRecon', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('StoreSingleMotifRecon', false, @(x)islogical(x)&&isscalar(x));
p.parse(varargin{:});

trialsPerChunk = p.Results.TrialsPerChunk;
imSize = p.Results.ImSize;
storeRecon = p.Results.StoreRecon;
storeSingleMotifRecon = p.Results.StoreSingleMotifRecon;

nChunks = numel(dffC);
nTrials = numel(tbytDat);

if numel(hC) ~= nChunks || numel(wC) ~= nChunks
    error('dffC, hC, and wC must have the same number of chunks.');
end

% infer K from first non-empty hC
K = [];
for i = 1:nChunks
    if ~isempty(hC{1,i})
        K = size(hC{1,i}, 1);
        break
    end
end
if isempty(K)
    error('Could not infer K from hC.');
end

% outputs
pevChunkFull = nan(1, nChunks);
pevChunkLOO  = nan(K, nChunks);
pevChunkLoss = nan(K, nChunks);

pevTrialFull = nan(1, nTrials);
pevTrialLOO  = nan(K, nTrials);
pevTrialLoss = nan(K, nTrials);

validChunk = false(1, nChunks);
validTrial = false(1, nTrials);

if storeRecon
    reconC = cell(1, nChunks);
    reconTrialC = cell(1, nTrials);
    dffTrialC = cell(1, nTrials);
else
    reconC = [];
    reconTrialC = [];
    dffTrialC = [];
end

if storeSingleMotifRecon
    singleMotifChunkPEV = nan(K, nChunks);
else
    singleMotifChunkPEV = [];
end

for i = 1:nChunks
    if isempty(dffC{1,i}) || isempty(hC{1,i}) || isempty(wC{1,i})
        continue
    end

    hCh = hC{1,i};   % [K x T]
    wCh = wC{1,i};   % [Pvalid x K x L]

    if size(hCh,1) ~= K
        error('Chunk %d has K=%d, expected K=%d.', i, size(hCh,1), K);
    end
    if size(wCh,2) ~= K
        error('Chunk %d has size(wCh,2)=%d, expected K=%d.', i, size(wCh,2), K);
    end

    % -------- full reconstruction in valid-pixel space
    whFull_valid = tensor_convolve(wCh, hCh);   % [Pvalid x T]
    reconFull = i_reconstruct_full(whFull_valid, nanpxs, nanpxsC, i, imSize);

    if isempty(reconFull)
        continue
    end

    dffChunk = dffC{1,i};
    if ~isequal(size(dffChunk), size(reconFull))
        warning('Chunk %d: reconstructed chunk size does not match dffChunk. Skipping chunk.', i);
        continue
    end

    validChunk(i) = true;

    if storeRecon
        reconC{1,i} = reconFull;
    end

    % -------- full chunk PEV
    pevChunkFull(i) = i_compute_pev(dffChunk, reconFull);

    % -------- precompute all single-motif reconstructions in full image space
    single_full = cell(1, K);
    for k = 1:K
        wk = wCh(:,k,:);   % [Pvalid x 1 x L]
        hk = hCh(k,:);     % [1 x T]
        single_valid_k = tensor_convolve(wk, hk);    % [Pvalid x T]
        single_full{1,k} = i_reconstruct_full(single_valid_k, nanpxs, nanpxsC, i, imSize);
    end

    % -------- chunk-level leave-one-out PEV
    for k = 1:K
        reconSingle = single_full{1,k};
        if isempty(reconSingle)
            continue
        end
        if ~isequal(size(reconSingle), size(reconFull))
            continue
        end

        reconDrop = reconFull - reconSingle;

        pevChunkLOO(k,i) = i_compute_pev(dffChunk, reconDrop);
        pevChunkLoss(k,i) = pevChunkFull(i) - pevChunkLOO(k,i);

        if storeSingleMotifRecon
            singleMotifChunkPEV(k,i) = i_compute_pev(dffChunk, reconSingle);
        end
    end

    % -------- trial-level stats
    trStart = (i-1)*trialsPerChunk + 1;
    trEnd   = min(i*trialsPerChunk, nTrials);

    for tr = trStart:trEnd
        if ~isfield(tbytDat(tr), 'blueLEDPulsesOfTrain')
            continue
        end

        frameEdge = tbytDat(tr).blueLEDPulsesOfTrain;
        if isempty(frameEdge)
            continue
        end
        if iscell(frameEdge)
            frameEdge = cell2mat(frameEdge);
        end
        if isempty(frameEdge)
            continue
        end

        frameEdge = frameEdge(:)';
        frameStart = frameEdge(1);
        frameStop = frameEdge(end);

        if ~isfinite(frameStart) || ~isfinite(frameStop)
            continue
        end

        frameStart = max(1, round(frameStart));
        frameStop  = min(size(reconFull,3), round(frameStop));

        if frameStop < frameStart
            continue
        end

        frameIdx = frameStart:frameStop;

        reconTr = reconFull(:,:,frameIdx);
        dffTr = dffChunk(:,:,frameIdx);

        if ~isequal(size(reconTr), size(dffTr))
            continue
        end

        pevTrialFull(tr) = i_compute_pev(dffTr, reconTr);
        validTrial(tr) = true;

        if storeRecon
            reconTrialC{tr} = reconTr;
            dffTrialC{tr} = dffTr;
        end

        for k = 1:K
            reconSingle = single_full{1,k};
            if isempty(reconSingle)
                continue
            end
            if size(reconSingle,3) < frameStop
                continue
            end

            reconSingleTr = reconSingle(:,:,frameIdx);
            reconDropTr = reconTr - reconSingleTr;

            pevTrialLOO(k,tr) = i_compute_pev(dffTr, reconDropTr);
            pevTrialLoss(k,tr) = pevTrialFull(tr) - pevTrialLOO(k,tr);
        end
    end
end

% -------- package output
rez = struct();
rez.pevChunkFull = pevChunkFull;
rez.pevChunkLOO  = pevChunkLOO;
rez.pevChunkLoss = pevChunkLoss;

rez.pevTrialFull = pevTrialFull;
rez.pevTrialLOO  = pevTrialLOO;
rez.pevTrialLoss = pevTrialLoss;

rez.validChunk = validChunk;
rez.validTrial = validTrial;

rez.K = K;
rez.nChunks = nChunks;
rez.nTrials = nTrials;
rez.trialsPerChunk = trialsPerChunk;
rez.imSize = imSize;

if storeRecon
    rez.reconC = reconC;
    rez.reconTrialC = reconTrialC;
    rez.dffTrialC = dffTrialC;
end

if storeSingleMotifRecon
    rez.singleMotifChunkPEV = singleMotifChunkPEV;
end

end


function reconFull = i_reconstruct_full(whCh, nanpxs, nanpxsC, iChunk, imSize)
% Convert valid-pixel x time reconstruction back to full image x time stack.
% Logic:
%   1) try chunk-specific mask nanpxsC{iChunk} if pixel count matches
%   2) otherwise try global nanpxs if pixel count matches
%   3) otherwise return []

reconFull = [];
nPixFull = prod(imSize);

% ---- first try chunk-specific mask
if ~isempty(nanpxsC) && numel(nanpxsC) >= iChunk && ~isempty(nanpxsC{iChunk})
    maskCh = nanpxsC{iChunk};

    if islogical(maskCh)
        nBad = nnz(maskCh);
        maskIdx = find(maskCh);
    else
        nBad = numel(maskCh);
        maskIdx = maskCh(:);
    end

    if size(whCh,1) == (nPixFull - nBad)
        reconFull = conditionDffMat(whCh', maskIdx);
        return
    end
end

% ---- fallback to global mask
if ~isempty(nanpxs)
    if islogical(nanpxs)
        nBad = nnz(nanpxs);
        maskIdx = find(nanpxs);
    else
        nBad = numel(nanpxs);
        maskIdx = nanpxs(:);
    end

    if size(whCh,1) == (nPixFull - nBad)
        reconFull = conditionDffMat(whCh', maskIdx);
        return
    end
end

end


function pev = i_compute_pev(y, yhat)
% Fraction of variance explained:
%   1 - SSE / SST

pev = nan;

if isempty(y) || isempty(yhat)
    return
end

if ~isequal(size(y), size(yhat))
    return
end

y = double(y);
yhat = double(yhat);

valid = isfinite(y) & isfinite(yhat);
if ~any(valid(:))
    return
end

yv = y(valid);
yhv = yhat(valid);

ss_res = sum((yv - yhv).^2);
ss_tot = sum((yv - mean(yv)).^2);

if ss_tot <= 0
    return
end

pev = 1 - ss_res / ss_tot;

end