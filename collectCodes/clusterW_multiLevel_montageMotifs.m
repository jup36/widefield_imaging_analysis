%This script visualizes basis motifs identified by ClusterW_multiLevels.m
fileSaveDir = '/Volumes/buschman/Rodent Data/dualImaging_parkj/collectData'; 
if ispc
    fileSaveDir = 'Z:\Rodent Data\dualImaging_parkj\collectData'; 
end

%% load the clusterW_output
load(fullfile(fileSaveDir, 'clusterW_output_DAmotifs'), 'W_basisC', 'kval', 'ovr_q', 'cluster_idxC', 'idx_knn', 'tcorr_mat', 'lag_mat', 'lags', 'nanpxs', 'shiftC');

%% montage and print
W = W_basisC{1,2}; 

montageMotifsPrint(W, nanpxs, 8, 1)

%% similarity matrix
ci = cluster_idxC{1, 3}; 
Plot_OrderedSimilarityMatrixPrint(tcorr_mat,ci,1);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function montageMotifsPrint(motif_w, nanpxs, motifsPerFig, printLogic, varargin)
%   montageMotifs(motif_w, nanpxs)
%   montageMotifs(motif_w, nanpxs, motif_Id)
%
%   INPUTS:
%     motif_w   – [P × K × L] array of spatiotemporal motifs (pixels × motifs × frames)
%     nanpxs    – [P × 1] logical array indicating NaN pixels to mask in each frame
%     motif_Id  – (optional) vector of motif indices to display (defaults to all)
%
%   Displays separate figures if motifs exceed `motifsPerFig`.

saveFigDir = '/Volumes/buschman/Rodent Data/dualImaging_parkj/collectFigure/montage';
if ispc
    saveFigDir = 'Z:\Rodent Data\dualImaging_parkj\collectFigure\montage';
end

%motifsPerFig = 14;  % Number of motifs per figure

% Optional: select subset of motifs
if ~isempty(varargin)
    motif_Id = varargin{1};
    motif_w = motif_w(:, motif_Id, :);
end

[P, nMotif, nFrames] = size(motif_w);
nFigures = ceil(nMotif / motifsPerFig);

for figIdx = 1:nFigures
    % Determine motif range for this figure
    startIdx = (figIdx - 1) * motifsPerFig + 1;
    endIdx = min(figIdx * motifsPerFig, nMotif);
    motifsThisFig = endIdx - startIdx + 1;

    allMotifs = zeros(64, 64, 1, motifsThisFig * nFrames); % 1 grayscale channel

    for i = 1:motifsThisFig
        motifIdx = startIdx + i - 1;
        if P == 64 * 64
            motif = reshape(squeeze(motif_w(:, motifIdx, :)), 64, 64, []);
        else
            motif = conditionDffMat(squeeze(motif_w(:, motifIdx, :))', nanpxs);
        end
        allMotifs(:,:,:, (i-1)*nFrames + (1:nFrames)) = motif;
    end

    % Show montage
    figure;
    montage(allMotifs, 'Size', [motifsThisFig, nFrames], 'DisplayRange', [0 0.3]);
    colormap magma;
    title(sprintf('Motifs %d–%d', startIdx, endIdx));

    % Print montage
    % Print montage
    if printLogic == 1
        timestampStr = datestr(now, 'mmddyy_HHMMSS');  % Date and time string
        figSaveName = sprintf('motifs_%d–%d_%s', startIdx, endIdx, timestampStr);
        print(fullfile(saveFigDir, figSaveName), '-dpdf', '-painters', '-bestfit')
    end

end
end


