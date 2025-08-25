function montageMotifs(motif_w, nanpxs, varargin)
%   montageMotifs(motif_w, nanpxs)
%   montageMotifs(motif_w, nanpxs, motif_Id)
%
%   INPUTS:
%     motif_w   – [P × K × L] array of spatiotemporal motifs (pixels × motifs × frames)
%     nanpxs    – [P × 1] logical array indicating NaN pixels to mask in each frame
%     motif_Id  – (optional) vector of motif indices to display (defaults to all)
%
%   Displays separate figures if motifs exceed `motifsPerFig`.

motifsPerFig = 14;  % Number of motifs per figure

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
end
end


% function montageMotifs(motif_w, nanpxs, varargin)
% %   montageMotifs(motif_w, nanpxs)
% %   montageMotifs(motif_w, nanpxs, motif_Id)
% %
% %   INPUTS:
% %     motif_w   – [P × K × L] array of spatiotemporal motifs (pixels × motifs × frames)
% %     nanpxs    – [P × 1] logical array indicating NaN pixels to mask in each frame
% %     motif_Id  – (optional) vector of motif indices to display (defaults to all)
% %
% %   This function reshapes each motif into a 64×64 image per frame,
% %   applies the nan pixel mask, and displays all frames from all motifs
% %   in a single montage using the 'magma' colormap.
% %
% %   Each row in the montage corresponds to a motif, and each column
% %   corresponds to a frame in time.
% %
% %   EXAMPLE:
% %     montageMotifs(W, nanpxs);         % show all motifs
% %     montageMotifs(W, nanpxs, [1 4]);  % show only motifs 1 and 4
% 
% if ~isempty(varargin)
%     motif_Id = varargin{1}; 
%     motif_w = motif_w(:, motif_Id, :); 
% end
% 
% [~, nMotif, nFrames] = size(motif_w); 
% allMotifs = zeros(64, 64, 1, 5 * nMotif); % 1 channel for grayscale
% 
% % Fill the allMotifs array
% for i = 1:nMotif
%     if size(motif_w, 1)==64*64
%         motif = reshape(squeeze(motif_w(:, i, :)), 64, 64, 10);
%     else
%         motif = conditionDffMat(squeeze(motif_w(:,i,:))', nanpxs); % Get 64x64x10 motif        
%     end
%     allMotifs(:,:,:, (i-1)*nFrames+ (1:nFrames)) = motif; % Stack frames in correct order
% end
% 
% % Create the figure and visualize all motifs in one montage
% figure;
% montage(allMotifs, 'Size', [nMotif, nFrames], 'DisplayRange', [0 0.3]); 
% colormap magma;
% %colorbar;
% title('All Motifs');
% 
% end