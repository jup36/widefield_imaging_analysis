function montageMotifs(motif_w, nanpxs, varargin)
%   montageMotifs(motif_w, nanpxs)
%   montageMotifs(motif_w, nanpxs, motif_Id)
%
%   INPUTS:
%     motif_w   – [P × K × L] array of spatiotemporal motifs (pixels × motifs × frames)
%     nanpxs    – [P × 1] logical array indicating NaN pixels to mask in each frame
%     motif_Id  – (optional) vector of motif indices to display (defaults to all)
%
%   This function reshapes each motif into a 64×64 image per frame,
%   applies the nan pixel mask, and displays all frames from all motifs
%   in a single montage using the 'magma' colormap.
%
%   Each row in the montage corresponds to a motif, and each column
%   corresponds to a frame in time.
%
%   EXAMPLE:
%     montageMotifs(W, nanpxs);         % show all motifs
%     montageMotifs(W, nanpxs, [1 4]);  % show only motifs 1 and 4

if ~isempty(varargin)
    motif_Id = varargin{1}; 
    motif_w = motif_w(:, motif_Id, :); 
end

[~, nMotif, nFrames] = size(motif_w); 
allMotifs = zeros(64, 64, 1, 5 * nMotif); % 1 channel for grayscale

% Fill the allMotifs array
for i = 1:nMotif
    motif = conditionDffMat(squeeze(motif_w(:,i,:))', nanpxs); % Get 64x64x10 motif
    allMotifs(:,:,:, (i-1)*nFrames+ (1:nFrames)) = motif; % Stack frames in correct order
end

% Create the figure and visualize all motifs in one montage
figure;
montage(allMotifs, 'Size', [nMotif, nFrames], 'DisplayRange', [0 0.5]); 
colormap magma;
%colorbar;
title('All Motifs');

end