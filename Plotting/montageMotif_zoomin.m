function montageMotif_zoomin(motif_w, nanpxs, lags, numFrames, varargin)
% motif_w: pixels-by-timelag (expected)

if mod(numFrames, 2)==0
    numFrames = numFrames+1; 
end
framesAdjacent = (numFrames-1)/2; 

if ~isempty(varargin)
    motif_Id = varargin{1}; 
    motif_w = motif_w(:, motif_Id, :); 
    lags = lags(motif_Id,1); 
end

nMotif = size(motif_w,2); 
assert(length(lags)==nMotif); 
allMotifs = zeros(64, 64, 1, 5 * nMotif); % 1 channel for grayscale

% Fill the allMotifs array
for i = 1:nMotif
    motif = conditionDffMat(squeeze(motif_w(:,i,:))', nanpxs); % Get 64x64x10 motif
    allMotifs(:,:,:, (i-1)*numFrames + (1:numFrames)) = motif(:,:,lags(i)-framesAdjacent:lags(i)+framesAdjacent); % Stack frames in correct order
end

% Create the figure and visualize all motifs in one montage
figure;
montage(allMotifs, 'Size', [nMotif, numFrames], 'DisplayRange', [0 0.2]); 
colormap hot;
%colorbar;
title('All Motifs');

end