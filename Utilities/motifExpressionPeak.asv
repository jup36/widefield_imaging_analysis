function [h_motifs_shift, h_lags] = motifExpressionPeak(w_motifs, h_motifs)
% w_motifs: spatiotemporal motifs (N x K x L), i.e., pixels x #motifs x lags 
% h_motifs: temporal weightings (K x T), i.e., #motifs x time (frames)
% This function performs temporal shifting of the temporal weightings
%  matrix (H) fitted by fpCNMF. As can be seen in the code snippet below,
%  H(t) contributes the reconstruction tensor X(:, t+tau-1), shifting the
%  effect forward (later in time) by tau-1 frames. Intuitively, this means
%  that a peak in H(i, :) marks the start of motif expression. Thus,
%  temporally shifting tao frames, at which W peaks, allows us to better
%  interpret H as the expression of a motif at its peak rather than
%  initiation. 

%for tau = 1:L 
%    X = X + W(:, :, tau) * circshift(H,[0,tau-1]); % H(t) contributes X(:, t+tau-1), shifting the effect forward (later in time) by tau-1 frames
%end

[~, h_lags] = max(squeeze(nanmean(w_motifs, 1)), [], 2);

%montageMotif(squeeze(w_motifs(:,12,:)), nanpxs); 

h_motifs_shift = zeros(size(h_motifs)); 

for i = 1:size(h_motifs_shift, 1)
    h_motifs_shift(i, :) = circshift(h_motifs(i, :), [0 h_lags(i)-1]);  
    h_motifs_shift(i, 1:h_lags(i)-1) = 0; 
end


end