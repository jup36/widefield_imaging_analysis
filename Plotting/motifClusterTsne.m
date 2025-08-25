function Y = motifClusterTsne(Ws) 
% Run t-SNE on the TEMPORALLY-ALIGNED motifs for low-dimensional embedding
% of individual motifs on a t-SNE space (2-D)
% Ws         – [P × nMotifs × L] matrix (temporally aligned motifs)

perplexity = 30;
nMotifs = size(Ws, 2);
[P, ~, L] = size(Ws);
% Reshape motifs into [nMotifs × (P×L)] vectors
Ws_vec = reshape(permute(Ws, [2, 1, 3]), nMotifs, P * L);

rng(100);

% Compute t-SNE
Y = tsne(Ws_vec, 'Perplexity', perplexity, 'NumDimensions', 2);
end