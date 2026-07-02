function Y = motifClusterTsne(Ws, varargin)
% Run t-SNE on temporally aligned motifs.
%
% Ws : [P x nMotifs x L]
%
% Each motif is vectorized into one row:
%   Ws_vec : [nMotifs x P*L]
%
% Default distance:
%   correlation
%
% Other useful options:
%   'euclidean'
%   'cosine'
%   'correlation'

p = inputParser;
p.addParameter('perplexity', 30, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('distance', 'correlation', @(x) ischar(x) || isstring(x));
p.addParameter('rngSeed', 100, @(x) isnumeric(x) && isscalar(x));
p.parse(varargin{:});

perplexity = p.Results.perplexity;
distanceMetric = char(p.Results.distance);
rngSeed = p.Results.rngSeed;

nMotifs = size(Ws, 2);
[P, ~, L] = size(Ws);

Ws_vec = reshape(permute(Ws, [2, 1, 3]), nMotifs, P * L);

rng(rngSeed);

Y = tsne(Ws_vec, ...
    'Perplexity', perplexity, ...
    'NumDimensions', 2, ...
    'Distance', distanceMetric);

end


% function Y = motifClusterTsne(Ws) 
% % Run t-SNE on the TEMPORALLY-ALIGNED motifs for low-dimensional embedding
% % of individual motifs on a t-SNE space (2-D)
% % Ws         – [P × nMotifs × L] matrix (temporally aligned motifs)
% 
% perplexity = 30;
% nMotifs = size(Ws, 2);
% [P, ~, L] = size(Ws);
% % Reshape motifs into [nMotifs × (P×L)] vectors
% Ws_vec = reshape(permute(Ws, [2, 1, 3]), nMotifs, P * L);
% 
% rng(100);
% 
% % Compute t-SNE
% Y = tsne(Ws_vec, 'Perplexity', perplexity, 'NumDimensions', 2);
% end