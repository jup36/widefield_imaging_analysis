function Y = motifClusterUmap(Ws, varargin)
% Run UMAP on TEMPORALLY-ALIGNED motifs for low-dimensional embedding.
%
% Inputs:
%   Ws – [P × nMotifs × L] matrix
%        P = valid pixels
%        nMotifs = number of discovered motifs
%        L = time lags
%
% Optional:
%   'n_neighbors' – neighborhood size, default 30
%   'min_dist'    – minimum embedded distance, default 0.1
%   'metric'      – distance metric, default 'correlation'
%
% Output:
%   Y – [nMotifs × 2] UMAP embedding

opts.n_neighbors = 30;
opts.min_dist    = 0.1;
opts.metric      = 'correlation';

opts = ParseOptionalInputs(opts, varargin);

nMotifs = size(Ws, 2);
[P, ~, L] = size(Ws);

% Reshape motifs into [nMotifs × (P × L)]
Ws_vec = reshape(permute(Ws, [2, 1, 3]), nMotifs, P * L);

rng(100);

[Y, ~, ~] = run_umap(Ws_vec, ...
    'n_components', 2, ...
    'n_neighbors', opts.n_neighbors, ...
    'min_dist', opts.min_dist, ...
    'metric', opts.metric, ...
    'verbose', 'none', ...
    'cluster_output', 'none', ...
    'randomize', true);

end