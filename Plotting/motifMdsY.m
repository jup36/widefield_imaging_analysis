function Y = motifMdsY(tcorr_mat)
% Performs multi-dimensional scaling on the temporal correlation (similarity) matrix 

% Convert to dissimilarity matrix
D = 1 - tcorr_mat;
D(isnan(D)) = 1;                     % Treat NaN as max dissimilarity
D = min(max(D, 0), 2);               % Clip to [0, 2]
D(1:size(D,1)+1:end) = 0;            % Force diagonal to zero
D = (D + D') / 2;                    % Ensure symmetry

% Classical MDS
Y = mdscale(D, 2, 'Criterion', 'metricstress');
end
