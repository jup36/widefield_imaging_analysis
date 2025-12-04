function plotBasisKernels(basis)
% plotBasisKernels  Visualize raised-cosine temporal basis functions
%
%   plotBasisKernels(basis)
%
%   Input:
%       basis – structure with fields:
%           • B        : [nTime x nBases] matrix of basis functions
%           • lagsSec  : [nTime x 1] vector of time lags (in seconds)
%
%   The function displays the basis set using imagesc:
%       • X-axis = time (s)
%       • Y-axis = basis kernel index (bottom-up)
%       • White background
%
%   Example:
%       plotBasisKernels(basisToneOn)

    % --- Validate input
    if ~isstruct(basis) || ~isfield(basis,'B') || ~isfield(basis,'lagsSec')
        error('Input must be a struct with fields "B" and "lagsSec".');
    end

    B = basis.B;
    t = basis.lagsSec(:);

    if size(B,1) ~= numel(t)
        error('Number of time points in B must match length of lagsSec.');
    end

    % --- Create figure
    figure('Color','w');
    imagesc(t, 1:size(B,2), B.');   % transpose but DO NOT flip vertically
    axis tight;
    colormap(parula);
    set(gca,'YDir','normal','Box','off','TickDir','out');

    % --- Labels & aesthetics
    xlabel('Time (s)','FontSize',12);
    ylabel('Basis kernel','FontSize',12);
    title('Raised-cosine temporal basis','FontWeight','normal','FontSize',12);
    colorbar;

end
