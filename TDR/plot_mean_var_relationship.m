function plot_mean_var_relationship(Y, varargin)
% plot_mean_var_relationship
%   Visualize the mean–variance relationship for each motif’s activity H.
%
% INPUT
%   Y : [T × K] matrix (time/trial bins × motifs)
%
% OPTIONAL
%   'nbins'   : number of mean bins (default 30)
%   'smooth'  : moving window (in samples) if you prefer running stats (default 0 → binned)
%
% EXAMPLE
%   plot_mean_var_relationship(Y, 'nbins', 40);
%
% OUTPUT
%   Makes a tiled plot showing how variance scales with mean per motif.
%
%   If variance rises sharply with mean → try sqrt/log transform.

% -------------------------------------------------------------------------
p = inputParser;
p.addParameter('nbins', 30, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('smooth', 0, @(x)isnumeric(x)&&isscalar(x));
p.parse(varargin{:});
nbins  = p.Results.nbins;
smooth = p.Results.smooth;

[T,K] = size(Y);
figure('Color','w','Position',[200 100 1000 600]);
tiledlayout('flow','TileSpacing','tight','Padding','compact');

for k = 1:K
    y = Y(:,k);
    y = y(isfinite(y));     % remove NaNs
    if isempty(y)
        nexttile; title(sprintf('Motif %d (empty)', k)); continue;
    end

    if smooth > 0
        % running stats
        mu = movmean(y, smooth, 'omitnan');
        va = movvar(y, smooth, 'omitnan');
    else
        % binned stats
        edges = linspace(min(y), max(y), nbins+1);
        [~,~,bin] = histcounts(y, edges);
        mu = accumarray(bin, y, [], @mean, NaN);
        va = accumarray(bin, y, [], @var,  NaN);
    end

    nexttile;
    scatter(mu, va, 12, 'k', 'filled', 'MarkerFaceAlpha', 0.4);
    hold on;

    % Fit a simple power-law (variance ~ mean^p)
    good = isfinite(mu) & isfinite(va);
    if sum(good)>5
        X = [ones(sum(good),1) log(mu(good))];
        b = X \ log(va(good));
        predMu = linspace(min(mu(good)), max(mu(good)), 100);
        predVa = exp(b(1)) * predMu.^b(2);
        plot(predMu, predVa, 'r-', 'LineWidth', 1.5);
        txt = sprintf('p = %.2f', b(2));
        text(0.05,0.9,txt,'Units','normalized','Color','r','FontSize',8);
    end

    xlabel('Mean(H)');
    ylabel('Var(H)');
    title(sprintf('Motif %d', k));
    grid on;
end

sgtitle('Mean–Variance Relationship Across Motifs (H)');
end
