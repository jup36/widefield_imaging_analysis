function plotLambdaByAnimalSessions(lambdaC, varargin)
% plotLambdaByAnimalSessions(trainTestRez.lambda, 'animalID', mIdC, ...
%     'Jitter', 0.18, 'AlphaRange', [0.06 0.98], 'PointSize', 28, 'LogY', true)
%
% lambdaC : {nAnimals x nSessions} cell. Each cell contains a vector OR a
%           1xN cell of numeric λ values (e.g., per chunk). Empty cells ok.
%
% Name–Value options:
%   'animalID'   : cellstr of length nAnimals (legend labels). Default: 'Animal i'
%   'Jitter'     : horizontal jitter half-width (default 0.18)
%   'AlphaRange' : [amin amax] session transparency (default [0.06 0.98])
%   'PointSize'  : marker size (default 28)
%   'LogY'       : logical, use semilogy-style y scale (default true)

% ---------- parse inputs
p = inputParser;
p.addParameter('animalID', {}, @(x) iscellstr(x) || (iscell(x) && all(cellfun(@ischar,x))));
p.addParameter('Jitter', 0.18, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('AlphaRange', [0.06 0.98], @(x) isnumeric(x) && numel(x)==2);
p.addParameter('PointSize', 28, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('LogY', true, @(x) islogical(x) || ismember(x,[0 1]));
p.parse(varargin{:});
opt = p.Results;

% ---------- sizes, colors, alpha
nAnimals  = size(lambdaC, 1);
nSessions = size(lambdaC, 2);
animalCols = lines(max(nAnimals,7));
alphaVals  = linspace(opt.AlphaRange(1), opt.AlphaRange(2), max(nSessions,1)).^0.65;

% legend labels
if isempty(opt.animalID)
    labels = arrayfun(@(i) sprintf('Animal %d', i), 1:nAnimals, 'UniformOutput', false);
else
    labels = opt.animalID(:).';
end

figure('Color','w'); hold on;

% Track animals that actually contributed points
hasData = false(1, nAnimals);
yAll = []; % collect all plotted values for ylim later

% ---------- scatter all lambdas
for a = 1:nAnimals
    col = animalCols(a,:);
    for s = 1:nSessions
        v = lambdaC{a,s};
        if isempty(v), continue; end

        % accept either numeric vector or 1xN cell of numerics
        if iscell(v), v = [v{:}]; end
        v = v(:);
        v = v(isfinite(v));          % remove NaN/Inf
        if isempty(v), continue; end
        if opt.LogY
            v = v(v>0);              % log-scale needs positive
            if isempty(v), continue; end
        end

        x  = a * ones(size(v));
        xj = x + (rand(size(x))*2 - 1) * opt.Jitter;

        scatter(xj, v, opt.PointSize, ...
            'MarkerFaceColor', col, ...
            'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', alphaVals(min(s,numel(alphaVals))), ...
            'MarkerEdgeAlpha', min(1, alphaVals(min(s,numel(alphaVals))) + 0.15));
        hasData(a) = true;
        yAll = [yAll; v]; %#ok<AGROW>
    end
end

plot(1:nAnimals, ones(nAnimals).*0.0005, 'r:')

% keep only animals with data
valid = find(hasData);
xticks(1:numel(valid));
xticklabels(labels(valid));
xtickangle(30);

xlim([0.5, numel(valid)+0.5]);
xlabel('Animal');
ylabel('\lambda');
title('\lambda by Animal (color) and Session (alpha)');

if opt.LogY
    set(gca,'YScale','log');
end

% ---------- add y-margin (10%)
if ~isempty(yAll)
    yl = [min(yAll) max(yAll)];
    margin = 0.1 * diff(log10(yl)); % margin in log space
    if opt.LogY
        ylim([10^(log10(yl(1))-margin), 10^(log10(yl(2))+margin)]);
    else
        margin = 0.1 * diff(yl);
        ylim([yl(1)-margin, yl(2)+margin]);
    end
end

set(gca, 'TickDir','out', 'Layer','top', 'Box','off');

% ---------- legend: opaque color chips (only for animals with data)
if ~isempty(valid)
    hLeg = gobjects(numel(valid),1);
    for k = 1:numel(valid)
        a = valid(k);
        hLeg(k) = scatter(nan, nan, opt.PointSize, ...
            'MarkerFaceColor', animalCols(a,:), ...
            'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', 1);
    end
    legend(hLeg, labels(valid), 'Location','northeastoutside');
end
hold off;
end
