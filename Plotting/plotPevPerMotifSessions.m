function plotPevPerMotifSessions(pevPerMotifC, varargin)
% plotPevPerMotifSessions(pevPerMotifC, 'animalID', mIdC, 'Jitter',0.22, ...
%                         'AlphaRange',[.05 .98], 'PointSize',30)
%
% pevPerMotifC : {nAnimals x nSessions} cell; each entry is 1xN PEV row.
%
% Name–Value options:
%   'animalID'   : cellstr of length nAnimals (labels for legend)
%   'Jitter'     : horizontal jitter half-width (default 0.22)
%   'AlphaRange' : [a_min a_max] session transparency (default [.05 .98])
%   'PointSize'  : marker size (default 30)

% ---------- parse inputs
p = inputParser;
p.addParameter('animalID', {}, @(x) iscellstr(x) || (iscell(x) && all(cellfun(@ischar,x))));
p.addParameter('Jitter', 0.22, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('AlphaRange', [.05 .98], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('PointSize', 15, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('yLim', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
p.parse(varargin{:});
opt = p.Results;

% ---------- sizes & colors
nAnimals  = size(pevPerMotifC,1);
nSessions = size(pevPerMotifC,2);

nonEmpty   = ~cellfun(@isempty, pevPerMotifC);
if any(nonEmpty(:))
    nMotifs = max(cellfun(@numel, pevPerMotifC(nonEmpty)));
else
    warning('No data to plot.'); return
end

animalCols = lines(max(nAnimals,7)); % color map
alphaVals = linspace(opt.AlphaRange(1), opt.AlphaRange(2), max(nSessions,1)).^0.6;

% legend labels
if isempty(opt.animalID)
    labels = arrayfun(@(a) sprintf('Animal %d',a), 1:nAnimals, 'UniformOutput',false);
else
    labels = opt.animalID(:).';
    if numel(labels) ~= nAnimals
        labels = labels(1:min(end,nAnimals));
        if numel(labels) < nAnimals
            labels(end+1:nAnimals) = {''};
        end
    end
end

% ---------- plot
% Get the current position
fig = gcf;
pos = get(fig, 'Position');

% Double the width (the 3rd element)
pos(3) = pos(3) * 1.5;

% Set the new position
set(fig, 'Position', pos);
hold on; 

% track which animals had data
hasData = false(1,nAnimals);

for a = 1:nAnimals
    col = animalCols(a,:);
    for s = 1:nSessions
        v = pevPerMotifC{a,s};
        if isempty(v), continue; end

        x  = 1:numel(v);
        xj = x + (rand(size(x))*2 - 1) * opt.Jitter;

        scatter(xj, v, opt.PointSize, ...
            'MarkerFaceColor', col, ...
            'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', alphaVals(s), ...
            'MarkerEdgeAlpha', min(1, alphaVals(s)+0.15));
        hasData(a) = true;
    end
end

% ---------- axes
xlabel('Motif #');
ylabel('PEV');
title('Motif-wise Percent Explained Variance (sessions overlaid)');
xlim([0.5, nMotifs+0.5]);

if isempty(opt.yLim)
    yl = ylim;
    ylim([min(0, yl(1)), 0.3]);
else
    ylim(opt.yLim);
end
set(gca, 'Layer','top', 'TickDir','out', 'Box','off');

% ---------- legend: create dummy opaque markers
valid = hasData & ~cellfun(@isempty, labels);
if any(valid)
    hLegend = gobjects(sum(valid),1);
    legLabs = labels(valid);
    cols    = animalCols(valid,:);
    for k = 1:sum(valid)
        hLegend(k) = scatter(nan,nan,opt.PointSize, ...
            'MarkerFaceColor', cols(k,:), ...
            'MarkerEdgeColor','none', ...
            'MarkerFaceAlpha',1); % full opacity
    end
    legend(hLegend, legLabs, 'Location','northeastoutside');
end
end
