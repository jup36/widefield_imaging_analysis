function [hFig, h] = plotMeanWithSemShade(meanC, semC, colorC, varargin)
% plotMeanWithSemShade  
% Plot mean ± SEM curves with translucent shaded error regions.
%
% [hFig, h] = plotMeanWithSemShade(meanC, semC, colorC, ...)
%
% This version automatically creates a new figure (white background).
%
% Inputs
%   meanC   : T×N numeric OR 1×N cell, each element a T×1 vector
%   semC    : T×N numeric OR 1×N cell, same structure as meanC
%   colorC  : 1×N cell array, each {1×3} RGB triplet
%
% Name–Value options
%   'X'         : x-axis vector (default 1:T)
%   'LineWidth' : line width (default 1.5)
%   'Alpha'     : shaded patch alpha (default 0.30)
%
% Outputs
%   hFig : figure handle (new figure created)
%   h    : struct with fields .patch and .line for each series
%

% ---------- parse options ----------
p = inputParser;
p.addParameter('X', [], @(v) isnumeric(v) && isvector(v));
p.addParameter('LineWidth', 1.5, @(v) isnumeric(v) && isscalar(v) && v>0);
p.addParameter('Alpha', 0.30, @(v) isnumeric(v) && isscalar(v) && v>=0 && v<=1);
p.parse(varargin{:});
opt = p.Results;

% ---------- create figure ----------
hFig = figure('Color','white');  % <-- required modification
ax = axes('Parent', hFig); 
hold(ax, 'on');
set(ax, 'Layer','top');

% ---------- coerce inputs to cell columns ----------
[mc, sc] = coerceToCellCols(meanC, semC);

N = numel(mc);
assert(iscell(colorC) && numel(colorC)==N, ...
    'colorC must be a 1×N cell array of RGB triplets.');

for i = 1:N
    assert(isnumeric(colorC{i}) && numel(colorC{i})==3, ...
        'Each entry of colorC must be a 1×3 RGB vector.');
end

T = numel(mc{1});
if isempty(opt.X)
    x = (1:T).';
else
    x = opt.X(:);
    assert(numel(x)==T, 'Length of X must match lengths in meanC.');
end

h.patch = gobjects(1,N);
h.line  = gobjects(1,N);

% ---------- plot each series ----------
for i = 1:N
    m = mc{i}(:);
    s = sc{i}(:);
    c = colorC{i};

    finiteMask = isfinite(m) & isfinite(s) & isfinite(x);
    segs = maskToSegments(finiteMask);

    for k = 1:size(segs,1)
        idx = segs(k,1):segs(k,2);

        xi = x(idx);
        ui = m(idx) + s(idx);
        li = m(idx) - s(idx);

        % shaded region
        h.patch(i) = patch(ax, ...
            [xi; flipud(xi)], [ui; flipud(li)], c, ...
            'EdgeColor','none', 'FaceAlpha', opt.Alpha);

        % mean line
        h.line(i) = plot(ax, xi, m(idx), '-', ...
            'Color', c, 'LineWidth', opt.LineWidth);
    end
end

xlabel(ax, 'Time');
ylabel(ax, 'Mean ± SEM');
box(ax, 'off');
grid on; 
set(gca, 'TickDir', 'out')

end % main

% ---------- helpers ----------
function [mc, sc] = coerceToCellCols(m, s)
    if isnumeric(m)
        [T,N] = size(m);
        assert(all(size(s)==[T N]), ...
            'meanC and semC must be same size.');
        mc = arrayfun(@(j) m(:,j), 1:N, 'UniformOutput', false);
        sc = arrayfun(@(j) s(:,j), 1:N, 'UniformOutput', false);
    else
        assert(iscell(m) && iscell(s) && numel(m)==numel(s));
        N = numel(m);
        T = numel(m{1});
        for j = 1:N
            assert(numel(m{j})==T && numel(s{j})==T, ...
                'All cell elements must have same length.');
        end
        mc = m; sc = s;
    end
end

function segs = maskToSegments(mask)
    d = diff([false; mask(:); false]);
    starts = find(d==1);
    ends   = find(d==-1) - 1;
    segs = [starts ends];
end
