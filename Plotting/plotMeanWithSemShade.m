function [hFig, h] = plotMeanWithSemShade(meanC, semC, colorC, varargin)
% plotMeanWithSemShade  
% Plot mean ± SEM curves with translucent shaded error regions.
% Optionally overlays reference means and handles separate legends.
%
% [hFig, h] = plotMeanWithSemShade(meanC, semC, colorC, ...)
%
% Inputs
%   meanC   : T×N numeric OR 1×N cell, each element a T×1 vector
%   semC    : T×N numeric OR 1×N cell, same structure as meanC
%   colorC  : 1×N cell array, each {1×3} RGB triplet
%
% Name–Value options
%   'X'          : x-axis vector (default 1:T)
%   'LineWidth'  : line width (default 1.5)
%   'Alpha'      : shaded patch alpha (default 0.30)
%   'refMeanDat' : Reference data (same size/format as meanC). 
%   'title'      : String for figure title (Interpreter = none)
%   'legends'    : Cell array of strings (Labels for Main Mean Lines)
%   'refLegends' : Cell array of strings (Labels for Ref Dotted Lines)
%   'figVisible' : Logical scalar (true=visible, false=invisible/headless)
%
% Outputs
%   hFig : figure handle
%   h    : struct with fields .patch, .line, and .refLine
%

% ---------- parse options ----------
p = inputParser;
p.addParameter('X', [], @(v) isnumeric(v) && isvector(v));
p.addParameter('LineWidth', 1.5, @(v) isnumeric(v) && isscalar(v) && v>0);
p.addParameter('Alpha', 0.30, @(v) isnumeric(v) && isscalar(v) && v>=0 && v<=1);
p.addParameter('refMeanDat', [], @(v) isnumeric(v) || iscell(v));
p.addParameter('title', '', @(v) ischar(v) || isstring(v));
p.addParameter('legends', {}, @(v) iscell(v) || isstring(v));    
p.addParameter('refLegends', {}, @(v) iscell(v) || isstring(v)); 
p.addParameter('figVisible', true, @(v) islogical(v) && isscalar(v)); % <-- NEW OPTION
p.parse(varargin{:});
opt = p.Results;

% ---------- create figure ----------
% Use the figVisible option to set the 'Visible' property
visibility = 'on';
if ~opt.figVisible
    visibility = 'off';
end

hFig = figure('Color','white', 'Visible', visibility); % <-- MODIFIED
ax = axes('Parent', hFig); 
hold(ax, 'on');
set(ax, 'Layer','top');

% ---------- coerce inputs to cell columns ----------
[mc, sc] = coerceToCellCols(meanC, semC);

N = numel(mc);
assert(iscell(colorC) && numel(colorC)==N, ...
    'colorC must be a 1×N cell array of RGB triplets.');

% Validate Ref Data if present
hasRef = ~isempty(opt.refMeanDat);
if hasRef
    rc = coerceRefToCell(opt.refMeanDat, N);
else
    rc = cell(1, N);
end

for i = 1:N
    assert(isnumeric(colorC{i}) && numel(colorC{i})==3, ...
        'Each entry of colorC must be a 1×3 RGB vector.');
    if hasRef
        assert(numel(rc{i}) == numel(mc{i}), ...
            sprintf('Ref data for series %d length mismatch.', i));
    end
end

T = numel(mc{1});
if isempty(opt.X)
    x = (1:T).';
else
    x = opt.X(:);
    assert(numel(x)==T, 'Length of X must match lengths in meanC.');
end

% Initialize output handles
h.patch   = gobjects(1,N);
h.line    = gobjects(1,N);
h.refLine = gobjects(1,N);

% ---------- plot each series ----------
for i = 1:N
    m = mc{i}(:);
    s = sc{i}(:);
    c = colorC{i};

    % --- 1. Plot Main Mean + SEM ---
    finiteMask = isfinite(m) & isfinite(s) & isfinite(x);
    segs = maskToSegments(finiteMask);

    for k = 1:size(segs,1)
        idx = segs(k,1):segs(k,2);
        xi = x(idx);
        ui = m(idx) + s(idx);
        li = m(idx) - s(idx);

        % Shaded region
        h.patch(i) = patch(ax, ...
            [xi; flipud(xi)], [ui; flipud(li)], c, ...
            'EdgeColor','none', 'FaceAlpha', opt.Alpha);
        
        % EXCLUDE SEM FROM LEGEND (Best practice)
        h.patch(i).Annotation.LegendInformation.IconDisplayStyle = 'off';

        % Mean line
        h.line(i) = plot(ax, xi, m(idx), '-', ...
            'Color', c, 'LineWidth', opt.LineWidth);
    end
    
    % --- 2. Plot Ref Data (if exists) ---
    if hasRef
        r = rc{i}(:);
        refMask = isfinite(r) & isfinite(x);
        refSegs = maskToSegments(refMask);
        
        for k = 1:size(refSegs, 1)
            idx = refSegs(k,1):refSegs(k,2);
            h.refLine(i) = plot(ax, x(idx), r(idx), ':', ...
                'Color', c, 'LineWidth', opt.LineWidth);
        end
    end
end

% ---------- Apply Title ----------
if ~isempty(opt.title)
    title(ax, opt.title, 'Interpreter', 'none');
end

% ---------- Apply Legends (Separate Inputs) ----------
legHandles = [];
legLabels  = {};

% 1. Process Main Legends
if ~isempty(opt.legends)
    assert(numel(opt.legends) == N, ...
        sprintf('Length of ''legends'' (%d) must match number of series (%d).', numel(opt.legends), N));
    
    legHandles = [legHandles, h.line];
    legLabels  = [legLabels, opt.legends];
end

% 2. Process Ref Legends (Only if Ref data exists)
if ~isempty(opt.refLegends)
    if ~hasRef
        warning('plotMeanWithSemShade:UnusedRefLegends', ...
            '''refLegends'' provided but ''refMeanDat'' is empty. Ref legends ignored.');
    else
        assert(numel(opt.refLegends) == N, ...
            sprintf('Length of ''refLegends'' (%d) must match number of series (%d).', numel(opt.refLegends), N));
        
        legHandles = [legHandles, h.refLine];
        legLabels  = [legLabels, opt.refLegends];
    end
end

% 3. Create Legend if we have any handles to show
if ~isempty(legHandles)
    legend(ax, legHandles, legLabels, 'Interpreter', 'none');
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

function rc = coerceRefToCell(r, expectedN)
    if isnumeric(r)
        [~, N] = size(r);
        assert(N == expectedN, 'refMeanDat columns must match meanC columns.');
        rc = arrayfun(@(j) r(:,j), 1:N, 'UniformOutput', false);
    else
        assert(iscell(r) && numel(r) == expectedN, ...
            'refMeanDat cell array must have same numel as meanC.');
        rc = r;
    end
end

function segs = maskToSegments(mask)
    d = diff([false; mask(:); false]);
    starts = find(d==1);
    ends   = find(d==-1) - 1;
    segs = [starts ends];
end