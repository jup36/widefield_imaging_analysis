function h = plot_beta_with_labels(beta, X_names, motifId, varargin)
% PLOT_BETA_WITH_LABELS  Bar plot of β for one motif with compact x labels.
%
% Adds:
%   - Color/style by tone predictor type:
%       tOnG   : bluish filled
%       tOffG  : bluish edge-only (no fill)
%       tOnNg  : redish filled
%       tOffNg : redish edge-only (no fill)
%   - If a predictor type has >4 members (e.g., tOnNg1..tOnNg9),
%     x-axis labels show only odd indices (2,4,6,... are blank labels).

% -------- parse args --------
p = inputParser;
p.addParameter('sort', 'none', @(s) any(strcmpi(s, {'none','absdesc'})));
p.addParameter('title', '', @(s) ischar(s) || isstring(s));
p.addParameter('predictorI', [], @(x) isempty(x) || islogical(x));
p.addParameter('figureScaleFactor', 2, @(x) isnumeric(x) && x>0);
p.addParameter('figureWidthFactor', 3, @(x) isnumeric(x) && x>0);
p.addParameter('visible', 'on', @(s) any(strcmpi(s, {'on', 'off'})));

% saving options (as you requested previously)
p.addParameter('figSaveDir', {}, @(x) isempty(x) || ischar(x) || isstring(x) || iscell(x));
p.addParameter('header', '', @(s) ischar(s) || isstring(s));
p.addParameter('figSaveKeyword', '', @(s) ischar(s) || isstring(s));

p.parse(varargin{:});
opt = p.Results;

% normalize figSaveDir
figSaveDir = opt.figSaveDir;
if iscell(figSaveDir)
    if isempty(figSaveDir), figSaveDir = ""; else, figSaveDir = string(figSaveDir{1}); end
else
    figSaveDir = string(figSaveDir);
end
header = string(opt.header);
figSaveKeyword = string(opt.figSaveKeyword);

% -------- sanity checks --------
[P, K] = size(beta);
assert(iscellstr(X_names) && numel(X_names)==P, 'X_names must be 1xP cellstr.');
assert(isscalar(motifId) && motifId>=1 && motifId<=K, 'motifId out of range.');

% Validate predictorI when provided
if ~isempty(opt.predictorI)
    assert(islogical(opt.predictorI) && numel(opt.predictorI)==P, 'predictorI must be logical(1xP).');
    assert(any(opt.predictorI), 'predictorI selects zero predictors.');
    sel = opt.predictorI(:).';  % row logical mask
else
    sel = true(1,P);
end

bAll  = beta(:, motifId);
b     = bAll(sel);

% -------- build compact labels --------
abbrAll = cellfun(@abbrev_name, X_names, 'uni', 0);
abbr    = abbrAll(sel);

% optional sorting by |β| (within the selected subset)
nSel = numel(b);
ord  = 1:nSel;
if strcmpi(opt.sort, 'absdesc')
    [~, ord] = sort(abs(b), 'descend');
end

bPlot    = b(ord);
abbrPlot = abbr(ord);

% -------- classify predictor type + index from abbreviated label --------
% expected: tOnNg7, tOffG3, etc.
[typeKey, typeIdxNum] = parse_abbr_type_(abbrPlot);

% how many of each type in the plotted set?
typeList = unique(typeKey(~strcmp(typeKey,"other")));
typeCount = struct();
for iT = 1:numel(typeList)
    t = typeList(iT);
    typeCount.(char(t)) = sum(typeKey == t);
end

% build x tick labels with odd-only rule when crowded
xtlbl = abbrPlot;
for j = 1:numel(xtlbl)
    t = typeKey(j);
    idxNum = typeIdxNum(j);
    if t ~= "other" && isfield(typeCount, char(t)) && typeCount.(char(t)) > 4
        if ~isnan(idxNum) && mod(idxNum,2)==0
            xtlbl{j} = ''; % blank even indices for that type
        end
    end
end

% -------- plot --------
h.fig = figure('Color','w', 'Visible',opt.visible);
set(h.fig, 'Units', 'normalized');
pos = get(h.fig, 'Position');
pos(3:4) = pos(3:4) * opt.figureScaleFactor;
set(h.fig, 'Position', pos);

h.ax = axes('Parent',h.fig); hold(h.ax,'on');

% draw each bar individually so we can mix filled vs edge-only
h.bars = gobjects(nSel,1);
x = 1:nSel;

% ---- color definitions (define ONCE) ----
blueish = [0.30 0.50 0.90];
redish  = [0.90 0.30 0.30];
grayish = [0.55 0.55 0.55];   % neutral monotone gray

% ---- styling knobs ----
filledEdgeColor   = 'none';   % keep filled bars clean
filledLineWidth   = 0.5;      % ignored if EdgeColor='none'
edgeOnlyLineWidth = 2;      % <<< make empty bars bold (adjust as you like)
barWidth          = 0.85;

for j = 1:nSel
    bj = bar(h.ax, x(j), bPlot(j), barWidth);
    h.bars(j) = bj;

    % default: non-tone predictors (gray filled)
    faceC = grayish;
    edgeC = filledEdgeColor;
    lw    = filledLineWidth;

    switch typeKey(j)

        % ---------- Go ----------
        case "tOnG"
            faceC = blueish;
            edgeC = filledEdgeColor;
            lw    = filledLineWidth;

        case "tOffG"
            faceC = 'none';
            edgeC = blueish;
            lw    = edgeOnlyLineWidth;

        % ---------- NoGo ----------
        case "tOnNg"
            faceC = redish;
            edgeC = filledEdgeColor;
            lw    = filledLineWidth;

        case "tOffNg"
            faceC = 'none';
            edgeC = redish;
            lw    = edgeOnlyLineWidth;

        % ---------- Non-tone predictors ----------
        otherwise
            faceC = grayish;
            edgeC = filledEdgeColor;   % set to grayish if you want outlines too
            lw    = filledLineWidth;
    end

    set(bj, 'FaceColor', faceC, 'EdgeColor', edgeC, 'LineWidth', lw);
end


% zero line
plot(h.ax, [0.5, nSel+0.5], [0 0], 'k-', 'LineWidth', 0.8);

% ticks & labels
xticks(h.ax, 1:nSel);
xticklabels(h.ax, xtlbl);
xtickangle(h.ax, 45);
set(h.ax, 'TickLabelInterpreter','none');

xlabel(h.ax, 'Predictors (abbreviated)');
ylabel(h.ax, sprintf('\\beta (motif %d)', motifId));

if strlength(string(opt.title))>0
    title(h.ax, opt.title, 'Interpreter','none');
else
    title(h.ax, sprintf('Motif %d: predictor weights', motifId), 'Interpreter','none');
end

box(h.ax,'off');
grid(h.ax,'on');
set(gca, 'TickDir', 'out')

% Aspect & tight layout
r = pbaspect;
pbaspect([opt.figureWidthFactor*r(1) r(2) r(3)]);
axis tight

% -------- save figure (optional) --------
if strlength(figSaveDir) > 0
    if ~isfolder(figSaveDir)
        mkdir(figSaveDir);
    end
    dateStr = char(datetime("today","Format","MMddyy"));

    parts = strings(0,1);
    if strlength(header) > 0,          parts(end+1,1) = header; end
    parts(end+1,1) = "glmBeta";
    if strlength(figSaveKeyword) > 0,  parts(end+1,1) = figSaveKeyword; end
    parts(end+1,1) = "motif" + string(motifId);
    parts(end+1,1) = dateStr;

    figSaveName = strjoin(parts, "_");
    print(h.fig, fullfile(figSaveDir, figSaveName), '-dpdf', '-painters', '-bestfit');
end

end % function


% ===== helper: make names compact and consistent =====
function s = abbrev_name(name)
idxTok = regexp(name, '_rc(\d+)$', 'tokens', 'once');
if ~isempty(idxTok)
    idx = str2double(idxTok{1});
    base = regexprep(name, '_rc\d+$', '');
    idxStr = num2str(idx);
else
    base = name;
    idxStr = '';
end

base = strrep(base, 'toneOnNoGo',  'tOnNg');
base = strrep(base, 'toneOffNoGo', 'tOffNg');
base = strrep(base, 'toneOnGo',    'tOnG');
base = strrep(base, 'toneOffGo',   'tOffG');

base = strrep(base, 'airpuff',     'pun');
base = strrep(base, 'water',       'rwd');
base = strrep(base, 'combinedLicksVid', 'lick');
base = strrep(base, 'postToneLicksVid', 'lickP');
base = strrep(base, 'periToneLicksVid', 'lickB');
base = strrep(base, 'locomVel',    'loc');
base = strrep(base, 'nosetip',     'nose');
base = strrep(base, 'whisker',     'whisk');
base = strrep(base, 'pupil',       'pupil');

base = regexprep(base, '_+', '');

if ~isempty(idxStr)
    s = sprintf('%s%s', base, idxStr);
else
    s = base;
end
end


function [typeKey, idxNum] = parse_abbr_type_(abbrC)
% Parse abbreviated labels like tOnNg7, tOffG3, etc.
% Returns:
%   typeKey : string array (tOnG/tOffG/tOnNg/tOffNg/other)
%   idxNum  : numeric index (NaN if none)
n = numel(abbrC);
typeKey = repmat("other", n, 1);
idxNum  = nan(n,1);

for i = 1:n
    s = string(abbrC{i});
    tok = regexp(s, '^(tOnG|tOffG|tOnNg|tOffNg)(\d+)$', 'tokens', 'once');
    if ~isempty(tok)
        typeKey(i) = string(tok{1});
        idxNum(i)  = str2double(tok{2});
    end
end
end
