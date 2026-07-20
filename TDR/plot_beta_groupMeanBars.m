function h = plot_beta_groupMeanBars(beta, X_names, motifId, varargin)
%PLOT_BETA_GROUPMEANBARS  Plot mean β grouped by predictor type (collapsed bars).
%
% Groups:
%   - tOnG   : mean over all toneOnGo_rc##
%   - tOffG  : mean over all toneOffGo_rc##
%   - tOnNg  : mean over all toneOnNoGo_rc##
%   - tOffNg : mean over all toneOffNoGo_rc##
%   - other  : (default) mean over all remaining predictors
%
% Adds (per your request):
%   - No "n=x" annotations on bars
%   - xticklabels rotated 45 degrees
%   - 'plotTonePredictors' (default=true). If false:
%       exclude tone predictors entirely and plot only "other" bars.
%
% ADDED (for reuse with synthetic/group-mean data, e.g.
% plotMotifOtherBetaByGroup.m -- see 'yLim' and 'saveMotifId' below):
%   - 'yLim'         : [min max] override for the y-axis, applied BEFORE
%                      the figure is saved (so harmonizing the y-axis
%                      across multiple calls/panels actually shows up in
%                      the saved PDF, not just the on-screen figure).
%                      Default: [] (auto-scaled via axis tight, as before).
%   - 'saveMotifId'  : override for the motif number used in the SAVED
%                      FILENAME only (data indexing/plotting still use
%                      the real motifId argument). Use this when beta is
%                      synthetic/single-column data standing in for some
%                      other real motif (e.g. a group-mean profile), so
%                      the file is tagged with the motif it actually
%                      represents instead of its column index.
%                      Default: [] (use motifId, as before).

% ---------------- parse ----------------
p = inputParser;
p.addParameter('title','', @(s)ischar(s)||isstring(s));
p.addParameter('predictorI', [], @(x) isempty(x) || (islogical(x) && isvector(x)));
p.addParameter('SplitOtherByBase', false, @(x)islogical(x)&&isscalar(x));
p.addParameter('plotTonePredictors', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('visible','on', @(s) any(strcmpi(s, {'on','off'})));

p.addParameter('figSaveDir', {}, @(x) isempty(x) || ischar(x) || isstring(x) || iscell(x));
p.addParameter('header', '', @(s) ischar(s) || isstring(s));
p.addParameter('figSaveKeyword', '', @(s) ischar(s) || isstring(s));

p.addParameter('yLim', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
p.addParameter('saveMotifId', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));

p.parse(varargin{:});
opt = p.Results;

% Single source of truth for "which motif number to show/save this as" --
% used by BOTH the y-axis label AND the saved filename. Falls back to
% motifId (the real column-index argument) unless overridden.
motifForDisplay = motifId;
if ~isempty(opt.saveMotifId)
    motifForDisplay = opt.saveMotifId;
end

% normalize figSaveDir
figSaveDir = opt.figSaveDir;
if iscell(figSaveDir)
    if isempty(figSaveDir), figSaveDir = ""; else, figSaveDir = string(figSaveDir{1}); end
else
    figSaveDir = string(figSaveDir);
end
header         = string(opt.header);
figSaveKeyword = string(opt.figSaveKeyword);

% ---------------- sanity ----------------
[P,K] = size(beta);
assert(iscellstr(X_names) && numel(X_names)==P, 'X_names must be 1xP cellstr.');
assert(isscalar(motifId) && motifId>=1 && motifId<=K, 'motifId out of range.');

if ~isempty(opt.predictorI)
    sel = opt.predictorI(:).';
    assert(numel(sel)==P, 'predictorI must be logical(1xP).');
    assert(any(sel), 'predictorI selects zero predictors.');
else
    sel = true(1,P);
end

bAll     = beta(:, motifId);
bSel     = bAll(sel);
namesSel = X_names(sel);

% ---------------- classify each predictor ----------------
[typeKey, baseKey] = classify_predictor_names_(namesSel);

% optionally exclude tone predictors entirely
if ~opt.plotTonePredictors
    keep = (typeKey == "other");
    bSel     = bSel(keep);
    typeKey  = typeKey(keep);
    baseKey  = baseKey(keep);
    namesSel = namesSel(keep); %#ok<NASGU>
end

% ---------------- define which groups to plot ----------------
toneGroups = ["tOnG","tOffG","tOnNg","tOffNg"];

groupNames = strings(0,1);
groupMeans = [];
groupNs    = []; %#ok<NASGU>  % kept for output bookkeeping, not plotted

% (1) tone groups (only if enabled)
if opt.plotTonePredictors
    for g = 1:numel(toneGroups)
        tg = toneGroups(g);
        m = (typeKey == tg);
        if any(m)
            groupNames(end+1,1) = tg; %#ok<AGROW>
            groupMeans(end+1,1) = mean(bSel(m), 'omitnan'); %#ok<AGROW>
            groupNs(end+1,1)    = sum(m); %#ok<AGROW>
        end
    end
end

% (2) other group(s)
mOther = (typeKey == "other");
if any(mOther)
    if ~opt.SplitOtherByBase
        groupNames(end+1,1) = "other";
        groupMeans(end+1,1) = mean(bSel(mOther), 'omitnan');
        groupNs(end+1,1)    = sum(mOther);
    else
        u = unique(baseKey(mOther), 'stable');
        for i = 1:numel(u)
            mm = mOther & (baseKey == u(i));
            groupNames(end+1,1) = u(i); %#ok<AGROW>
            groupMeans(end+1,1) = mean(bSel(mm), 'omitnan'); %#ok<AGROW>
            groupNs(end+1,1)    = sum(mm); %#ok<AGROW>
        end
    end
end

assert(~isempty(groupNames), 'No predictors survived selection to plot (check predictorI / plotTonePredictors).');

% ---------------- plot ----------------
h = struct();
h.fig = figure('Color','w', 'Visible', opt.visible);
h.ax  = axes('Parent',h.fig); hold(h.ax,'on');

grayFill = [0.55 0.55 0.55];

barWidth = 0.35; 
h.bars = bar(h.ax, groupMeans, barWidth, 'FaceColor', grayFill, 'EdgeColor', 'none');

% zero line
plot(h.ax, [0.5, numel(groupMeans)+0.5], [0 0], 'k-', 'LineWidth', 0.8);

xticks(h.ax, 1:numel(groupNames));
xticklabels(h.ax, cellstr(groupNames));
set(h.ax, 'TickLabelInterpreter','none');
xtickangle(h.ax, 45);   % <<< rotated as requested

xlabel(h.ax, 'Predictor group (collapsed)');
ylabel(h.ax, sprintf('\\beta mean (motif %d)', motifForDisplay));

if strlength(string(opt.title))>0
    title(h.ax, opt.title, 'Interpreter','none');
else
    title(h.ax, sprintf('Motif %d: mean \\beta by predictor type', motifForDisplay), 'Interpreter','none');
end

box(h.ax,'off'); grid(h.ax,'on'); set(h.ax,'TickDir','out'); axis tight;

% yLim override, applied BEFORE saving so a harmonized range across
% multiple calls actually ends up in the saved PDF
if ~isempty(opt.yLim)
    ylim(h.ax, opt.yLim);
end

% ---------------- save (optional) ----------------
if strlength(figSaveDir) > 0
    if ~isfolder(figSaveDir), mkdir(figSaveDir); end
    dateStr = char(datetime("today","Format","MMddyy"));

    parts = strings(0,1);
    if strlength(header)>0, parts(end+1,1) = header; end
    parts(end+1,1) = "glmBetaGroupMean";
    if ~opt.plotTonePredictors
        parts(end+1,1) = "noTone";
    end
    if strlength(figSaveKeyword)>0, parts(end+1,1) = figSaveKeyword; end

    parts(end+1,1) = "motif" + string(motifForDisplay);
    parts(end+1,1) = dateStr;

    figSaveName = strjoin(parts, "_");
    print(h.fig, fullfile(figSaveDir, figSaveName), '-dpdf', '-painters', '-bestfit');
end

% ---------------- outputs ----------------
h.groupNames = groupNames;
h.groupMeans = groupMeans;
h.groupNs    = groupNs;
h.opt        = opt;

end

% ======================================================================
function [typeKey, baseKey] = classify_predictor_names_(nameC)
% Classify predictors into tone types or "other".
n = numel(nameC);
typeKey = repmat("other", n, 1);
baseKey = repmat("other", n, 1);

for i = 1:n
    nm = string(nameC{i});

    if contains(nm, "toneOnGo", "IgnoreCase", true)
        typeKey(i) = "tOnG";   baseKey(i) = "tOnG";
    elseif contains(nm, "toneOffGo", "IgnoreCase", true)
        typeKey(i) = "tOffG";  baseKey(i) = "tOffG";
    elseif contains(nm, "toneOnNoGo", "IgnoreCase", true)
        typeKey(i) = "tOnNg";  baseKey(i) = "tOnNg";
    elseif contains(nm, "toneOffNoGo", "IgnoreCase", true)
        typeKey(i) = "tOffNg"; baseKey(i) = "tOffNg";
    else
        % base extraction: strip trailing _rc##
        base = regexprep(nm, '_rc\d+$', '');
        base = regexprep(base, '_+', '');
        baseKey(i) = lower(base);
    end
end
end