function h = plotMotifXcorrMatrix(xcorrMat, varargin)
% PLOTMOTIFXCORRMATRIX
%   Visualize a motif-by-motif cross-correlogram matrix (imagesc) with
%   options for color scaling, axis orientation, diagonal masking, and saving.
%
%   h = plotMotifXcorrMatrix(xcorrMat, 'Name', value, ...)
%
% REQUIRED INPUT
%   xcorrMat : [K x K] symmetric cross-correlogram matrix
%
% NAME–VALUE PAIRS
%   Display
%   'cLim'         : [lo hi] color axis limits (default: symmetric auto)
%   'colormap'     : colormap name or Mx3 array (default: 'parula')
%   'yDir'         : 'normal' | 'reverse' (default: 'normal')
%   'showColorbar' : true/false (default: true)
%   'nanDiagonal'  : true/false mask diagonal as NaN (default: true)
%
%   Labels
%   'xLabel'       : string (default: 'Motif')
%   'yLabel'       : string (default: 'Motif')
%
%   Figure
%   'figPos'       : [x y w h] (default: [])
%
%   Saving
%   'printFig'     : true/false (default: false)
%   'figSaveDir'   : directory for saving figures (default: '')
%   'baseName'     : base filename prefix (default: 'motif_xcorr_matrix')
%   'header'       : session header string (e.g., 'm1045_122425') (default: '')
%   'trialTypeTag' : 'both' | 'go' | 'nogo' (default: 'both')
%   'shuffleTag'   : true/false; if true append '_shuffle' to filename (default: false)
%   'reprint'      : true/false overwrite protection (default: false)
%
% Saving rule (UPDATED):
%   If header is non-empty, figure is saved to:
%     <figSaveDir>/<header>/<baseName>_<header>_<trialTypeTag>[_shuffle].pdf
%   If header is empty:
%     <figSaveDir>/<baseName>_<trialTypeTag>[_shuffle].pdf
%
% Built-in title rule:
%   base: "motif xcorr (pos-lag)"
%   prefix depends on trialTypeTag:
%     both -> "Go and Nogo"
%     go   -> "Go"
%     nogo -> "Nogo"
%   final title: "<prefix> motif xcorr (pos-lag)"
%
% OUTPUT
%   h : struct with handles + save info

% -------------------- parse --------------------
p = inputParser;

p.addRequired('xcorrMat', @(x) isnumeric(x) && ndims(x)==2 && size(x,1)==size(x,2));

% display
p.addParameter('cLim', [], @(x) isempty(x) || (isnumeric(x)&&numel(x)==2));
p.addParameter('colormap', 'parula');
p.addParameter('yDir', 'normal', @(s) any(strcmpi(s,{'normal','reverse'})));
p.addParameter('showColorbar', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('nanDiagonal', true, @(x)islogical(x)&&isscalar(x));

% labels
p.addParameter('xLabel', 'Motif', @(s)ischar(s)||isstring(s));
p.addParameter('yLabel', 'Motif', @(s)ischar(s)||isstring(s));

% figure
p.addParameter('figPos', [], @(x) isempty(x) || (isnumeric(x)&&numel(x)==4));

% saving
p.addParameter('printFig', false, @(x)islogical(x)&&isscalar(x));
p.addParameter('figSaveDir', '', @(s)ischar(s)||isstring(s));
p.addParameter('baseName', 'motif_xcorr_matrix', @(s)ischar(s)||isstring(s));
p.addParameter('header', '', @(s)ischar(s)||isstring(s));
p.addParameter('trialTypeTag', 'both', @(s) any(strcmpi(string(s), ["both","go","nogo"])));
p.addParameter('shuffleTag', false, @(x)islogical(x)&&isscalar(x));
p.addParameter('reprint', false, @(x)islogical(x)&&isscalar(x));

p.parse(xcorrMat, varargin{:});
opt = p.Results;

% -------------------- prep matrix (mask diagonal) --------------------
M = xcorrMat;
if opt.nanDiagonal
    M(1:size(M,1)+1:end) = NaN;
end

% -------------------- figure --------------------
h.fig = figure('Color','w');
if ~isempty(opt.figPos)
    set(h.fig, 'Position', opt.figPos);
end

h.ax = axes('Parent', h.fig);
imagesc(h.ax, M);
axis(h.ax, 'square');
set(h.ax, 'YDir', opt.yDir);

pbaspect(h.ax, [1 1 1]);

colormap(h.ax, opt.colormap);

% color limits (ignore NaNs)
if isempty(opt.cLim)
    mx = max(abs(M(:)), [], 'omitnan');
    if ~isfinite(mx) || mx==0, mx = 1; end
    caxis(h.ax, [-mx mx]);
else
    caxis(h.ax, opt.cLim);
end

% built-in title
trialTag = lower(string(opt.trialTypeTag));
switch trialTag
    case "both"
        prefix = "Go and Nogo";
    case "go"
        prefix = "Go";
    case "nogo"
        prefix = "Nogo";
end
title(h.ax, prefix + " motif xcorr (pos-lag)", 'Interpreter','none');

xlabel(h.ax, opt.xLabel);
ylabel(h.ax, opt.yLabel);

if opt.showColorbar
    h.cb = colorbar(h.ax);

    % ticks every 0.1
    clim = caxis(h.ax);
    t0 = ceil(clim(1)*10)/10;
    t1 = floor(clim(2)*10)/10;
    if t1 >= t0
        h.cb.Ticks = t0:0.1:t1;
    else
        h.cb.Ticks = clim;
    end
else
    h.cb = gobjects(0);
end

% -------------------- save --------------------
h.save = struct('didSave', false, 'file', '', 'figSaveDir', '', 'subDir', '', 'baseName', '');

if opt.printFig
    figSaveDir = char(string(opt.figSaveDir));
    assert(~isempty(figSaveDir), 'figSaveDir must be specified when printFig=true.');

    % ensure base dir exists
    if ~exist(figSaveDir,'dir')
        mkdir(figSaveDir);
    end

    header = strtrim(string(opt.header));

    % UPDATED: save under <figSaveDir>/<header>/ if header provided
    if strlength(header) > 0
        subDir = fullfile(figSaveDir, char(header));
        if ~exist(subDir, 'dir')
            mkdir(subDir);
        end
        outDir = subDir;
    else
        outDir = figSaveDir;
        subDir = '';
    end

    baseName = string(opt.baseName);

    parts = string.empty(1,0);
    parts(end+1) = baseName; %#ok<AGROW>
    if strlength(header) > 0
        parts(end+1) = header; %#ok<AGROW>
    end
    parts(end+1) = trialTag; %#ok<AGROW>

    fnBase = strjoin(parts, "_");
    if opt.shuffleTag
        fnBase = fnBase + "_shuffle";
    end

    outFile = fullfile(outDir, char(fnBase + ".pdf"));

    if exist(outFile,'file') && ~opt.reprint
        % skip
    else
        set(h.fig, 'InvertHardcopy','off');
        print(h.fig, outFile, '-dpdf', '-painters', '-bestfit');
        h.save.didSave    = true;
        h.save.file       = outFile;
        h.save.figSaveDir = figSaveDir;
        h.save.subDir     = outDir;
        h.save.baseName   = char(fnBase);
    end
end

end
