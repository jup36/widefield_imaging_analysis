function fig = plotFastSlowAbsSlopeHeatmap(dataFast, dataSlow, timeX, idsFast, idsSlow, varargin)
%PLOTFASTSLOWABSSLOPEHEATMAP
% Two-panel heatmap (fast learners on top, slow learners on bottom) using
% TWO DIFFERENT COLORMAPS -- one hue ramp per group -- so that group
% membership is visually encoded by color/hue, not just row position.
%
% This exploits per-axes colormaps (supported natively since MATLAB
% R2014b): each stacked axes gets its own colormap() call, so a single
% figure can show two different hue ramps while sharing the same x-axis
% (time) and the same color LIMITS (clims), keeping magnitude directly
% comparable across groups -- only the hue differs.
%
% USAGE
%   plotFastSlowAbsSlopeHeatmap(mean_muNoGo_absSlope_onNoGoAxes_fasts_ds, ...
%       mean_muNoGo_absSlope_onNoGoAxes_slows_ds, timeX_ds, fastC, slowC, ...
%       'clims', [0 1], ...
%       'colorFast', [0.10 0.45 1.00], ...   % vivid blue
%       'colorSlow', [1.00 0.55 0.10], ...   % vivid orange
%       'titleStr', '|muNoGo learning slope| on NoGo-related axes (correct trials)', ...
%       'figSaveLogic', true, ...
%       'figSaveDir', figSaveDirSlope, ...
%       'figName', "muNoGo_absSlope_onNoGoAxes_correctTrials_fastSlowColor");
%
% REQUIRED
%   dataFast : [nFast x nTime] matrix, one row per fast-learner mouse
%   dataSlow : [nSlow x nTime] matrix, one row per slow-learner mouse
%   timeX    : [1 x nTime] time vector (shared by both groups)
%   idsFast  : cellstr/string, length nFast, mouse ID labels for dataFast rows
%   idsSlow  : cellstr/string, length nSlow, mouse ID labels for dataSlow rows
%
% NAME-VALUE
%   'clims'         : [0 1] (default) -- shared color limits for BOTH panels
%   'colorFast'     : [0.10 0.45 1.00] (default, vivid blue) -- RGB triplet,
%                     the "full intensity" end of the fast-group ramp
%   'colorSlow'     : [1.00 0.55 0.10] (default, vivid orange) -- RGB triplet,
%                     the "full intensity" end of the slow-group ramp
%   'nColors'       : 256 (default) -- colormap resolution
%   'titleStr'      : '' (default)
%   'xlabelStr'     : 'Time (s)' (default)
%   'panelHeightRatio' : [] (default) -- if empty, panel heights are
%                     proportional to nFast/nSlow row counts; otherwise
%                     pass e.g. [0.5 0.5] to force equal height
%   'showColorbar'  : true (default) -- shows two colorbars, one per group,
%                     labeled "Fast" / "Slow"
%   'figSaveLogic'  : false (default)
%   'figSaveDir'    : '' (default; required if figSaveLogic=true)
%   'figName'       : "fastSlowHeatmap" (default)
%
% OUTPUT
%   fig : figure handle

%% -------------------- parse --------------------
ip = inputParser;
ip.FunctionName = mfilename;

addRequired(ip, 'dataFast', @(x) isnumeric(x) && ismatrix(x));
addRequired(ip, 'dataSlow', @(x) isnumeric(x) && ismatrix(x));
addRequired(ip, 'timeX', @(x) isnumeric(x) && isvector(x));
addRequired(ip, 'idsFast', @(x) iscell(x) || isstring(x));
addRequired(ip, 'idsSlow', @(x) iscell(x) || isstring(x));

addParameter(ip, 'clims', [0 1], @(x) isnumeric(x) && numel(x)==2);
addParameter(ip, 'colorFast', [0.10 0.45 1.00], @(x) isnumeric(x) && numel(x)==3);
addParameter(ip, 'colorSlow', [1.00 0.55 0.10], @(x) isnumeric(x) && numel(x)==3);
addParameter(ip, 'nColors', 256, @(x) isnumeric(x) && isscalar(x) && x>=2);
addParameter(ip, 'titleStr', '', @(x) ischar(x) || isstring(x));
addParameter(ip, 'xlabelStr', 'Time (s)', @(x) ischar(x) || isstring(x));
addParameter(ip, 'panelHeightRatio', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
addParameter(ip, 'showColorbar', true, @(x) islogical(x) && isscalar(x));
addParameter(ip, 'figSaveLogic', false, @(x) islogical(x) && isscalar(x));
addParameter(ip, 'figSaveDir', '', @(x) ischar(x) || isstring(x));
addParameter(ip, 'figName', "fastSlowHeatmap", @(x) ischar(x) || isstring(x));

parse(ip, dataFast, dataSlow, timeX, idsFast, idsSlow, varargin{:});
P = ip.Results;

idsFast = cellstr(string(idsFast(:)));
idsSlow = cellstr(string(idsSlow(:)));

nFast = size(dataFast,1);
nSlow = size(dataSlow,1);
assert(numel(idsFast)==nFast, 'idsFast length must match rows of dataFast.');
assert(numel(idsSlow)==nSlow, 'idsSlow length must match rows of dataSlow.');

if isempty(P.panelHeightRatio)
    total = nFast + nSlow;
    heightRatio = [nFast/total, nSlow/total];
else
    heightRatio = P.panelHeightRatio(:)' / sum(P.panelHeightRatio);
end

%% -------------------- build colormaps: black -> group color --------------------
cmapFast = blackToColorRamp_(P.colorFast, P.nColors);
cmapSlow = blackToColorRamp_(P.colorSlow, P.nColors);

%% -------------------- layout: two stacked axes, no gap --------------------
fig = figure('Color', 'w');

leftMargin   = 0.12;
rightMargin  = 0.18;   % room for two colorbars
bottomMargin = 0.12;
topMargin    = 0.10;
panelWidth   = 1 - leftMargin - rightMargin;
panelHeightTotal = 1 - bottomMargin - topMargin;

hSlow = panelHeightTotal * heightRatio(2);
hFast = panelHeightTotal * heightRatio(1);

posSlow = [leftMargin, bottomMargin,          panelWidth, hSlow];
posFast = [leftMargin, bottomMargin + hSlow,  panelWidth, hFast];

% ---- Fast panel (top) ----
axFast = axes('Parent', fig, 'Position', posFast);
imagesc(axFast, timeX, 1:nFast, dataFast, P.clims);
colormap(axFast, cmapFast);
set(axFast, 'YTick', 1:nFast, 'YTickLabel', idsFast, 'XTickLabel', []);
box(axFast, 'off');
if ~isempty(P.titleStr)
    title(axFast, P.titleStr, 'Interpreter', 'none');
end

% ---- Slow panel (bottom) ----
axSlow = axes('Parent', fig, 'Position', posSlow);
imagesc(axSlow, timeX, 1:nSlow, dataSlow, P.clims);
colormap(axSlow, cmapSlow);
set(axSlow, 'YTick', 1:nSlow, 'YTickLabel', idsSlow);
xlabel(axSlow, P.xlabelStr);
box(axSlow, 'off');

linkaxes([axFast, axSlow], 'x');
xlim(axFast, [min(timeX) max(timeX)]);

%% -------------------- colorbars (one per group) --------------------
if P.showColorbar
    cbFast = colorbar(axFast);
    cbFast.Position = [leftMargin + panelWidth + 0.03, posFast(2), 0.03, posFast(4)];
    cbFast.Label.String = 'Fast |slope|';

    cbSlow = colorbar(axSlow);
    cbSlow.Position = [leftMargin + panelWidth + 0.03, posSlow(2), 0.03, posSlow(4)];
    cbSlow.Label.String = 'Slow |slope|';
end

%% -------------------- save --------------------
if P.figSaveLogic
    if isempty(P.figSaveDir)
        error('figSaveDir must be provided when figSaveLogic=true.');
    end
    if ~exist(P.figSaveDir, 'dir')
        mkdir(P.figSaveDir);
    end
    fpath = fullfile(P.figSaveDir, char(string(P.figName)) + ".png");
    exportgraphics(fig, fpath, 'Resolution', 300);

    fpathFig = fullfile(P.figSaveDir, char(string(P.figName)) + ".fig");
    savefig(fig, fpathFig);
end

end

%% ======================================================================
function cmap = blackToColorRamp_(colorRGB, n)
% Linear ramp from black [0 0 0] to the given RGB triplet, n x 3.
colorRGB = colorRGB(:)';
cmap = [linspace(0, colorRGB(1), n)', ...
        linspace(0, colorRGB(2), n)', ...
        linspace(0, colorRGB(3), n)'];
end