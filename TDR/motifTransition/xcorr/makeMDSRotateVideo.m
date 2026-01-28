function outFile = makeMDSRotateVideo(h, outDir, baseName, startView, varargin)
% makeMDSRotateVideo
% Create an MP4 (or GIF) that rotates an EXISTING open MDS figure.
%
% PHASES:
%   1) Horizontal rotation: azimuth sweep, elevation fixed
%   2) Vertical tilt: elevation sweep, azimuth fixed (at end of horizontal)
%
% REQUIRED:
%   h         : handle or struct returned by your plotting function (optional; can be [])
%   outDir    : output directory
%   baseName  : filename without extension
%   startView : [az el] starting view (required)
%
% RECOMMENDED:
%   Pass 'ax', gca (after clicking the correct axes)
%
% EXAMPLE:
% ax = gca; % click inside 3D axes first
% outFile = makeMDSRotateVideo(h_full, figSaveDir, ...
%   'xcorr_mds_rotate', [-111 28], ...
%   'ax', ax, ...
%   'fps', 30, ...
%   'horizSeconds', 6, 'horizDegrees', 360, ...
%   'vertSeconds', 4, 'vertRange', [-5 45], ...
%   'renderer','opengl', 'figPos',[100 100 1200 900]);

% -------------------- Parse inputs --------------------
p = inputParser;
p.addRequired('h'); % may be [] or struct or graphics handle
p.addRequired('outDir', @(x)ischar(x)||isstring(x));
p.addRequired('baseName', @(x)ischar(x)||isstring(x));
p.addRequired('startView', @(x)isnumeric(x)&&numel(x)==2);

p.addParameter('ax', [], @(x) isempty(x) || isgraphics(x,'axes'));
p.addParameter('format', 'mp4', @(s)ischar(s)||isstring(s));     % 'mp4' or 'gif'
p.addParameter('fps', 30, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('renderer', 'opengl', @(s)ischar(s)||isstring(s));
p.addParameter('figPos', [], @(x)isnumeric(x)&&(isempty(x)||numel(x)==4));
p.addParameter('lockVis3d', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('ease', 'cosine', @(s)ischar(s)||isstring(s));    % 'cosine' or 'linear'

% Horizontal rotation
p.addParameter('horizDegrees', 360, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('horizSeconds', 6, @(x)isnumeric(x)&&isscalar(x)&&x>0);

% Vertical rotation
p.addParameter('vertRange', [-5 45], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('vertSeconds', 4, @(x)isnumeric(x)&&isscalar(x)&&x>0);

p.parse(h, outDir, baseName, startView, varargin{:});
opt = p.Results;

format = lower(string(opt.format));
ease   = lower(string(opt.ease));
fps    = opt.fps;

% -------------------- Resolve fig & ax (figure must be open) --------------------
[figGuess, axGuess] = localResolveFigAx(opt.h);

ax = opt.ax;
if isempty(ax), ax = axGuess; end
if isempty(ax) || ~isgraphics(ax,'axes')
    ax = gca; % last resort
end
if isempty(ax) || ~isgraphics(ax,'axes')
    error('makeMDSRotateVideo:NoAxes', 'Could not resolve a valid axes handle. Pass ''ax'', gca explicitly.');
end

fig = ancestor(ax,'figure');
if isempty(fig) || ~isgraphics(fig,'figure')
    if ~isempty(figGuess) && isgraphics(figGuess,'figure')
        fig = figGuess;
    else
        error('makeMDSRotateVideo:NoFigure', 'No valid open figure found. Keep the figure open.');
    end
end

% -------------------- Configure figure BEFORE writing frames --------------------
set(fig, 'Renderer', char(opt.renderer));
set(fig, 'Units','pixels');
if ~isempty(opt.figPos)
    set(fig, 'Position', opt.figPos);
end
drawnow;

if opt.lockVis3d
    try, axis(ax,'vis3d'); catch, end
end

startView = opt.startView(:).';
view(ax, startView);
drawnow;

% -------------------- Output path --------------------
outDir = char(opt.outDir);
if ~exist(outDir,'dir'), mkdir(outDir); end
baseName = char(opt.baseName);

switch format
    case "mp4"
        outFile = fullfile(outDir, baseName + ".mp4");
    case "gif"
        outFile = fullfile(outDir, baseName + ".gif");
    otherwise
        error('format must be ''mp4'' or ''gif''.');
end

% -------------------- Trajectories --------------------
% Horizontal (az sweep)
horizFrames = max(2, round(opt.horizSeconds * fps));
tH  = linspace(0,1,horizFrames);
teH = localEase(tH, ease);

az0 = startView(1);
el0 = startView(2);

azH = az0 + opt.horizDegrees * teH;
elH = el0 * ones(size(azH));

% Vertical (el sweep; az fixed at end of horizontal)
vertFrames = max(2, round(opt.vertSeconds * fps));
tV  = linspace(0,1,vertFrames);
teV = localEase(tV, ease);

azV = azH(end) * ones(size(teV));
elV = opt.vertRange(1) + (opt.vertRange(2)-opt.vertRange(1))*teV;

% -------------------- Write --------------------
switch format
    case "mp4"
        vw = VideoWriter(char(outFile), 'MPEG-4');
        vw.FrameRate = fps;
        open(vw);

        % --- Grab first frame, define target size (even dims) ---
        fr0 = getframe(ax);
        [H0, W0, ~] = size(fr0.cdata);
        Ht = H0 - mod(H0,2);   % make even
        Wt = W0 - mod(W0,2);   % make even

        % write first frame (cropped to even target)
        writeVideo(vw, localFrameToEven(fr0, Ht, Wt));

        % Phase 1: horizontal
        for i = 1:numel(azH)
            if ~isgraphics(fig,'figure') || ~isgraphics(ax,'axes')
                error('Figure/axes was closed during writing. Keep it open.');
            end
            view(ax, [azH(i) elH(i)]);
            drawnow;
            fr = getframe(ax);
            writeVideo(vw, localFrameToEven(fr, Ht, Wt));
        end

        % Phase 2: vertical
        for i = 1:numel(azV)
            if ~isgraphics(fig,'figure') || ~isgraphics(ax,'axes')
                error('Figure/axes was closed during writing. Keep it open.');
            end
            view(ax, [azV(i) elV(i)]);
            drawnow;
            fr = getframe(ax);
            writeVideo(vw, localFrameToEven(fr, Ht, Wt));
        end

        close(vw);


    case "gif"
        delay = 1 / fps;
        frameCount = 0;

        % warmup
        localWriteGifFrame(ax, outFile, delay, frameCount);
        frameCount = frameCount + 1;

        for i = 1:numel(azH)
            view(ax, [azH(i) elH(i)]);
            drawnow;
            localWriteGifFrame(ax, outFile, delay, frameCount);
            frameCount = frameCount + 1;
        end

        for i = 1:numel(azV)
            view(ax, [azV(i) elV(i)]);
            drawnow;
            localWriteGifFrame(ax, outFile, delay, frameCount);
            frameCount = frameCount + 1;
        end
end

fprintf('Saved: %s\n', outFile);

end % ===== main =====


% ======================= helpers =======================

function [fig, ax] = localResolveFigAx(h)
fig = [];
ax  = [];

if isstruct(h)
    try
        if isfield(h,'base')
            if isfield(h.base,'fig') && isgraphics(h.base.fig,'figure'), fig = h.base.fig; end
            if isfield(h.base,'ax')  && isgraphics(h.base.ax,'axes'),    ax  = h.base.ax;  end
        end
        if isempty(fig) && isfield(h,'fig') && isgraphics(h.fig,'figure'), fig = h.fig; end
        if isempty(ax)  && isfield(h,'ax')  && isgraphics(h.ax,'axes'),    ax  = h.ax;  end
    catch
    end
elseif isgraphics(h)
    try
        fig = ancestor(h,'figure');
        ax  = ancestor(h,'axes');
    catch
    end
end

% fallbacks
if isempty(fig) || ~isgraphics(fig,'figure')
    try, fig = gcf; catch, fig = []; end
end
if isempty(ax) || ~isgraphics(ax,'axes')
    try, ax = gca; catch, ax = []; end
end
end

function te = localEase(t, ease)
switch lower(string(ease))
    case "cosine"
        te = 0.5 - 0.5*cos(pi*t);
    case "linear"
        te = t;
    otherwise   
        te = t;
end
end

function localWriteGifFrame(ax, outFile, delay, frameCount)
fr = getframe(ax);
[im, map] = rgb2ind(fr.cdata, 256);
if frameCount == 0
    imwrite(im, map, char(outFile), 'gif', 'LoopCount', inf, 'DelayTime', delay);
else
    imwrite(im, map, char(outFile), 'gif', 'WriteMode', 'append', 'DelayTime', delay);
end
end

function frOut = localFrameToEven(frIn, Ht, Wt)
% Ensure frIn.cdata is exactly Ht-by-Wt (and both even), by center-cropping or padding.
im = frIn.cdata;
[H, W, C] = size(im);

% If larger, center-crop
if H > Ht
    y0 = floor((H - Ht)/2) + 1;
    im = im(y0:y0+Ht-1, :, :);
elseif H < Ht
    padTop = floor((Ht - H)/2);
    padBot = (Ht - H) - padTop;
    im = padarray(im, [padTop 0], 0, 'pre');
    im = padarray(im, [padBot 0], 0, 'post');
end

if W > Wt
    x0 = floor((W - Wt)/2) + 1;
    im = im(:, x0:x0+Wt-1, :);
elseif W < Wt
    padLeft  = floor((Wt - W)/2);
    padRight = (Wt - W) - padLeft;
    im = padarray(im, [0 padLeft], 0, 'pre');
    im = padarray(im, [0 padRight], 0, 'post');
end

% sanity
im = im(1:Ht, 1:Wt, 1:C);

frOut = frIn;
frOut.cdata = im;
end
