function outFile = makeMDSRevealThenRotateVideo(h, outDir, baseName, defaultView, varargin)
% makeMDSRevealThenRotateVideo
%
% PHASES:
%   (A) BLANK start (all revealables hidden)
%   (B) REVEAL session-by-session at defaultView:
%         at step k, each trajectory shows its first k sessions (aligned column)
%   (C) One horizontal rotation
%   (D) One vertical rotation
%
% REQUIRED:
%   defaultView : [az el]
%
% STRONGLY RECOMMENDED:
%   'ax', gca   % pass axes explicitly
%
% KEY OPTIONS:
%   'revealSeconds'  (default 6)
%   'horizSeconds'   (default 6), 'horizDegrees' (default 360)
%   'vertSeconds'    (default 4), 'vertRange' (default [-5 45])
%   'fps'            (default 30)
%   'figPos'         (default []), must be set to lock frame size
%   'includeLegend'  (default true) include legend in capture
%   'hideSurfacesDuringReveal' (default false) if true, hide ellipsoid until reveal begins
%
% OUTPUT:
%   outFile (mp4)

% -------------------- Parse inputs --------------------
p = inputParser;
p.addRequired('h');
p.addRequired('outDir', @(x)ischar(x)||isstring(x));
p.addRequired('baseName', @(x)ischar(x)||isstring(x));
p.addRequired('defaultView', @(x)isnumeric(x)&&numel(x)==2);

p.addParameter('ax', [], @(x) isempty(x) || isgraphics(x));
p.addParameter('fps', 30, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('renderer', 'opengl', @(s)ischar(s)||isstring(s));
p.addParameter('figPos', [], @(x)isnumeric(x)&&(isempty(x)||numel(x)==4));
p.addParameter('lockVis3d', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('ease', 'cosine', @(s)ischar(s)||isstring(s));

% Reveal
p.addParameter('revealSeconds', 6, @(x)isnumeric(x)&&isscalar(x)&&x>0);

% Horizontal
p.addParameter('horizDegrees', 360, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('horizSeconds', 6, @(x)isnumeric(x)&&isscalar(x)&&x>0);

% Vertical
p.addParameter('vertRange', [-5 45], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('vertSeconds', 4, @(x)isnumeric(x)&&isscalar(x)&&x>0);

% Capture behavior
p.addParameter('includeLegend', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('hideSurfacesDuringReveal', false, @(x)islogical(x)&&isscalar(x));

p.parse(h, outDir, baseName, defaultView, varargin{:});
opt = p.Results;

fps  = opt.fps;
ease = lower(string(opt.ease));

% -------------------- Resolve fig & ax --------------------
[fig, axGuess] = localResolveFigAx(h);

if isempty(fig) || ~isgraphics(fig,'figure')
    error('No live figure found. Keep the figure open while generating the video.');
end

ax = opt.ax;
if isempty(ax), ax = axGuess; end
if isempty(ax) || ~isgraphics(ax)
    ax = gca;
end
if isempty(ax) || ~isgraphics(ax)
    error('Could not resolve axes. Pass ''ax'', gca explicitly.');
end

% -------------------- Configure figure (MUST happen BEFORE VideoWriter locks size) --------------------
set(fig, 'Renderer', char(opt.renderer));
set(fig, 'Units', 'pixels');
if ~isempty(opt.figPos)
    set(fig, 'Position', opt.figPos);
end
drawnow;  % freeze layout

% lock vis3d if possible
if opt.lockVis3d
    try, axis(ax,'vis3d'); catch, end
end

defaultView = opt.defaultView(:).';
localSetView(ax, defaultView);
drawnow;

% -------------------- Fix capture rectangle (prevents frame-size errors) --------------------
rect = localCaptureRect(fig, ax, opt.includeLegend);
drawnow;

% -------------------- Ensure output dir exists --------------------
outDir = char(opt.outDir);
if ~exist(outDir,'dir'), mkdir(outDir); end
baseName = char(opt.baseName);
outFile  = fullfile(outDir, baseName + ".mp4");

% -------------------- Collect reveal targets --------------------
G = localCollectRevealTargets(ax);

% Cache originals + hide everything revealable (blank start)
G = localCacheAndHide(G, opt.hideSurfacesDuringReveal);

% determine how many "sessions" exist (aligned columns)
maxK = localMaxValidSessions(G);
if maxK < 2
    warning('Reveal: maxK<2 (not enough session points found). Will still write rotations.');
end

% -------------------- Create writer (AFTER size is locked) --------------------
vw = VideoWriter(char(outFile), 'MPEG-4');
vw.FrameRate = fps;
open(vw);

% -------------------- PHASE 0: write an explicit BLANK start frame --------------------
localSetView(ax, defaultView);
localApplyRevealBySession(G, 0, opt.hideSurfacesDuringReveal); % k=0 => blank
drawnow;
writeVideo(vw, getframe(fig, rect));  % locks frame size reliably

% -------------------- PHASE 1: Reveal (session-by-session, gradual) --------------------
revealFrames = max(2, round(opt.revealSeconds * fps));
tR  = linspace(0, 1, revealFrames);
teR = localEase(tR, ease);
kSeq = floor(maxK * teR);        % goes 0..maxK gradually
kSeq = localCumMax(kSeq);        % monotonic

for ii = 1:numel(kSeq)
    localSetView(ax, defaultView);
    localApplyRevealBySession(G, kSeq(ii), opt.hideSurfacesDuringReveal);
    drawnow;
    writeVideo(vw, getframe(fig, rect));
end

% ensure fully revealed before rotations
localApplyRevealBySession(G, maxK, opt.hideSurfacesDuringReveal);
drawnow;
writeVideo(vw, getframe(fig, rect));

% -------------------- PHASE 2: Horizontal rotation --------------------
horizFrames = max(2, round(opt.horizSeconds * fps));
tH  = linspace(0,1,horizFrames);
teH = localEase(tH, ease);

az0 = defaultView(1);
el0 = defaultView(2);
azH = az0 + opt.horizDegrees * teH;
elH = el0 * ones(size(azH));

for i = 1:numel(azH)
    localSetView(ax, [azH(i) elH(i)]);
    drawnow;
    writeVideo(vw, getframe(fig, rect));
end

% -------------------- PHASE 3: Vertical rotation --------------------
vertFrames = max(2, round(opt.vertSeconds * fps));
tV  = linspace(0,1,vertFrames);
teV = localEase(tV, ease);

azV = azH(end) * ones(size(teV));
elV = opt.vertRange(1) + (opt.vertRange(2)-opt.vertRange(1))*teV;

for i = 1:numel(azV)
    localSetView(ax, [azV(i) elV(i)]);
    drawnow;
    writeVideo(vw, getframe(fig, rect));
end

close(vw);

% restore figure to original look
localRestore(G);

fprintf('Saved: %s\n', outFile);

end % ===== END MAIN =====


% ======================= Helpers =======================

function localSetView(ax, v)
try
    view(ax, v);
catch
    view(v);
end
end

function rect = localCaptureRect(fig, ax, includeLegend)
% Fixed capture rectangle to avoid VideoWriter frame size mismatch.
% If includeLegend=true, capture full figure area; else capture axes box.

set(fig, 'Units','pixels');
figPos = get(fig,'Position'); % [x y w h]

if includeLegend
    rect = [1 1 figPos(3) figPos(4)];
else
    try
        set(ax,'Units','pixels');
        ap = get(ax,'Position'); % within figure
        rect = round(ap);
    catch
        rect = [1 1 figPos(3) figPos(4)];
    end
end
rect = round(rect);
end

function te = localEase(t, ease)
switch lower(string(ease))
    case "cosine"
        te = 0.5 - 0.5*cos(pi*t);
    otherwise
        te = t;
end
end

function y = localCumMax(x)
y = x;
for ii = 2:numel(y)
    if y(ii) < y(ii-1)
        y(ii) = y(ii-1);
    end
end
end

function [fig, ax] = localResolveFigAx(h)
fig = []; ax = [];

if isstruct(h)
    if isfield(h,'base')
        if isfield(h.base,'fig') && isgraphics(h.base.fig,'figure'), fig = h.base.fig; end
        if isfield(h.base,'ax')  && isgraphics(h.base.ax),           ax  = h.base.ax;  end
    end
    if isempty(fig) && isfield(h,'fig') && isgraphics(h.fig,'figure'), fig = h.fig; end
    if isempty(ax)  && isfield(h,'ax')  && isgraphics(h.ax),          ax  = h.ax;  end
end

if isempty(fig) || ~isgraphics(fig,'figure')
    fig = gcf;
end
end

function tf = localHasXYZ(h)
tf = false;
try
    xd = get(h,'XData'); yd = get(h,'YData'); zd = get(h,'ZData');
    tf = isnumeric(xd) && isnumeric(yd) && isnumeric(zd) && ~isempty(xd) && ~isempty(yd) && ~isempty(zd);
catch
end
end

function G = localCollectRevealTargets(ax)
% We reveal any object in the axes that has X/Y/Z:
% - line (including marker-only dots)
% - scatter
% We do NOT partially reveal surface/patch; those are handled by Visible toggle.
objs = findall(ax);

G = struct('h',{},'type',{},'X0',{},'Y0',{},'Z0',{},'validIdx',{},'nValid',{},'Vis0',{},'isVector',{});
idx = 0;

for i = 1:numel(objs)
    hi = objs(i);
    if ~isgraphics(hi), continue; end
    if isequal(hi, ax), continue; end

    % Surfaces/patches might have X/Y/Z too; we will store them for visibility toggles
    t = "";
    try, t = string(get(hi,'Type')); catch, end
    t = lower(t);

    if any(t == ["line","scatter","surface","patch"]) && localHasXYZ(hi)
        idx = idx+1;
        G(idx).h = hi;
        G(idx).type = char(t);
    end
end
end

function G = localCacheAndHide(G, hideSurfacesDuringReveal)
for i = 1:numel(G)
    hi = G(i).h;

    % cache visibility
    try, G(i).Vis0 = get(hi,'Visible'); catch, G(i).Vis0 = 'on'; end

    % surface/patch: hide/show only
    if any(strcmp(G(i).type, {'surface','patch'}))
        if hideSurfacesDuringReveal
            try, set(hi,'Visible','off'); catch, end
        end
        continue;
    end

    X0 = get(hi,'XData'); Y0 = get(hi,'YData'); Z0 = get(hi,'ZData');

    % only handle vectors for session-by-session reveal; otherwise hide by Visible
    if ~(isvector(X0) && isvector(Y0) && isvector(Z0))
        G(i).isVector = false;
        G(i).X0 = X0; G(i).Y0 = Y0; G(i).Z0 = Z0;
        try, set(hi,'Visible','off'); catch, end
        continue;
    end

    G(i).isVector = true;
    X0 = X0(:).'; Y0 = Y0(:).'; Z0 = Z0(:).';
    G(i).X0 = X0; G(i).Y0 = Y0; G(i).Z0 = Z0;

    validMask = isfinite(X0) & isfinite(Y0) & isfinite(Z0);
    G(i).validIdx = find(validMask);
    G(i).nValid   = numel(G(i).validIdx);

    % hide everything by NaN masking
    Xh = X0; Yh = Y0; Zh = Z0;
    Xh(G(i).validIdx) = NaN;
    Yh(G(i).validIdx) = NaN;
    Zh(G(i).validIdx) = NaN;

    set(hi,'XData',Xh,'YData',Yh,'ZData',Zh,'Visible',G(i).Vis0);
end
end

function maxK = localMaxValidSessions(G)
maxK = 0;
for i = 1:numel(G)
    if isfield(G(i),'nValid') && ~isempty(G(i).nValid)
        maxK = max(maxK, G(i).nValid);
    end
end
end

function localApplyRevealBySession(G, k, hideSurfacesDuringReveal)
% k is session index; k=0 => blank.
for i = 1:numel(G)
    hi = G(i).h;
    if ~isgraphics(hi), continue; end

    % surfaces/patch: show only once reveal begins (optional)
    if any(strcmp(G(i).type, {'surface','patch'}))
        if hideSurfacesDuringReveal
            if k==0
                try, set(hi,'Visible','off'); catch, end
            else
                try, set(hi,'Visible',G(i).Vis0); catch, end
            end
        end
        continue;
    end

    if ~isfield(G(i),'isVector') || ~G(i).isVector, continue; end
    if ~isfield(G(i),'nValid') || G(i).nValid==0, continue; end

    kk = max(0, min(k, G(i).nValid));
    vidx = G(i).validIdx;

    X = G(i).X0; Y = G(i).Y0; Z = G(i).Z0;
    X(vidx) = NaN; Y(vidx) = NaN; Z(vidx) = NaN;

    if kk > 0
        showIdx = vidx(1:kk);
        X(showIdx) = G(i).X0(showIdx);
        Y(showIdx) = G(i).Y0(showIdx);
        Z(showIdx) = G(i).Z0(showIdx);
    end

    set(hi,'XData',X,'YData',Y,'ZData',Z,'Visible',G(i).Vis0);
end
end

function localRestore(G)
for i = 1:numel(G)
    hi = G(i).h;
    if ~isgraphics(hi), continue; end
    if isfield(G(i),'X0')
        try
            set(hi,'XData',G(i).X0,'YData',G(i).Y0,'ZData',G(i).Z0);
        catch
        end
    end
    if isfield(G(i),'Vis0')
        try, set(hi,'Visible',G(i).Vis0); catch, end
    end
end
end
