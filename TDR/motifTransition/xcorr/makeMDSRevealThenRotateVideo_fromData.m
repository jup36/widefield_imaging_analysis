function outFile = makeMDSRevealThenRotateVideo_fromData( ...
    Y_full, sessInfo, cmap, ellipsoid_full, outDir, baseName, defaultView, varargin)
% makeMDSRevealThenRotateVideo_fromData
%
% Rebuilds the MDS+ellipsoid plot internally and generates a video:
%   0) Blank first frame
%   1) Reveal session-by-session at defaultView
%   2) One horizontal rotation
%   3) One vertical rotation
%
% REQUIRED:
%   Y_full, sessInfo, cmap, ellipsoid_full : inputs to plotMDSWithEllipsoid
%   outDir, baseName
%   defaultView : [az el]
%
% EXAMPLE:
% outFile = makeMDSRevealThenRotateVideo_fromData( ...
%   Y_full, sessInfo, cmap9, ellipsoid_full, figSaveDir, ...
%   'xcorr_mds_reveal_then_rotate_bothTrials', [-111 28], ...
%   'revealSeconds', 6, 'horizSeconds', 6, 'horizDegrees', 360, ...
%   'vertSeconds', 4, 'vertRange', [-5 45], ...
%   'fps', 30, 'renderer','opengl', 'figPos',[100 100 1200 900]);

% -------------------- Parse inputs --------------------
p = inputParser;
p.addRequired('Y_full');
p.addRequired('sessInfo');
p.addRequired('cmap');
p.addRequired('ellipsoid_full');
p.addRequired('outDir', @(x)ischar(x)||isstring(x));
p.addRequired('baseName', @(x)ischar(x)||isstring(x));
p.addRequired('defaultView', @(x)isnumeric(x)&&numel(x)==2);

% Pass-through to plotMDSWithEllipsoid
p.addParameter('dims', [1 2 3], @(x)isnumeric(x)&&numel(x)==3);
p.addParameter('view0', [35 20], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('ellAlpha', 0.10, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('ellN', 50, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('minAlpha', 0.25, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('maxAlpha', 1.0, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('minDotSize', 24, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('dotScaleMax', 3.0, @(x)isnumeric(x)&&isscalar(x));

% Video params
p.addParameter('fps', 30, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('renderer', 'opengl', @(s)ischar(s)||isstring(s));
p.addParameter('figPos', [100 100 1200 900], @(x)isnumeric(x)&&numel(x)==4);
p.addParameter('lockVis3d', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('ease', 'cosine', @(s)ischar(s)||isstring(s));

% Reveal / rotate
p.addParameter('revealSeconds', 6, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('horizDegrees', 360, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('horizSeconds', 6, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('vertRange', [-5 45], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('vertSeconds', 4, @(x)isnumeric(x)&&isscalar(x)&&x>0);

% Behavior toggles
p.addParameter('hideSurfacesDuringReveal', false, @(x)islogical(x)&&isscalar(x));

p.parse(Y_full, sessInfo, cmap, ellipsoid_full, outDir, baseName, defaultView, varargin{:});
opt = p.Results;

fps  = opt.fps;
ease = lower(string(opt.ease));
defaultView = opt.defaultView(:).';

% -------------------- Ensure output dir exists --------------------
outDir = char(opt.outDir);
if ~exist(outDir,'dir'), mkdir(outDir); end
baseName = char(opt.baseName);
outFile  = fullfile(outDir, baseName + ".mp4");

% -------------------- Build figure internally --------------------
fig = figure('Color','w');
set(fig,'Renderer',char(opt.renderer));
set(fig,'Units','pixels');
set(fig,'Position',opt.figPos);

% Make plot
h = plotMDSWithEllipsoid(Y_full, sessInfo, opt.cmap, opt.ellipsoid_full, ...
    'view',        opt.view0, ...
    'ellAlpha',    opt.ellAlpha, 'ellN', opt.ellN, ...
    'minAlpha',    opt.minAlpha, 'maxAlpha', opt.maxAlpha, ...
    'minDotSize',  opt.minDotSize, 'dotScaleMax', opt.dotScaleMax, ...
    'dims',        opt.dims);

drawnow;

% Resolve axes reliably
ax = gca;
try
    if isstruct(h) && isfield(h,'base') && isfield(h.base,'ax') && isgraphics(h.base.ax,'axes')
        ax = h.base.ax;
    end
catch
end

if opt.lockVis3d
    try, axis(ax,'vis3d'); catch, end
end

% Lock to desired default view for the animation
view(ax, defaultView);
drawnow;

% -------------------- Collect reveal targets --------------------
G = localCollectRevealTargets(ax);
G = localCacheAndHide(G, opt.hideSurfacesDuringReveal);

maxK = localMaxValidSessions(G);
if maxK < 2
    warning('Reveal: found maxK<2. Reveal may be trivial; rotations will still run.');
end

% -------------------- VideoWriter --------------------
vw = VideoWriter(char(outFile), 'MPEG-4');
vw.FrameRate = fps;
open(vw);

% -------------------- PHASE 0: blank first frame --------------------
view(ax, defaultView);
localApplyRevealBySession(G, 0, opt.hideSurfacesDuringReveal);
drawnow;
writeVideo(vw, getframe(ax));   % axes capture avoids toolbar rect warnings

% -------------------- PHASE 1: reveal --------------------
revealFrames = max(2, round(opt.revealSeconds * fps));
tR  = linspace(0, 1, revealFrames);
teR = localEase(tR, ease);
kSeq = floor(maxK * teR);
kSeq = localCumMax(kSeq);

for ii = 1:numel(kSeq)
    view(ax, defaultView);
    localApplyRevealBySession(G, kSeq(ii), opt.hideSurfacesDuringReveal);
    drawnow;
    writeVideo(vw, getframe(ax));
end

% ensure fully revealed
localApplyRevealBySession(G, maxK, opt.hideSurfacesDuringReveal);
drawnow;
writeVideo(vw, getframe(ax));

% -------------------- PHASE 2: horizontal rotation --------------------
horizFrames = max(2, round(opt.horizSeconds * fps));
tH  = linspace(0, 1, horizFrames);
teH = localEase(tH, ease);
az0 = defaultView(1); el0 = defaultView(2);
azH = az0 + opt.horizDegrees * teH;
elH = el0 * ones(size(azH));

for i = 1:numel(azH)
    view(ax, [azH(i) elH(i)]);
    drawnow;
    writeVideo(vw, getframe(ax));
end

% -------------------- PHASE 3: vertical rotation --------------------
vertFrames = max(2, round(opt.vertSeconds * fps));
tV  = linspace(0, 1, vertFrames);
teV = localEase(tV, ease);
azV = azH(end) * ones(size(teV));
elV = opt.vertRange(1) + (opt.vertRange(2)-opt.vertRange(1))*teV;

for i = 1:numel(azV)
    view(ax, [azV(i) elV(i)]);
    drawnow;
    writeVideo(vw, getframe(ax));
end

close(vw);

% Restore data (nice to leave the figure intact)
localRestore(G);

fprintf('Saved: %s\n', outFile);

end % ===== end main =====


% ======================= locals =======================

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
    if y(ii) < y(ii-1), y(ii) = y(ii-1); end
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
% Collect common 3D graphics we want to reveal:
% line (including marker-only dots), scatter.
% Also include surface/patch for optional hide/show.
objs = findall(ax);

G = struct('h',{},'type',{},'X0',{},'Y0',{},'Z0',{},'validIdx',{},'nValid',{},'Vis0',{},'isVector',{});
idx = 0;

for i = 1:numel(objs)
    hi = objs(i);
    if ~isgraphics(hi), continue; end
    if isequal(hi, ax), continue; end

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

    try, G(i).Vis0 = get(hi,'Visible'); catch, G(i).Vis0 = 'on'; end

    if any(strcmp(G(i).type, {'surface','patch'}))
        if hideSurfacesDuringReveal
            try, set(hi,'Visible','off'); catch, end
        end
        continue;
    end

    X0 = get(hi,'XData'); Y0 = get(hi,'YData'); Z0 = get(hi,'ZData');

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
for i = 1:numel(G)
    hi = G(i).h;
    if ~isgraphics(hi), continue; end

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
