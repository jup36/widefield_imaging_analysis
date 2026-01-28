function out = makeMDSRevealAllAnimals_thenRotate_fromData( ...
    Y, sessInfo, cmap, ellipsoid_full, outDir, baseName, varargin)
% makeMDSRevealAllAnimals_thenRotate_fromData
%
% Reveal sessions ACROSS ALL animals together (session index 1..max),
% with ALL trajectory line segments shown from the beginning (graded faint->solid; fixed),
% then do a horizontal rotate (starting at defaultView),
% then reset to defaultView and do a tilt (starting at defaultView).
%
% Key features kept:
%   - All lines appear from the beginning (graded faint->solid; fixed)
%   - Points reveal over time, fixed faint->solid (early stays faint)
%   - Ellipsoid overlay + mu marker
%   - Black edge for points inside ellipsoid
%   - Legend overlaid INSIDE axes (line handles only)
%   - Title removed
%   - NO live-window clipping: pixel-margin axes placement every frame
%   - SaveMethod supports: 'print' | 'exportgraphics' | 'getframe'
%
% REQUIRED:
%   Y              : [N x >=3] coords
%   sessInfo       : table/struct with mouseId or parsable header + ordering
%   cmap           : [nAnimals x 3]
%   ellipsoid_full : struct mu,Sigma (optional confLevel) or []
%   outDir, baseName
%
% OPTIONS (Name/Value):
%   'defaultView'     : [az el] REQUIRED
%   'dims'            : [1 2 3]
%   'orderField'      : 'sessWithin' etc (recommended)
%
%   'fps'             : 30
%   'revealSeconds'   : 6
%   'horizSeconds'    : 4
%   'tiltSeconds'     : 3
%   'horizAzRange'    : [0 90]  relative to default az
%   'tiltElRange'     : [0 45]  relative to default el
%
%   'figPos'          : [x y w h] default [100 100 1200 900]
%   'figWidthScale'   : 1.0
%   'figScale'        : 0.5
%   'renderer'        : 'opengl'
%
%   'imgExt'          : 'png'
%   'dpi'             : 150
%   'saveMethod'      : 'print' | 'exportgraphics' | 'getframe'
%
% Legend:
%   'legendMode'      : 'inside' | 'off'
%   'legendPos'       : [] or [x y w h] axes-normalized overlay
%   'legendLoc'       : fallback ('northeast')
%
% Styling:
%   'minAlpha','maxAlpha','minDotSize','maxDotSize'
%
% Axes layout (THIS FIXES LIVE CLIPPING):
%   'axesMarginPx'    : [L B R T] in pixels
%                       default = [140 90 40 40]
%                       (increase L/B if labels/ticks clip)
%
% OUTPUT:
%   out.frameDir, out.videoFile, out.nFrames

% -------------------- Parse(toggle) --------------------
p = inputParser;
p.addRequired('Y', @(x)isnumeric(x)&&size(x,1)>0);
p.addRequired('sessInfo');
p.addRequired('cmap', @(x)isnumeric(x)&&size(x,2)==3);
p.addRequired('ellipsoid_full');
p.addRequired('outDir', @(x)ischar(x)||isstring(x));
p.addRequired('baseName', @(x)ischar(x)||isstring(x));

p.addParameter('defaultView', [], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('dims', [1 2 3], @(x)isnumeric(x)&&numel(x)==3);
p.addParameter('orderField', '', @(s)ischar(s)||isstring(s));

p.addParameter('fps', 30, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('revealSeconds', 6, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('horizSeconds', 4, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('tiltSeconds', 3, @(x)isnumeric(x)&&isscalar(x)&&x>0);

p.addParameter('horizAzRange', [0 90], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('tiltElRange',  [0 45], @(x)isnumeric(x)&&numel(x)==2);

p.addParameter('figPos', [100 100 1200 900], @(x)isnumeric(x)&&numel(x)==4);
p.addParameter('figWidthScale', 1.0, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('figScale', 0.5, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('renderer', 'opengl', @(s)ischar(s)||isstring(s));

p.addParameter('imgExt', 'png', @(s)ischar(s)||isstring(s));
p.addParameter('dpi', 150, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('saveMethod', 'print', @(s)ischar(s)||isstring(s)); % print|exportgraphics|getframe

p.addParameter('legendMode', 'inside', @(s)ischar(s)||isstring(s));
p.addParameter('legendLoc',  'northeast', @(s)ischar(s)||isstring(s));
p.addParameter('legendPos',  [], @(x) isempty(x) || (isnumeric(x)&&numel(x)==4));

p.addParameter('minAlpha', 0.15, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('maxAlpha', 1.00, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('minDotSize', 40, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('maxDotSize', 90, @(x)isnumeric(x)&&isscalar(x)&&x>0);

% Ellipsoid
p.addParameter('ellAlpha', 0.12, @(v)isnumeric(v)&&isscalar(v)&&v>=0&&v<=1);
p.addParameter('ellEdgeAlpha', 0.05, @(v)isnumeric(v)&&isscalar(v)&&v>=0&&v<=1);
p.addParameter('ellN', 40, @(v)isnumeric(v)&&isscalar(v)&&v>=8);
p.addParameter('ellLineWidth', 0.5, @(v)isnumeric(v)&&isscalar(v)&&v>=0);
p.addParameter('ellColor', [0 0 0], @(v)isnumeric(v)&&numel(v)==3);
p.addParameter('highlightInside', true, @(v)islogical(v)||ismember(v,[0 1]));
p.addParameter('insideEdgeColor', 'k');
p.addParameter('insideLineWidth', 2, @(v)isnumeric(v)&&isscalar(v)&&v>=0);

% Pixel margins (L B R T) to prevent LIVE clipping
p.addParameter('axesMarginPx', [140 90 40 40], @(x)isnumeric(x)&&numel(x)==4);

p.parse(Y, sessInfo, cmap, ellipsoid_full, outDir, baseName, varargin{:});
opt = p.Results;

assert(~isempty(opt.defaultView), 'You must pass ''defaultView'', [az el].');

dims = opt.dims(:).';
defaultView = opt.defaultView(:).';
az0 = defaultView(1); el0 = defaultView(2);

outDir = char(opt.outDir);
if ~exist(outDir,'dir'), mkdir(outDir); end

imgExt = lower(char(opt.imgExt));
saveMethod = lower(string(opt.saveMethod));

% -------------------- Unique frame folder --------------------
frameDir = localMakeUniqueDir(outDir, [char(baseName) '_frames']);

% -------------------- Extract IDs + ordering --------------------
[mIdPerRow, orderKey] = localExtractIdAndOrder(sessInfo, opt.orderField);

allAnimals = unique(mIdPerRow, 'stable');
nAnimals = numel(allAnimals);

rowsByAnimal = cell(nAnimals,1);
for a = 1:nAnimals
    idx = find(strcmp(mIdPerRow, allAnimals{a}));
    [~,ord] = sort(orderKey(idx), 'ascend');
    rowsByAnimal{a} = idx(ord);
end
maxSess = max(cellfun(@numel, rowsByAnimal));
fprintf('[Reveal-all] animals=%d, maxSessAcrossAnimals=%d, totalRows=%d\n', ...
    nAnimals, maxSess, size(Y,1));

% -------------------- Figure / axes (NO LIVE CLIPPING) --------------------
fig = figure('Color','w', 'Renderer', char(opt.renderer));
pos = opt.figPos;
pos(3) = round(pos(3) * opt.figWidthScale);
pos(3:4) = round(pos(3:4) * opt.figScale);
set(fig, 'Position', pos);
set(fig, 'InvertHardcopy','off');

% remove figure UI chrome that can mess with available canvas
try, fig.MenuBar = 'none'; end %#ok<TRYNC>
try, fig.ToolBar = 'none'; end %#ok<TRYNC>

ax = axes('Parent',fig);
hold(ax,'on'); grid(ax,'on');
axis(ax,'vis3d');
view(ax, defaultView);

% kill the little axes toolbar overlay (top-right icons)
try, ax.Toolbar.Visible = 'off'; end %#ok<TRYNC>
try, disableDefaultInteractivity(ax); end %#ok<TRYNC>

xlabel(ax, sprintf('MDS%d', dims(1)));
ylabel(ax, sprintf('MDS%d', dims(2)));
zlabel(ax, sprintf('MDS%d', dims(3)));
title(ax, '');

camproj(ax,'perspective');
set(ax,'CameraViewAngleMode','manual');

drawnow;
localApplyAxesMarginPx(fig, ax, opt.axesMarginPx);
drawnow;

% -------------------- Ellipsoid overlay + inside test --------------------
hEll = struct('surf', gobjects(1), 'muMarker', gobjects(1));
insideMask = false(size(Y,1),1);

if ~isempty(ellipsoid_full) && isstruct(ellipsoid_full) ...
        && isfield(ellipsoid_full,'mu') && isfield(ellipsoid_full,'Sigma')
    [hEll, insideMask] = localPlotEllipsoidAndInside(ax, Y(:,dims), ellipsoid_full, opt);
end

% -------------------- Build handles --------------------
ptH  = cell(nAnimals,1);
lineProxyH = gobjects(nAnimals,1);

for a = 1:nAnimals
    idxList = rowsByAnimal{a};
    nA = numel(idxList);
    baseC = cmap(1+mod(a-1,size(cmap,1)),:);

    % Legend proxy line (line-only legend)
    lineProxyH(a) = plot3(ax, nan, nan, nan, '-', 'Color', baseC, 'LineWidth', 2);
    lineProxyH(a).HandleVisibility = 'on';
    lineProxyH(a).DisplayName = allAnimals{a};

    % Segments: all visible immediately (graded, fixed)
    for k = 1:max(nA-1,0)
        strength = localLerp(opt.minAlpha, opt.maxAlpha, (k+1)/max(2,nA));
        cSeg = localMixWithWhite(baseC, 1 - strength);
        r1 = idxList(k); r2 = idxList(k+1);
        h = plot3(ax, ...
            Y([r1 r2],dims(1)), Y([r1 r2],dims(2)), Y([r1 r2],dims(3)), ...
            '-', 'Color', cSeg, 'LineWidth', 2);
        h.HandleVisibility = 'off';
    end

    % Points: fixed alpha/color; start hidden; revealed over time
    ptH{a} = gobjects(nA,1);
    for s = 1:nA
        strength = localLerp(opt.minAlpha, opt.maxAlpha, s/max(1,nA));
        cPt = localMixWithWhite(baseC, 1 - strength);
        sz  = localLerp(opt.minDotSize, opt.maxDotSize, s/max(1,nA));
        r = idxList(s);

        ptH{a}(s) = scatter3(ax, Y(r,dims(1)), Y(r,dims(2)), Y(r,dims(3)), ...
            sz, cPt, 'filled', 'MarkerFaceAlpha', strength);
        ptH{a}(s).HandleVisibility = 'off';
        ptH{a}(s).Visible = 'off';

        if opt.highlightInside && insideMask(r)
            ptH{a}(s).MarkerEdgeColor = opt.insideEdgeColor;
            ptH{a}(s).LineWidth = opt.insideLineWidth;
            if isprop(ptH{a}(s), 'MarkerEdgeAlpha')
                ptH{a}(s).MarkerEdgeAlpha = 0.6;
            end
        else
            ptH{a}(s).MarkerEdgeColor = 'none';
            ptH{a}(s).LineWidth = 0.4;
        end
    end
end

% Ellipsoid behind
try
    if isgraphics(hEll.surf), uistack(hEll.surf,'bottom'); end
catch
end

% -------------------- Legend: overlay inside axes --------------------
lgd = gobjects(1);
legMode = lower(string(opt.legendMode));
if legMode ~= "off"
    lgd = legend(ax, lineProxyH, allAnimals, 'Interpreter','none');
    lgd.Box = 'off';
    lgd.Color = 'none';
    lgd.Units = 'normalized';

    if ~isempty(opt.legendPos)
        localPlaceLegendInsideAxes(ax, lgd, opt.legendPos);
    else
        lgd.Location = char(opt.legendLoc);
        drawnow;
        localClampLegendToAxes(ax, lgd, 0.01);
    end
end
drawnow;

% IMPORTANT: re-apply pixel margin after legend settles
localApplyAxesMarginPx(fig, ax, opt.axesMarginPx);
if isgraphics(lgd)
    if ~isempty(opt.legendPos), localPlaceLegendInsideAxes(ax, lgd, opt.legendPos);
    else,                      localClampLegendToAxes(ax, lgd, 0.01);
    end
end
drawnow;

% -------------------- Timing --------------------
fps = opt.fps;
revealFrames = max(1, round(opt.revealSeconds * fps));
hFrames      = max(1, round(opt.horizSeconds  * fps));
tFrames      = max(1, round(opt.tiltSeconds   * fps));

frameIdx = 0;

% -------------------- Phase 1: reveal --------------------
for f = 1:revealFrames
    frac = (f-1) / max(1,(revealFrames-1));
    sShow = floor(1 + frac*(maxSess-1));

    for a = 1:nAnimals
        nA = numel(ptH{a});
        nVis = min(sShow, nA);
        for s = 1:nA
            ptH{a}(s).Visible = ternary(s<=nVis, 'on', 'off');
        end
    end

    view(ax, defaultView);
    localPreSaveFix(fig, ax, lgd, opt.legendPos, opt.axesMarginPx);
    frameIdx = frameIdx + 1;
    localSaveFrame(fig, ax, frameDir, frameIdx, imgExt, opt.dpi, saveMethod);
end

% -------------------- Phase 2: horizontal rotation --------------------
view(ax, defaultView); drawnow;

azStart = az0 + opt.horizAzRange(1);
azEnd   = az0 + opt.horizAzRange(2);
az = linspace(azStart, azEnd, hFrames);
el = el0 * ones(size(az));

for i = 1:hFrames
    view(ax, [az(i) el(i)]);
    localPreSaveFix(fig, ax, lgd, opt.legendPos, opt.axesMarginPx);
    frameIdx = frameIdx + 1;
    localSaveFrame(fig, ax, frameDir, frameIdx, imgExt, opt.dpi, saveMethod);
end

% -------------------- Phase 3: tilt --------------------
view(ax, defaultView); drawnow;

elStart = el0 + opt.tiltElRange(1);
elEnd   = el0 + opt.tiltElRange(2);
el = linspace(elStart, elEnd, tFrames);
az = az0 * ones(size(el));

for i = 1:tFrames
    view(ax, [az(i) el(i)]);
    localPreSaveFix(fig, ax, lgd, opt.legendPos, opt.axesMarginPx);
    frameIdx = frameIdx + 1;
    localSaveFrame(fig, ax, frameDir, frameIdx, imgExt, opt.dpi, saveMethod);
end

% -------------------- Assemble video --------------------
profile = 'Motion JPEG AVI';  % 'MPEG-4' or 'Motion JPEG AVI'
quality = 100;

switch lower(profile)
    case 'mpeg-4'
        videoFile = fullfile(outDir, [char(baseName) '.mp4']);
    case 'motion jpeg avi'
        videoFile = fullfile(outDir, [char(baseName) '.avi']);
    otherwise
        error('Unsupported profile: %s', profile);
end

localWriteMp4FromFrames(frameDir, imgExt, videoFile, fps, ...
    'profile', profile, 'quality', quality);

% or:
% localWriteMp4FromFrames(frameDir, imgExt, videoFile, fps, 'profile','Motion JPEG AVI', 'quality', 100);


out.frameDir = frameDir;
out.videoFile = videoFile;
out.nFrames = frameIdx;

fprintf('Saved frames: %s\n', frameDir);
fprintf('Saved video : %s\n', videoFile);

end

% =====================================================================
%                               SUBFUNCTIONS
% =====================================================================

function frameDir = localMakeUniqueDir(parent, name)
frameDir = fullfile(parent, name);
if ~exist(frameDir,'dir')
    mkdir(frameDir);
    return;
end
k = 1;
while true
    cand = sprintf('%s_%02d', frameDir, k);
    if ~exist(cand,'dir')
        mkdir(cand);
        frameDir = cand;
        return;
    end
    k = k + 1;
end
end

function localApplyAxesMarginPx(fig, ax, marginPx)
% marginPx = [L B R T] in pixels
marginPx = double(marginPx(:).');
L = marginPx(1); B = marginPx(2); R = marginPx(3); T = marginPx(4);

fig.Units = 'pixels';
ax.Units  = 'pixels';
fp = fig.Position;  % [x y w h] in px
W = fp(3); H = fp(4);

% guard
W = max(W, 10); H = max(H, 10);
L = min(L, W-5); R = min(R, W-5);
B = min(B, H-5); T = min(T, H-5);

axPosPx = [L, B, max(10, W - L - R), max(10, H - B - T)];

% lock axes box to what we set (avoid auto outer-position meddling)
try
    ax.PositionConstraint = 'innerposition';
catch
end

ax.Position = axPosPx;
ax.Units = 'normalized'; % keep rest of code simple after this
end

function localPreSaveFix(fig, ax, lgd, legendPosAxes, axesMarginPx)
drawnow;
localApplyAxesMarginPx(fig, ax, axesMarginPx);

% Re-place legend after axes move
if isgraphics(lgd)
    if ~isempty(legendPosAxes)
        localPlaceLegendInsideAxes(ax, lgd, legendPosAxes);
    else
        localClampLegendToAxes(ax, lgd, 0.01);
    end
end

set(fig,'InvertHardcopy','off');
drawnow;
end

function localSaveFrame(fig, ax, frameDir, idx, ext, dpi, saveMethod)
fn = fullfile(frameDir, sprintf('frame_%06d.%s', idx, ext));
drawnow;
saveMethod = lower(string(saveMethod));

switch saveMethod
    case "getframe"
        % Captures EXACT on-screen content (best match to what you monitor)
        fr = getframe(fig);
        imwrite(fr.cdata, fn);

    case "exportgraphics"
        % Export figure canvas (not axes) so margins/legend overlay are respected
        try
            exportgraphics(fig, fn, 'Resolution', dpi, 'BackgroundColor','white');
        catch
            set(fig,'PaperPositionMode','auto');
            print(fig, fn, ['-d' ext], sprintf('-r%d', dpi));
        end

    otherwise % "print"
        set(fig,'PaperPositionMode','auto');
        print(fig, fn, ['-d' ext], sprintf('-r%d', dpi));
end
end

function localWriteMp4FromFrames(frameDir, ext, outFile, fps, varargin)
% localWriteMp4FromFrames(..., 'profile','MPEG-4', 'quality',100)

p = inputParser;
p.addParameter('profile', 'MPEG-4', @(s)ischar(s)||isstring(s));
p.addParameter('quality', 100, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=100);
p.parse(varargin{:});
profile = char(p.Results.profile);
quality = p.Results.quality;

files = dir(fullfile(frameDir, ['frame_*.' ext]));
assert(~isempty(files), 'No frame files found.');
[~,ord] = sort({files.name});
files = files(ord);

vw = VideoWriter(outFile, profile);
vw.FrameRate = fps;

% Quality is supported for MPEG-4 and Motion JPEG AVI in most MATLAB builds
try
    vw.Quality = quality;
catch
end

open(vw);

im0 = imread(fullfile(frameDir, files(1).name));
im0 = localPadToEven(im0);
targetH = size(im0,1);
targetW = size(im0,2);
writeVideo(vw, im0);

for i = 2:numel(files)
    im = imread(fullfile(frameDir, files(i).name));
    im = localPadToEven(im);
    if size(im,1) ~= targetH || size(im,2) ~= targetW
        im = imresize(im, [targetH targetW]);
    end
    writeVideo(vw, im);
end

close(vw);
end


function im = localPadToEven(im)
h = size(im,1); w = size(im,2);
newH = h + mod(h,2);
newW = w + mod(w,2);
if newH==h && newW==w, return; end
im2 = uint8(255*ones(newH, newW, size(im,3)));
im2(1:h,1:w,:) = im;
im = im2;
end

function [hEll, inside] = localPlotEllipsoidAndInside(ax, Y3, ell, opt)
hEll = struct('surf', gobjects(1), 'muMarker', gobjects(1));
inside = false(size(Y3,1),1);

mu = ell.mu(:)'; mu = mu(1:3);
Sigma = ell.Sigma(1:3,1:3);

confLevel = 0.95;
if isfield(ell,'confLevel') && ~isempty(ell.confLevel)
    confLevel = ell.confLevel;
end
thr = chi2inv(confLevel, 3);
scale = sqrt(thr);

Sigma = (Sigma + Sigma')/2;
Sigma = Sigma + 1e-8*eye(3);

[V,L] = eig(Sigma);
radii = scale * sqrt(max(diag(L), 0));

[xs, ys, zs] = sphere(opt.ellN);
U = [xs(:) ys(:) zs(:)]';
E = V * (diag(radii) * U);

Ex = reshape(E(1,:) + mu(1), size(xs));
Ey = reshape(E(2,:) + mu(2), size(ys));
Ez = reshape(E(3,:) + mu(3), size(zs));

hEll.surf = surf(ax, Ex, Ey, Ez, ...
    'FaceAlpha', opt.ellAlpha, ...
    'EdgeAlpha', opt.ellEdgeAlpha, ...
    'LineWidth', opt.ellLineWidth, ...
    'FaceColor', opt.ellColor, ...
    'EdgeColor', opt.ellColor);

hEll.muMarker = scatter3(ax, mu(1), mu(2), mu(3), 110, ...
    'Marker','p', ...
    'MarkerFaceColor', opt.ellColor, ...
    'MarkerEdgeColor', opt.ellColor);

% keep out of legend
try
    set(get(get(hEll.surf,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');
    set(get(get(hEll.muMarker,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');
catch
end

% inside test
X = Y3 - mu;
R = chol(Sigma);
Z = X / R;
md2 = sum(Z.^2, 2);
inside = (md2 <= thr);
end

function localPlaceLegendInsideAxes(ax, lgd, posAxesNorm)
% posAxesNorm is [x y w h] in axes-normalized coords
ax.Units = 'normalized';
P = ax.Position; % figure-normalized
posAxesNorm = double(posAxesNorm(:).');

posFig = [P(1)+posAxesNorm(1)*P(3), ...
          P(2)+posAxesNorm(2)*P(4), ...
          posAxesNorm(3)*P(3), ...
          posAxesNorm(4)*P(4)];
lgd.Location = 'none';
lgd.Units = 'normalized';
lgd.Position = posFig;
end

function localClampLegendToAxes(ax, lgd, pad)
if nargin < 3, pad = 0.01; end
ax.Units = 'normalized';
P = ax.Position;
lgd.Units = 'normalized';
L = lgd.Position;
L(1) = max(L(1), P(1)+pad);
L(2) = max(L(2), P(2)+pad);
L(1) = min(L(1), P(1)+P(3)-L(3)-pad);
L(2) = min(L(2), P(2)+P(4)-L(4)-pad);
lgd.Position = L;
lgd.Location = 'none';
end

function c = localMixWithWhite(baseC, whitenAmount)
whitenAmount = max(0,min(1,whitenAmount));
c = baseC*(1-whitenAmount) + 1.0*whitenAmount;
end

function y = localLerp(a,b,t)
t = max(0,min(1,t));
y = a + (b-a)*t;
end

function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end

function [mId, orderKey] = localExtractIdAndOrder(sessInfo, orderField)
if istable(sessInfo)
    vars = sessInfo.Properties.VariableNames;
    getcol = @(name) sessInfo.(name);
    n = height(sessInfo);
else
    vars = fieldnames(sessInfo);
    getcol = @(name) sessInfo.(name);
    n = numel(getcol(vars{1}));
end

idField = localPickFirstCI(vars, {'mouseId','mId','mouse','animal','animalId','mID'});
if ~isempty(idField)
    mId = cellstr(string(getcol(idField)));
else
    hField = localPickFirstCI(vars, {'header'});
    assert(~isempty(hField), 'sessInfo needs mouseId or header to parse mouseId.');
    hdr = string(getcol(hField));
    mId = cell(n,1);
    for ii = 1:n
        tok = regexp(hdr(ii), '(m\\d{3,5})', 'tokens','once');
        if isempty(tok)
            mId{ii} = sprintf('mouse_%03d', ii);
        else
            mId{ii} = tok{1};
        end
    end
end

orderField = char(string(orderField));
if ~isempty(orderField)
    assert(any(strcmp(vars, orderField)), 'orderField "%s" not found in sessInfo.', orderField);
    orderKey = localToNumericOrderKey(getcol(orderField));
    return;
end

cand = {'dateNum','datenum','date','sessDate','sessionDate', ...
        'sessWithin','sess','sessIdx','sessionIdx','trainingDay','dayNum','dayIndex'};
f = localPickFirstCI(vars, cand);
if ~isempty(f)
    orderKey = localToNumericOrderKey(getcol(f));
    return;
end

hField = localPickFirstCI(vars, {'header'});
assert(~isempty(hField), 'sessInfo needs an order field (e.g., sessWithin) or header containing date token.');
hdr = string(getcol(hField));
orderKey = nan(n,1);
for ii = 1:n
    tok = regexp(hdr(ii), '_(\\d{6})', 'tokens','once');
    if isempty(tok)
        orderKey(ii) = ii;
    else
        d6 = char(tok{1}); % assume MMDDYY
        mm = str2double(d6(1:2));
        dd = str2double(d6(3:4));
        yy = str2double(d6(5:6)) + 2000;
        orderKey(ii) = datenum(yy,mm,dd); %#ok<DATNM>
    end
end
end

function v = localToNumericOrderKey(x)
if isnumeric(x), v = double(x); return; end
if isdatetime(x), v = datenum(x); return; end %#ok<DATNM>
xs = string(x);
v = double(str2double(xs));
if any(~isfinite(v))
    try
        v = datenum(xs); %#ok<DATNM>
    catch
        error('Order key contains non-numeric / non-parseable values.');
    end
end
end

function f = localPickFirstCI(varNames, candidates)
f = '';
vLower = lower(string(varNames));
for i = 1:numel(candidates)
    c = lower(string(candidates{i}));
    hit = find(vLower == c, 1, 'first');
    if ~isempty(hit)
        f = varNames{hit};
        return;
    end
end
end
