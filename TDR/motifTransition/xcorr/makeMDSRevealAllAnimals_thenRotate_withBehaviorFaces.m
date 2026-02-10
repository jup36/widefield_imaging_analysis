function out = makeMDSRevealAllAnimals_thenRotate_withBehaviorFaces( ...
    Y, sessInfo, cmap, ellipsoid_full, behTbl, outDir, baseName, varargin)
% makeMDSRevealAllAnimals_thenRotate_withBehaviorFaces
%
% Reveal sessions ACROSS ALL animals together (session index 1..max),
% with ALL trajectory line segments shown from the beginning (graded faint->solid; fixed),
% then do a horizontal rotate (starting at defaultView),
% then reset to defaultView and do a tilt (starting at defaultView).
%
% NEW vs your previous function:
%   - Point FACE color encodes a behavioral metric (e.g., biasC).
%   - Point size scales with metric magnitude (or raw).
%   - NaN metric sessions are OPEN circles (no face fill).
%   - OPEN circles obey the SAME per-session alpha ramp (via MarkerEdgeAlpha).
%   - NO colorbar is drawn.
%
% OUTPUT:
%   out.frameDir, out.videoFile, out.nFrames

% -------------------- Parse --------------------
p = inputParser;
p.addRequired('Y', @(x)isnumeric(x)&&size(x,1)>0);
p.addRequired('sessInfo');
p.addRequired('cmap', @(x)isnumeric(x)&&size(x,2)==3);
p.addRequired('ellipsoid_full');
p.addRequired('behTbl', @(t) istable(t));
p.addRequired('outDir', @(x)ischar(x)||isstring(x));
p.addRequired('baseName', @(x)ischar(x)||isstring(x));

p.addParameter('defaultView', [], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('dims', [1 2 3], @(x)isnumeric(x)&&numel(x)==3);
p.addParameter('orderField', '', @(s)ischar(s)||isstring(s));

% behavior
p.addParameter('metric', 'biasC', @(s)ischar(s)||isstring(s));
p.addParameter('cLim', [], @(x) isempty(x) || (isnumeric(x)&&numel(x)==2));
p.addParameter('colormap', 'parula');
p.addParameter('flipColormap', true, @(x)islogical(x)&&isscalar(x));

p.addParameter('sizeBy', 'abs', @(s) any(strcmpi(string(s), ["abs","raw","none"])));
p.addParameter('baseSize', 40, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('sizeScale', 120, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('openIfNaN', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('openLineWidth', 1.5, @(x)isnumeric(x)&&isscalar(x)&&x>0);

% video timing
p.addParameter('fps', 30, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('revealSeconds', 6, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('horizSeconds', 4, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('tiltSeconds', 3, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('horizAzRange', [0 90], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('tiltElRange',  [0 45], @(x)isnumeric(x)&&numel(x)==2);

% figure/render/save
p.addParameter('figPos', [100 100 1200 900], @(x)isnumeric(x)&&numel(x)==4);
p.addParameter('figWidthScale', 1.0, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('figScale', 0.5, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('renderer', 'opengl', @(s)ischar(s)||isstring(s));
p.addParameter('imgExt', 'png', @(s)ischar(s)||isstring(s));
p.addParameter('dpi', 150, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('saveMethod', 'print', @(s)ischar(s)||isstring(s)); % print|exportgraphics|getframe

% legend
p.addParameter('legendMode', 'inside', @(s)ischar(s)||isstring(s));
p.addParameter('legendLoc',  'northeast', @(s)ischar(s)||isstring(s));
p.addParameter('legendPos',  [], @(x) isempty(x) || (isnumeric(x)&&numel(x)==4));

% styling
p.addParameter('minAlpha', 0.15, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('maxAlpha', 1.00, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('minDotSize', 40, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('maxDotSize', 90, @(x)isnumeric(x)&&isscalar(x)&&x>0);

% ellipsoid
p.addParameter('ellAlpha', 0.12, @(v)isnumeric(v)&&isscalar(v)&&v>=0&&v<=1);
p.addParameter('ellEdgeAlpha', 0.05, @(v)isnumeric(v)&&isscalar(v)&&v>=0&&v<=1);
p.addParameter('ellN', 40, @(v)isnumeric(v)&&isscalar(v)&&v>=8);
p.addParameter('ellLineWidth', 0.5, @(v)isnumeric(v)&&isscalar(v)&&v>=0);
p.addParameter('ellColor', [0 0 0], @(v)isnumeric(v)&&numel(v)==3);
p.addParameter('highlightInside', true, @(v)islogical(v)||ismember(v,[0 1]));
p.addParameter('insideEdgeColor', 'k');
p.addParameter('insideLineWidth', 2, @(v)isnumeric(v)&&isscalar(v)&&v>=0);

% axes margins
p.addParameter('axesMarginPx', [140 90 40 40], @(x)isnumeric(x)&&numel(x)==4);

p.parse(Y, sessInfo, cmap, ellipsoid_full, behTbl, outDir, baseName, varargin{:});
opt = p.Results;

assert(~isempty(opt.defaultView), 'You must pass ''defaultView'', [az el].');

dims = opt.dims(:).';
defaultView = opt.defaultView(:).';
az0 = defaultView(1); el0 = defaultView(2);

outDir = char(opt.outDir);
if ~exist(outDir,'dir'), mkdir(outDir); end

imgExt = lower(char(opt.imgExt));
saveMethod = lower(string(opt.saveMethod));

% -------------------- Align behavior metric to Y/sessInfo rows --------------------
assert(ismember('header', behTbl.Properties.VariableNames), 'behTbl must contain variable ''header''.');
hdrSess = string(localGetHeader(sessInfo));
[tf, loc] = ismember(hdrSess, string(behTbl.header));

metricName = localResolveMetricName(opt.metric, behTbl);

metricVal = nan(numel(hdrSess),1);
metricVal(tf) = behTbl{loc(tf), metricName};

% Determine CLim
if isempty(opt.cLim)
    mx = max(abs(metricVal), [], 'omitnan');
    if ~isfinite(mx) || mx==0, mx = 1; end
    cLim = [-mx mx];
else
    cLim = opt.cLim;
end

% Point sizes from metric
sz = localMetricToSize(metricVal, opt.baseSize, opt.sizeScale, opt.sizeBy);
sz = min(max(sz, opt.minDotSize), opt.maxDotSize);

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
fprintf('[Reveal-all + behFaces] animals=%d, maxSessAcrossAnimals=%d, totalRows=%d\n', ...
    nAnimals, maxSess, size(Y,1));

% -------------------- Figure / axes --------------------
fig = figure('Color','w', 'Renderer', char(opt.renderer));
pos = opt.figPos;
pos(3) = round(pos(3) * opt.figWidthScale);
pos(3:4) = round(pos(3:4) * opt.figScale);
set(fig, 'Position', pos);
set(fig, 'InvertHardcopy','off');
try, fig.MenuBar = 'none'; end %#ok<TRYNC>
try, fig.ToolBar = 'none'; end %#ok<TRYNC>

ax = axes('Parent',fig);
hold(ax,'on'); grid(ax,'on');
axis(ax,'vis3d');
view(ax, defaultView);

try, ax.Toolbar.Visible = 'off'; end %#ok<TRYNC>
try, disableDefaultInteractivity(ax); end %#ok<TRYNC>

xlabel(ax, sprintf('MDS%d', dims(1)));
ylabel(ax, sprintf('MDS%d', dims(2)));
zlabel(ax, sprintf('MDS%d', dims(3)));
title(ax, '');

camproj(ax,'perspective');
set(ax,'CameraViewAngleMode','manual');

% Metric colormap (NO colorbar)
cmapMetric = feval(char(opt.colormap), 256);
if opt.flipColormap, cmapMetric = flipud(cmapMetric); end
colormap(ax, cmapMetric);
caxis(ax, cLim);

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

    % Legend proxy
    lineProxyH(a) = plot3(ax, nan, nan, nan, '-', 'Color', baseC, 'LineWidth', 2);
    lineProxyH(a).HandleVisibility = 'on';
    lineProxyH(a).DisplayName = allAnimals{a};

    % Segments: all visible immediately (graded)
    for k = 1:max(nA-1,0)
        strength = localLerp(opt.minAlpha, opt.maxAlpha, (k+1)/max(2,nA));
        cSeg = localMixWithWhite(baseC, 1 - strength);
        r1 = idxList(k); r2 = idxList(k+1);
        h = plot3(ax, ...
            Y([r1 r2],dims(1)), Y([r1 r2],dims(2)), Y([r1 r2],dims(3)), ...
            '-', 'Color', cSeg, 'LineWidth', 2);
        h.HandleVisibility = 'off';
    end

    % Points: hidden initially; revealed over time
    ptH{a} = gobjects(nA,1);
    for s = 1:nA
        r = idxList(s);

        % session-dependent "faint -> solid" strength
        strength = localLerp(opt.minAlpha, opt.maxAlpha, s/max(1,nA));

        m  = metricVal(r);
        sZ = sz(r);

        if (~isfinite(m)) && opt.openIfNaN
            % =========================
            % OPEN circle (NaN metric)
            %   - edge color = animal baseC
            %   - edge alpha = strength   <-- THIS IS THE KEY FIX
            %   - no face
            % IMPORTANT: pass a real C argument (baseC) to avoid MATLAB parsing pitfalls
            % =========================
            ptH{a}(s) = scatter3(ax, ...
                Y(r,dims(1)), Y(r,dims(2)), Y(r,dims(3)), ...
                sZ, baseC, 'o');  % C argument set to RGB

            % enforce open style
            ptH{a}(s).MarkerFaceColor = 'none';
            ptH{a}(s).MarkerEdgeColor = baseC;
            ptH{a}(s).LineWidth = opt.openLineWidth;

            % edge alpha ramp (makes early sessions faint)
            if isprop(ptH{a}(s), 'MarkerEdgeAlpha')
                ptH{a}(s).MarkerEdgeAlpha = strength;
            end

        else
            % FILLED: face color encodes metric via colormap/caxis
            ptH{a}(s) = scatter3(ax, ...
                Y(r,dims(1)), Y(r,dims(2)), Y(r,dims(3)), ...
                sZ, m, 'filled'); % CData = m

            ptH{a}(s).MarkerFaceColor  = 'flat';
            ptH{a}(s).MarkerFaceAlpha  = strength; % session gradation
            ptH{a}(s).MarkerEdgeColor  = baseC;    % mouse identity on edge
            ptH{a}(s).LineWidth        = 0.75;
        end

        % Ellipsoid-inside edge emphasis (overrides if requested)
        if isgraphics(ptH{a}(s)) && opt.highlightInside && insideMask(r)
            ptH{a}(s).MarkerEdgeColor = opt.insideEdgeColor;
            ptH{a}(s).LineWidth = opt.insideLineWidth;
            if isprop(ptH{a}(s), 'MarkerEdgeAlpha')
                ptH{a}(s).MarkerEdgeAlpha = 0.6;
            end
        end

        ptH{a}(s).HandleVisibility = 'off';
        ptH{a}(s).Visible = 'off';
    end
end

% Ellipsoid behind
try
    if isgraphics(hEll.surf), uistack(hEll.surf,'bottom'); end
catch
end

% -------------------- Legend inside axes --------------------
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

% re-apply pixel margin after legend settles
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

out.frameDir = frameDir;
out.videoFile = videoFile;
out.nFrames = frameIdx;

fprintf('Saved frames: %s\n', frameDir);
fprintf('Saved video : %s\n', videoFile);

end

% =====================================================================
%                               HELPERS
% =====================================================================

function hdr = localGetHeader(sessInfo)
if istable(sessInfo)
    if ismember('header', sessInfo.Properties.VariableNames)
        hdr = sessInfo.header(:);
        return;
    end
else
    if isfield(sessInfo, 'header')
        hdr = sessInfo.header(:);
        return;
    end
end
error('sessInfo must contain a ''header'' field/variable for behavior alignment.');
end

function metricName = localResolveMetricName(metricReq, behTbl)
metricReq = lower(strtrim(string(metricReq)));
varNames  = string(behTbl.Properties.VariableNames);
varNamesL = lower(varNames);

if any(metricReq == ["bias","criterion","c"])
    metricReq = "biasc";
end

if any(metricReq == ["d'","dprime","dprm","d_prm"])
    cand = ["dprime","dprm","d_prm","d"];
    hit = cand(ismember(cand, varNamesL));
    if ~isempty(hit), metricReq = hit(1); end
end

if ~ismember(metricReq, varNamesL)
    error('Requested metric "%s" not found in behTbl. Available: %s', ...
        char(metricReq), strjoin(varNames, ", "));
end

metricName = varNames(varNamesL == metricReq);
metricName = metricName(1);
end

function sz = localMetricToSize(metricVal, baseSize, sizeScale, sizeBy)
modeSize = lower(string(sizeBy));
switch modeSize
    case "none"
        sz = baseSize * ones(size(metricVal));
    case "raw"
        m0 = metricVal;
        mx = max(abs(m0), [], 'omitnan'); if ~isfinite(mx) || mx==0, mx=1; end
        sz = baseSize + sizeScale * (m0./mx);
    otherwise % "abs"
        m0 = abs(metricVal);
        mx = max(m0, [], 'omitnan'); if ~isfinite(mx) || mx==0, mx=1; end
        sz = baseSize + sizeScale * (m0./mx);
end
sz(~isfinite(sz)) = baseSize;
end

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
marginPx = double(marginPx(:).');
L = marginPx(1); B = marginPx(2); R = marginPx(3); T = marginPx(4);

fig.Units = 'pixels';
ax.Units  = 'pixels';
fp = fig.Position;
W = fp(3); H = fp(4);

W = max(W, 10); H = max(H, 10);
L = min(L, W-5); R = min(R, W-5);
B = min(B, H-5); T = min(T, H-5);

axPosPx = [L, B, max(10, W - L - R), max(10, H - B - T)];

try, ax.PositionConstraint = 'innerposition'; catch, end
ax.Position = axPosPx;
ax.Units = 'normalized';
end

function localPreSaveFix(fig, ax, lgd, legendPosAxes, axesMarginPx)
drawnow;
localApplyAxesMarginPx(fig, ax, axesMarginPx);

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
        fr = getframe(fig);
        imwrite(fr.cdata, fn);

    case "exportgraphics"
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
try, vw.Quality = quality; catch, end

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
if isfield(ell,'confLevel') && ~isempty(ell.confLevel), confLevel = ell.confLevel; end
thr = chi2inv(confLevel, 3);
scale = sqrt(thr);

Sigma = (Sigma + Sigma')/2 + 1e-8*eye(3);
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
    'Marker','p', 'MarkerFaceColor', opt.ellColor, 'MarkerEdgeColor', opt.ellColor);

try
    set(get(get(hEll.surf,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');
    set(get(get(hEll.muMarker,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');
catch
end

X = Y3 - mu;
R = chol(Sigma);
Z = X / R;
md2 = sum(Z.^2, 2);
inside = (md2 <= thr);
end

function localPlaceLegendInsideAxes(ax, lgd, posAxesNorm)
ax.Units = 'normalized';
P = ax.Position;
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
        tok = regexp(hdr(ii), '(m\d{3,5})', 'tokens','once');
        if isempty(tok), mId{ii} = sprintf('mouse_%03d', ii);
        else,           mId{ii} = tok{1};
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
    tok = regexp(hdr(ii), '_(\d{6})', 'tokens','once');
    if isempty(tok)
        orderKey(ii) = ii;
    else
        d6 = char(tok{1}); % MMDDYY
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
