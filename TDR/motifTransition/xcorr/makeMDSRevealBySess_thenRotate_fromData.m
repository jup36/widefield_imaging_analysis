function out = makeMDSRevealBySess_thenRotate_fromData( ...
    Y_full, sessInfo, cmap, ellipsoid_full, outDir, baseName, varargin)
% makeMDSRevealBySess_thenRotate_fromData
% Reveals sessions step-by-step across animals (synced by session index),
% saves frames to disk, then compiles to MP4.
%
% IMPORTANT: plotMDSWithEllipsoid may create/choose its own figure.
% We therefore capture the figure/axes actually used, and always save that.

% -------------------- parse --------------------
p = inputParser;
p.addRequired('Y_full', @(x)isnumeric(x) && ndims(x)==2);
p.addRequired('sessInfo');
p.addRequired('cmap');
p.addRequired('ellipsoid_full');
p.addRequired('outDir', @(x)ischar(x)||isstring(x));
p.addRequired('baseName', @(x)ischar(x)||isstring(x));

% plot params -> plotMDSWithEllipsoid
p.addParameter('dims', [1 2 3], @(x)isnumeric(x)&&numel(x)==3);
p.addParameter('ellAlpha', 0.10, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('ellN', 50, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('minAlpha', 0.25, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('maxAlpha', 1.0, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('minDotSize', 24, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('dotScaleMax', 3.0, @(x)isnumeric(x)&&isscalar(x));

% frame/video params
p.addParameter('fps', 30, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('imgExt', 'png', @(s)ischar(s)||isstring(s)); % png/jpg/tif
p.addParameter('renderer', 'opengl', @(s)ischar(s)||isstring(s));
p.addParameter('figPos', [100 100 1200 900], @(x)isnumeric(x)&&numel(x)==4);
p.addParameter('figWidthScale', 1.5, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('lockVis3d', true, @(x)islogical(x)&&isscalar(x));

% views
p.addParameter('defaultView', [-111 28], @(x)isnumeric(x)&&numel(x)==2);

% rotations
p.addParameter('horizAzRange', [0 90], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('horizSeconds', 4, @(x)isnumeric(x)&&isscalar(x)&&x>0);

p.addParameter('tiltElRange', [-5 45], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('tiltSeconds', 3, @(x)isnumeric(x)&&isscalar(x)&&x>0);

% reveal
p.addParameter('revealSeconds', 6, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('revealHoldSeconds', 0.0, @(x)isnumeric(x)&&isscalar(x)&&x>=0);

% date / animal parsing
p.addParameter('dateField', '', @(s)ischar(s)||isstring(s));
p.addParameter('dateFormat', '', @(s)ischar(s)||isstring(s));
p.addParameter('animalField', '', @(s)ischar(s)||isstring(s));

p.addParameter('verbose', true, @(x)islogical(x)&&isscalar(x));

p.parse(Y_full, sessInfo, cmap, ellipsoid_full, outDir, baseName, varargin{:});
opt = p.Results;

N = size(Y_full,1);
nInfo = localSessInfoLen(sessInfo);
if nInfo ~= N
    error('sessInfo length (%d) must match size(Y_full,1) (%d).', nInfo, N);
end

% -------------------- output dirs --------------------
outDir = char(opt.outDir);
if ~exist(outDir,'dir'), mkdir(outDir); end

framesDir = fullfile(outDir, [char(opt.baseName) '_frames']);
if ~exist(framesDir,'dir'), mkdir(framesDir); end

imgExt = lower(string(opt.imgExt));
if imgExt == "jpeg", imgExt = "jpg"; end

% -------------------- extract dates + animals --------------------
dt     = localExtractDates(sessInfo, opt.dateField, opt.dateFormat);
animal = localExtractAnimalIDs(sessInfo, opt.animalField);

[uAnimals, ~, g] = unique(animal, 'stable');
nAnimals = numel(uAnimals);

idxByAnimal = cell(nAnimals,1);
for a = 1:nAnimals
    idx = find(g == a);
    dta = dt(idx);

    isBad = isnat(dta);
    key = dta;
    key(isBad) = datetime(3000,1,1);

    [~, ord] = sort(key, 'ascend');
    idxByAnimal{a} = idx(ord);
end
maxSess = max(cellfun(@numel, idxByAnimal));

if opt.verbose
    fprintf('[Reveal-all] animals=%d, totalSessions=%d, maxSessAcrossAnimals=%d\n', ...
        nAnimals, N, maxSess);
end

% -------------------- state for saving --------------------
frameFiles = strings(0,1);
frameIdx = 0;

% We'll hold onto the *actual* fig/ax used by plotMDSWithEllipsoid
figUsed = [];
axUsed  = [];

    function [figNow, axNow] = localPlotAndCapture(idxKeep)
        % Plot and capture the fig/ax that actually received the graphics.
        % We do this by snapshotting existing figures before the plot call,
        % then checking what exists after, and selecting the newest/active one.

        figsBefore = findall(0,'Type','figure');

        % Plot (this may create a new figure or reuse an existing)
        plotMDSWithEllipsoid(Y_full(idxKeep,:), localSubsetSessInfo(sessInfo, idxKeep), ...
            cmap, ellipsoid_full, ...
            'view', opt.defaultView, ...
            'ellAlpha', opt.ellAlpha, 'ellN', opt.ellN, ...
            'minAlpha', opt.minAlpha, 'maxAlpha', opt.maxAlpha, ...
            'minDotSize', opt.minDotSize, 'dotScaleMax', opt.dotScaleMax, ...
            'dims', opt.dims);

        drawnow;

        figsAfter = findall(0,'Type','figure');

        % Prefer: a new figure created by the call
        newFigs = setdiff(figsAfter, figsBefore);
        if ~isempty(newFigs)
            figNow = newFigs(1);
        else
            % Else: use current figure
            figNow = gcf;
        end

        % Configure it once
        set(figNow,'Renderer', char(opt.renderer));
        set(figNow,'Units','pixels');
        figPos = opt.figPos;
        figPos(3) = round(figPos(3) * opt.figWidthScale);
        set(figNow,'Position', figPos);
        set(figNow,'Visible','on');
        set(figNow,'InvertHardcopy','off');

        axNow = get(figNow,'CurrentAxes');
        if isempty(axNow) || ~isgraphics(axNow,'axes')
            axList = findall(figNow,'Type','axes');
            if ~isempty(axList), axNow = axList(1); end
        end
        if isempty(axNow) || ~isgraphics(axNow,'axes')
            error('Could not find axes after plotMDSWithEllipsoid call. The plotting function may not be creating axes.');
        end

        if opt.lockVis3d
            try, axis(axNow,'vis3d'); catch, end
        end
        view(axNow, opt.defaultView);
        drawnow;
    end

    function saveCurrentFrame()
        % Save the axes that actually contains the plot
        if isempty(figUsed) || ~isgraphics(figUsed,'figure') || isempty(axUsed) || ~isgraphics(axUsed,'axes')
            error('Internal fig/ax not set. Plot must run before saving.');
        end

        frameIdx = frameIdx + 1;
        fname = sprintf('%06d_%s.%s', frameIdx, char(opt.baseName), char(imgExt));
        fpath = fullfile(framesDir, fname);

        drawnow;
        pause(0.01); %#ok<PAUS>

        try
            exportgraphics(axUsed, fpath, 'Resolution', 200, 'BackgroundColor','white');
        catch
            % fallback: print the *figure* (not axes) to file
            try
                print(figUsed, fpath, ['-d' char(imgExt)], '-r200');
            catch
                fr = getframe(figUsed);
                imwrite(fr.cdata, fpath);
            end
        end

        frameFiles(end+1,1) = string(fpath);
    end

% ===================== PHASE A: REVEAL =====================
pausePerFrame = opt.revealSeconds / max(maxSess,1);

for t = 1:maxSess
    idxKeep = [];
    for a = 1:nAnimals
        nA = numel(idxByAnimal{a});
        kk = min(t, nA);
        if kk > 0
            idxKeep = [idxKeep; idxByAnimal{a}(1:kk)]; %#ok<AGROW>
        end
    end
    idxKeep = sort(idxKeep,'ascend');

    % Plot + capture the correct fig/ax
    [figUsed, axUsed] = localPlotAndCapture(idxKeep);

    % Save frame
    saveCurrentFrame();

    if pausePerFrame > 0
        pause(pausePerFrame); %#ok<PAUS>
    end
end

if opt.revealHoldSeconds > 0
    nHold = max(1, round(opt.revealHoldSeconds * opt.fps));
    for i = 1:nHold
        saveCurrentFrame();
    end
end

% ===================== PHASE B: HORIZONTAL ROTATION =====================
% Plot full once and reuse the same fig/ax
[figUsed, axUsed] = localPlotAndCapture( (1:N).' );

hFrames = max(2, round(opt.horizSeconds * opt.fps));
az = linspace(opt.horizAzRange(1), opt.horizAzRange(2), hFrames);
elFixed = opt.defaultView(2);

for i = 1:numel(az)
    view(axUsed, [az(i) elFixed]);
    drawnow;
    saveCurrentFrame();
end

% ===================== PHASE C: TILT =====================
tFrames = max(2, round(opt.tiltSeconds * opt.fps));
%el = linspace(opt.tiltElRange(1), opt.tiltElRange(2), tFrames);
el = linspace(opt.defaultView(2), opt.tiltElRange(2), tFrames);

azFixed = opt.horizAzRange(2);

for i = 1:numel(el)
    view(axUsed, [azFixed el(i)]);
    drawnow;
    saveCurrentFrame();
end

% ===================== PHASE D: COMPILE VIDEO =====================
mp4File = fullfile(outDir, [char(opt.baseName) '.mp4']);
vw = VideoWriter(mp4File, 'MPEG-4');
vw.FrameRate = opt.fps;
open(vw);

im0 = imread(frameFiles(1));
[Ht, Wt, ~] = size(im0);
Ht = Ht - mod(Ht,2);
Wt = Wt - mod(Wt,2);

for i = 1:numel(frameFiles)
    im = imread(frameFiles(i));
    im = localToFixedEven(im, Ht, Wt);
    writeVideo(vw, im);
end
close(vw);

out.mp4File    = mp4File;
out.framesDir  = framesDir;
out.frameFiles = frameFiles;

if opt.verbose
    fprintf('Saved frames: %s\n', framesDir);
    fprintf('Saved mp4:    %s\n', mp4File);
end

end % main


% ------------------------- helpers -------------------------

function n = localSessInfoLen(sessInfo)
if istable(sessInfo)
    n = height(sessInfo);
elseif isstruct(sessInfo)
    n = numel(sessInfo);
elseif iscell(sessInfo)
    n = numel(sessInfo);
else
    n = numel(sessInfo);
end
end

function sub = localSubsetSessInfo(sessInfo, idx)
if istable(sessInfo)
    sub = sessInfo(idx,:);
elseif isstruct(sessInfo)
    sub = sessInfo(idx);
elseif iscell(sessInfo)
    sub = sessInfo(idx);
else
    sub = sessInfo(idx);
end
end

function dt = localExtractDates(sessInfo, dateFieldOverride, dateFormatOverride)
n = localSessInfoLen(sessInfo);
dt = repmat(datetime(NaT), n, 1);

cands = ["date","Date","sessionDate","session_date","sessDate","sess_date", ...
         "day","Day","dt","datetime","time","Time"];

if strlength(string(dateFieldOverride)) > 0
    cands = [string(dateFieldOverride) cands];
end

for i = 1:n
    v = [];

    if istable(sessInfo)
        for c = cands
            if any(strcmp(sessInfo.Properties.VariableNames, c))
                v = sessInfo{i, c};
                break;
            end
        end
    elseif isstruct(sessInfo)
        for c = cands
            if isfield(sessInfo, c)
                v = sessInfo(i).(c);
                break;
            end
        end
    end

    if isempty(v), dt(i) = datetime(NaT); continue; end

    try
        if isdatetime(v)
            dt(i) = v;
        elseif isnumeric(v)
            if v > 7e5
                dt(i) = datetime(v, 'ConvertFrom','datenum');
            else
                dt(i) = datetime(NaT);
            end
        elseif isstring(v) || ischar(v)
            s = string(v);
            if strlength(string(dateFormatOverride)) > 0
                dt(i) = datetime(s, 'InputFormat', char(dateFormatOverride));
            else
                dt(i) = datetime(s);
            end
        elseif iscell(v)
            vv = v{1};
            if isdatetime(vv)
                dt(i) = vv;
            elseif isstring(vv) || ischar(vv)
                s = string(vv);
                if strlength(string(dateFormatOverride)) > 0
                    dt(i) = datetime(s, 'InputFormat', char(dateFormatOverride));
                else
                    dt(i) = datetime(s);
                end
            else
                dt(i) = datetime(NaT);
            end
        else
            dt(i) = datetime(NaT);
        end
    catch
        dt(i) = datetime(NaT);
    end
end
end

function animal = localExtractAnimalIDs(sessInfo, animalFieldOverride)
n = localSessInfoLen(sessInfo);
animal = strings(n,1);

cands = ["mId","mid","mouse","Mouse","mouseId","mouseID","animal","Animal","animalId","animalID","subject","Subject"];
if strlength(string(animalFieldOverride)) > 0
    cands = [string(animalFieldOverride) cands];
end

for i = 1:n
    v = [];

    if istable(sessInfo)
        for c = cands
            if any(strcmp(sessInfo.Properties.VariableNames, c))
                v = sessInfo{i, c};
                break;
            end
        end
    elseif isstruct(sessInfo)
        for c = cands
            if isfield(sessInfo, c)
                v = sessInfo(i).(c);
                break;
            end
        end
    end

    if isempty(v)
        animal(i) = "animal_" + string(i);
        continue;
    end

    if iscell(v), v = v{1}; end
    animal(i) = string(v);
end
end

function im = localToFixedEven(im, Ht, Wt)
[H,W,C] = size(im);
if C == 1
    im = repmat(im, [1 1 3]);
    [H,W,~] = size(im);
end

if H > Ht
    y0 = floor((H-Ht)/2)+1;
    im = im(y0:y0+Ht-1,:,:);
elseif H < Ht
    padTop = floor((Ht-H)/2);
    padBot = (Ht-H) - padTop;
    im = padarray(im, [padTop 0], 0, 'pre');
    im = padarray(im, [padBot 0], 0, 'post');
end

if W > Wt
    x0 = floor((W-Wt)/2)+1;
    im = im(:, x0:x0+Wt-1,:);
elseif W < Wt
    padLeft  = floor((Wt-W)/2);
    padRight = (Wt-W) - padLeft;
    im = padarray(im, [0 padLeft], 0, 'pre');
    im = padarray(im, [0 padRight], 0, 'post');
end

im = im(1:Ht, 1:Wt, 1:3);
end
