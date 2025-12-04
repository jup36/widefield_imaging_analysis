function MotifsToGifs(W, save_path, varargin)
% MotifToGif  Save a 3-D movie (H x W x T) to an animated GIF using frame capture.
%
% MotifToGif(W, save_path, ...
%   'colormap','parula', ...        % colormap name/handle
%   'CLim',[], ...                  % [lo hi]; [] => percentile
%   'percentile', [2 98], ...       % used only if CLim=[]
%   'gamma', 0.55, ...              % <1 brightens dim parts
%   'delay', 0.10, ...              % sec per frame
%   'loopCount', Inf, ...           % loops (Inf ok here)
%   'upsample', 1, ...              % integer upsample for display
%   'bg', 'w', ...                  % figure background color
%   'type','simple')                % kept for backward compat
%
% W: H x W x T numeric (movie)

% ---- parse ----
ip = inputParser;
ip.addRequired('W', @(x) isnumeric(x) && ndims(x)==3);
ip.addRequired('save_path', @(s) ischar(s)||isstring(s));
ip.addParameter('colormap','parula', @(c) ischar(c)||isstring(c)||isa(c,'function_handle'));
ip.addParameter('CLim', [], @(x) isempty(x) || (isnumeric(x)&&numel(x)==2));
ip.addParameter('percentile', [2 98], @(x) isnumeric(x)&&numel(x)==2);
ip.addParameter('gamma', 0.55, @(x) isnumeric(x)&&isscalar(x)&&x>0);
ip.addParameter('delay', 0.10, @(x) isnumeric(x)&&isscalar(x)&&x>=0.02);
ip.addParameter('loopCount', Inf, @(x) isnumeric(x)&&isscalar(x)&&x>=0);
ip.addParameter('upsample', 1, @(x) isnumeric(x)&&isscalar(x)&&x>=1);
ip.addParameter('bg', 'w', @(x) ischar(x)||isstring(x)|| (isnumeric(x)&&numel(x)==3));
ip.addParameter('type','simple', @(s) ischar(s)||isstring(s));
ip.parse(W, save_path, varargin{:});
opt = ip.Results;

[H,Wid,T] = size(W);
W = double(W);

% ---- decide CLim ----
if isempty(opt.CLim)
    pr = prctile(W(:), sort(opt.percentile));
    lo = pr(1); hi = pr(2);
    if ~isfinite(lo) || ~isfinite(hi) || lo==hi
        lo = min(W(:)); hi = max(W(:));
        if lo==hi, lo = hi-1; end
    end
else
    lo = opt.CLim(1); hi = opt.CLim(2);
end

% ---- figure for rendering (invisible) ----
% Use a fixed pixel size to keep frames identical.
scale = max(1, round(300 / min(H,Wid)));          % heuristic
M = H*opt.upsample*scale; N = Wid*opt.upsample*scale;

hFig = figure('Visible','off','Color',opt.bg,'Units','pixels','Position',[100 100 N M]);
ax = axes('Parent',hFig,'Position',[0 0 1 1]);    % full-bleed axes
axis(ax,'tight','ij'); axis(ax,'off');
colormap(ax, opt.colormap);
set(ax,'CLim',[lo hi]);

% ---- render frames -> RGB -> indexed -> write ----
if exist(save_path,'file'), delete(save_path); end
for t = 1:T
    A = squeeze(W(:,:,t));
    % gamma mapping in data space for better dim visibility
    A = (A - lo) ./ max(hi - lo, eps);            % [0,1]
    A = max(0,min(1,A)).^opt.gamma;

    if opt.upsample>1
        A = imresize(A, opt.upsample, 'nearest');
    end

    imagesc(ax, A, [0 1]);                        % display uses normalized range
    axis(ax,'image','off');
    drawnow;                                       %#ok<*DRAWNOW>

    frameRGB = frame2im(getframe(hFig));          % truecolor RGB
    [Aind,map] = rgb2ind(frameRGB, 256);          % per-frame palette (robust)

    if t==1
        imwrite(Aind, map, save_path, 'gif', 'LoopCount', opt.loopCount, 'DelayTime', opt.delay);
    else
        imwrite(Aind, map, save_path, 'gif', 'WriteMode','append', 'DelayTime', opt.delay);
    end
end

close(hFig);
end
