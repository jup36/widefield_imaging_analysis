function files = saveMotifsToGif(matFile, varargin)
% saveMotifsToGif  Save each spatiotemporal motif (k) in W_basis as a GIF.
%
% files = saveMotifsToGif(matFile, ...
%   'varName','W_basis', 'outDir','', 'colormap','parula', 'CLim',[], ...
%   'delay',0.10, 'loopCount',Inf, 'upsample',1, 'filePrefix','motif', ...
%   'disposalMethod','donotspecify')
%
% INPUT
%   matFile : path to .mat containing W_basis (P x K x T)
%
% OPTIONS
%   varName        : variable name inside MAT (default 'W_basis')
%   outDir         : output directory (default = matFile folder)
%   colormap       : colormap name/handle or Nx3 double (default 'parula')
%   CLim           : [lo hi] global per-motif range; [] => 2–98th pctile
%   delay          : scalar seconds between frames (default 0.10)
%   loopCount      : Inf by default
%   upsample       : integer (default 1)
%   filePrefix     : 'motif' by default
%   disposalMethod : 'donotspecify' (also accepts 'leaveinplace','restorebg',...)
%
% OUTPUT
%   files : 1xK cell array of saved GIF paths

% ---------- parse ----------
p = inputParser;
p.addRequired('matFile', @(s)ischar(s)||isstring(s));
p.addParameter('varName','W_basis', @(s)ischar(s)||isstring(s));
p.addParameter('outDir','', @(s)ischar(s)||isstring(s));
p.addParameter('colormap','parula', @(c)ischar(c)||isstring(c)||isa(c,'function_handle')|| ...
                                      (isnumeric(c)&&size(c,2)==3));
p.addParameter('CLim', [], @(x) isempty(x) || (isnumeric(x)&&numel(x)==2));
p.addParameter('delay', 0.10, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('loopCount', Inf, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('upsample', 1, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('filePrefix','motif', @(s)ischar(s)||isstring(s));
p.addParameter('disposalMethod','donotspecify', @(s)ischar(s)||isstring(s));
p.parse(matFile, varargin{:});
opt = p.Results;

% ---------- load ----------
S = load(matFile, char(opt.varName));
assert(isfield(S, char(opt.varName)), 'Variable "%s" not found in %s.', opt.varName, matFile);
W = S.(char(opt.varName));                 % P x K x T
assert(ndims(W)==3, 'W_basis must be P x K x T.');
[P,K,T] = size(W);

% infer spatial size (e.g., 4096 -> 64x64)
n = sqrt(P); assert(abs(n-round(n))<1e-9, 'P=%d must be a perfect square.', P);
n = round(n);

% ---------- output dir ----------
if isempty(opt.outDir)
    [folder,~,~] = fileparts(char(matFile)); outDir = folder;
else
    outDir = char(opt.outDir);
end
if exist(outDir,'dir')~=7, mkdir(outDir); end

% ---------- colormap (fixed 256×3 double in [0,1]) ----------
if isnumeric(opt.colormap)
    cmap = double(opt.colormap);
else
    if isa(opt.colormap,'function_handle'), cmap = feval(opt.colormap,256);
    else,                                   cmap = feval(char(opt.colormap),256);
    end
end
cmap = max(0,min(1,double(cmap)));
if size(cmap,1)~=256, cmap = imresize(cmap,[256 3],'nearest'); end

% ---------- per-motif export ----------
files = cell(1,K);
[~, base, ~] = fileparts(char(matFile));

for k = 1:K
    X = reshape(double(W(:,k,:)), n, n, T);

    % CLim per motif (robust)
    if isempty(opt.CLim)
        allv = X(:);
        lo = prctile(allv, 2);
        hi = prctile(allv, 98);
        if ~isfinite(lo) || ~isfinite(hi) || lo==hi
            lo = min(allv); hi = max(allv); if lo==hi, lo = hi-1; end
        end
    else
        lo = opt.CLim(1); hi = opt.CLim(2);
    end

    gifName = sprintf('%s_%s_%02d.gif', char(opt.filePrefix), base, k);
    gifPath = fullfile(outDir, gifName);
    if exist(gifPath,'file')==2, delete(gifPath); end   % important

    for t = 1:T
        A = (X(:,:,t) - lo) ./ max(hi - lo, eps);   % linear map to [0,1]
        A = max(0,min(1,A));
        % (optional) gamma boost for dim parts:
        A = A .^ 0.55;

        if opt.upsample>1
            A = imresize(A, opt.upsample, 'nearest');
        end

        Aind = gray2ind(A, size(cmap,1));          % uint8 indices

        if t==1
            imwrite(Aind, cmap, gifPath, 'gif', ...
                'LoopCount', opt.loopCount, ...
                'DelayTime', opt.delay, ...
                'DisposalMethod', char(opt.disposalMethod));
        else
            imwrite(Aind, cmap, gifPath, 'gif', ...
                'WriteMode','append', ...
                'DelayTime', opt.delay, ...
                'DisposalMethod', char(opt.disposalMethod));
        end
    end

    files{k} = gifPath;
    fprintf('Saved motif %d/%d -> %s\n', k, K, gifPath);
end
end
