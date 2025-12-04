function files = saveWbasisMotifsToGif(matFile, varargin)
% saveWbasisMotifsToGif  Load W_basis (P x K x T), reshape each to 64x64xT, save GIFs.
%
% files = saveWbasisMotifsToGif(matFile, ...
%   'varName','W_basis', 'outDir','', 'filePrefix','motif', ...
%   'colormap','magma', 'CLim',[], 'percentile',[2 98], ...
%   'gamma',0.55, 'delay',0.10, 'loopCount',Inf, 'upsample',1)
%
ip = inputParser;
ip.addRequired('matFile', @(s)ischar(s)||isstring(s));
ip.addParameter('varName','W_basis', @(s)ischar(s)||isstring(s));
ip.addParameter('outDir','', @(s)ischar(s)||isstring(s));
ip.addParameter('filePrefix','motif', @(s)ischar(s)||isstring(s));
ip.addParameter('colormap','magma', @(c) ischar(c)||isstring(c)||isa(c,'function_handle'));
ip.addParameter('CLim', [], @(x) isempty(x) || (isnumeric(x)&&numel(x)==2));
ip.addParameter('percentile', [2 98], @(x) isnumeric(x)&&numel(x)==2);
ip.addParameter('gamma', 0.55, @(x) isnumeric(x)&&isscalar(x)&&x>0);
ip.addParameter('delay', 0.10, @(x) isnumeric(x)&&isscalar(x)&&x>=0.02);
ip.addParameter('loopCount', Inf, @(x) isnumeric(x)&&isscalar(x)&&x>=0);
ip.addParameter('upsample', 1, @(x) isnumeric(x)&&isscalar(x)&&x>=1);
ip.parse(matFile, varargin{:});
opt = ip.Results;

S = load(matFile, char(opt.varName));
assert(isfield(S, char(opt.varName)), 'Variable "%s" not found in %s.', opt.varName, matFile);
W = S.(char(opt.varName));   % P x K x T
assert(ndims(W)==3, 'W_basis must be P x K x T.');
[P,K,T] = size(W);
n = sqrt(P); assert(abs(n-round(n))<1e-9, 'P=%d must be a perfect square.', P); n = round(n);

if isempty(opt.outDir)
    [folder,~,~] = fileparts(char(matFile)); outDir = folder;
else
    outDir = char(opt.outDir);
end
if exist(outDir,'dir')~=7, mkdir(outDir); end
[~, base, ~] = fileparts(char(matFile));

files = cell(1,K);
for k = 1:K
    Wk = reshape(W(:,k,:), n, n, T);
    gifPath = fullfile(outDir, sprintf('%s_%s_%02d.gif', char(opt.filePrefix), base, k));
    MotifsToGifs(Wk, gifPath, ...
        'colormap', opt.colormap, ...
        'CLim', opt.CLim, ...
        'percentile', opt.percentile, ...
        'gamma', opt.gamma, ...
        'delay', opt.delay, ...
        'loopCount', opt.loopCount, ...
        'upsample', opt.upsample, ...
        'bg', 'w');
    files{k} = gifPath;
end
end
