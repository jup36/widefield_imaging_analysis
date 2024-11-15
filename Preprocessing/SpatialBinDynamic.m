function imds = SpatialBinDynamic(imstack, vmask, IgnoreNan)

% spatially bins image to achieve 64x64 output
% LucasPinto 2017, modified version by Junchol Park (11/1/24)
if nargin < 2
    vmask = [];
end
if nargin < 3
    IgnoreNan = 1;
end

% Get the original dimensions
[nX, nY, nZ] = size(imstack);

% Initialize the output image stack
imds = zeros(64, 64, nZ);

% Create bin edges for x and y, ensuring they partition into 64 equal sections
xbins = round(linspace(1, nX + 1, 65));
ybins = round(linspace(1, nY + 1, 65));

% Check if input is gpuArray and make output agree
if isa(imstack, 'gpuArray')
    imds = gpuArray(imds);
end

% Downsample and ignore vasculature pixels
for zz = 1:nZ
    im = imstack(:,:,zz);
    if ~isempty(vmask)
        im(vmask == 0) = nan;
    end
    for xx = 1:64
        for yy = 1:64
            vals = im(xbins(xx):xbins(xx+1)-1, ybins(yy):ybins(yy+1)-1);
            if IgnoreNan
                imds(xx, yy, zz) = nansum(vals(:));
            else
                imds(xx, yy, zz) = sum(vals(:));
            end
        end
    end
end
end
