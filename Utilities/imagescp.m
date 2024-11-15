function h = imagescp(img2d, varargin)
if nargin < 2
    climit = [-2 2]; 
elseif nargin == 2
    climit = varargin{1}; 
end
h = figure;
imagesc(img2d); 
colormap('parula')
clim(climit); 
pbaspect([1 1 1]); 
colorbar; 


end