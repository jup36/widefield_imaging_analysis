function [u, v] = HS_flow(im1, im2, alpha, ite, uInitial, vInitial)
% Horn-Schunck optical flow method 
% Horn, B.K.P., and Schunck, B.G., Determining Optical Flow, AI(17), No.
% 1-3, August 1981, pp. 185-203 http://dspace.mit.edu/handle/1721.1/6337
%
% Usage:
% [u, v] = HS(im1, im2, alpha, ite, uInitial, vInitial, displayFlow)
% For an example, run this file from the menu Debug->Run or press (F5)
%
% -im1,im2 : two subsequent frames or images.
% -alpha : a parameter that reflects the influence of the smoothness term.
% -ite : number of iterations.
% -uInitial, vInitial : initial values for the flow. If available, the
% flow would converge faster and hence would need less iterations ; default is zero. 
% -displayFlow : 1 for display, 0 for no display ; default is 1.
% -displayImg : specify the image on which the flow would appear ( use an
% empty matrix "[]" for no image. )
%
% Author: Mohd Kharbat at Cranfield Defence and Security
% mkharbat(at)ieee(dot)org , http://mohd.kharbat.com
% Published under a Creative Commons Attribution-Non-Commercial-Share Alike
% 3.0 Unported Licence http://creativecommons.org/licenses/by-nc-sa/3.0/
%
% October 2008
% Rev: Jan 2009

%% Default parameters
if nargin<1 || nargin<2
    im1=imread('yos9.tif');
    im2=imread('yos10.tif');
end
if nargin<3
    alpha=1;
end
if nargin<4
    ite=100;
end
if nargin<5 || nargin<6
    uInitial = zeros(size(im1(:,:,1)));
    vInitial = zeros(size(im2(:,:,1)));
elseif size(uInitial,1) ==0 || size(vInitial,1)==0
    uInitial = zeros(size(im1(:,:,1)));
    vInitial = zeros(size(im2(:,:,1)));
end

 
%% Convert images to grayscale
if size(size(im1),2)==3
    im1=rgb2gray(im1);
end
if size(size(im2),2)==3
    im2=rgb2gray(im2);
end
im1=double(im1);
im2=double(im2);

im1=imgaussfilt(im1,[1 1]);
im2=imgaussfilt(im2,1.0);

%%
% Set initial value for the flow vectors
u = uInitial;
v = vInitial;

% Estimate spatiotemporal derivatives
[fx, fy, ft] = computeDerivatives(im1, im2);

% Averaging kernel
kernel_1=[1/12 1/6 1/12;1/6 0 1/6;1/12 1/6 1/12];
% kernel_2=0.25*[0 1 0;1 0 1;0 1 0];

%  hWaitBar = waitbar(0,sprintf('Processed HS iteration -'));
% Iterations
for i=1:ite
%     waitbar(i/ite,hWaitBar,sprintf('Processing HS iteration %d',i));
    % Compute local averages of the flow vectors
    uAvg=conv2(u,kernel_1,'same');
    vAvg=conv2(v,kernel_1,'same');
    % Compute flow vectors constrained by its local average and the optical flow constraints
    u= uAvg - ( fx .* ( ( fx .* uAvg ) + ( fy .* vAvg ) + ft ) ) ./ ( alpha^2 + fx.^2 + fy.^2); 
    v= vAvg - ( fy .* ( ( fx .* uAvg ) + ( fy .* vAvg ) + ft ) ) ./ ( alpha^2 + fx.^2 + fy.^2);
end
% close(hWaitBar) ;

u(isnan(u))=0;
v(isnan(v))=0;



end


function [fx, fy, ft] = computeDerivatives(im1, im2)

if size(im2,1)==0
    im2=zeros(size(im1));
end
% derivatives as in Barron
I = im1+im2;
fx= conv2(I/2,(1/12)*[-1 8 0 -8 1],'same');
fy= conv2(I/2,(1/12)*[-1 8 0 -8 1]','same');
ft = conv2(im1, 0.25*ones(2),'same') + conv2(im2, -0.25*ones(2),'same');
fx=-fx;fy=-fy;

end

