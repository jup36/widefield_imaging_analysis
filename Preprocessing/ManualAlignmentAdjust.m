function [opts] = ManualAlignmentAdjust(ref_img,opts)
ref_img = double(ref_img);

%show image
figure('name','Move rectangle and double-click to crop')
imagesc(ref_img); colormap gray
clim([min(ref_img(:)) max(min(ref_img(:)), max(ref_img(:)))])
%clim([min(ref_img(:)) max(min(ref_img(:)), max(ref_img(:))/2)]) % to adjust the overall brightness of the image
%imcontrast(gca) % interactive adjustment

%Set to correct aspect ratio and shift to ~center of 1080p screen
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1)-200, pos(2)-200, size(ref_img,2), size(ref_img,1)]);
hold on

% opts.crop_w = 448; % original: 540
% opts.crop_h = 448; % original: 540

%cropping rectangle
rect = imrect(gca,[0 0 opts.crop_w opts.crop_h]);
fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
setPositionConstraintFcn(rect,fcn); 
setFixedAspectRatioMode(rect,1);
setResizable(rect,0);
crop_position = wait(rect);
close

% Crop the ref image
ref_img = imcrop(ref_img,crop_position);

%Get angle of brain 
angle = setMidlineAngle(ref_img);

%Roate the brain to vertical 
if angle <= 90 %if angled to right of yaxis
    aligned_img = imrotate(ref_img,-angle);
else %if angled left of yaxis
    aligned_img = imrotate(ref_img,180-angle);
end

%pad to keep cropping within bounds across lots of 
%different image sizes after rotating
aligned_img = padarray(aligned_img,[60,60]); 
bregma = setBregma(aligned_img);

%get coordinates for conservative cropping. 
%This make all recordings identically positioned with respect to bregma
% opts.x_bregma_margin = 280-(540-448)/2; % original: 280 (adjusted for imaging data from PNI284)
% opts.y_bregma_margin = 230-(540-448)/2; % original: 230 (adjusted for imaging data from PNI284)

x_crop_cord = (bregma(1)-opts.x_bregma_margin);
y_crop_cord = (bregma(2)-opts.y_bregma_margin); 

%crop the alligned image
cropped_aligned_img = double(imcrop(aligned_img,[x_crop_cord,y_crop_cord,...
    opts.crop_w,opts.crop_h]));

%crop can lead to inconsistent # of pixels. Make equal. 
cropped_aligned_img = cropped_aligned_img([1:opts.crop_h],[1:opts.crop_w]); 

%shift bregma to new coords in cropped image
bregma = [bregma(1)-x_crop_cord,bregma(2)-y_crop_cord];
%figure; imagesc(cropped_aligned_img); hold on; plot(bregma(1), bregma(2), 'r*', 'MarkerSize', 10, 'LineWidth', 1.5); % to confirm the location of bregma on 'cropped_aligned_img'

%store variables in preprocessing_log structure
opts.bregma = bregma; 
opts.crop_cord = [x_crop_cord,y_crop_cord];
opts.aligned_img = aligned_img;
opts.cropped_aligned_img = cropped_aligned_img;
opts.angle = angle; 
opts.crop_position =crop_position; 


end





