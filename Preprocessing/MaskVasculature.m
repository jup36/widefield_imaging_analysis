function options = MaskVasculature(ref_img_aligned,options)
% ref_img_aligned: a reference image cropped and aligned (prepro_log.cropped_alligned_img). 
% options: prepro_log that contains the result from ManualAlignmentAdjust. 

if options.mask_vasc %mask vasc
    %Indentify vasculature by subtracting a large median filter  (e.g. smoothed image
    %form the base image. This excentuates dark areas that are surrounded by
    %light areas (e.g. vasculature). 
    smoothed = medfilt2(ref_img_aligned,[125 125]); %Use large neighboorhood e.g. 1/4 of size
    mask = ref_img_aligned-smoothed; 

    % Calculate the threshold to remove pixels lower than 'x' standard deviations below the mean
    meanMask = nanmean(mask(:)); % Calculate the mean over the entire mask
    stdMask = nanstd(mask(:)); % Calculate the standard deviation over the entire mask

    % Apply threshold: Keep pixels 'x' standard deviations below the mean
    x = options.vasc_std; % 'x' standard deviations specified by user
    threshold = meanMask - (x * stdMask); % Calculate the threshold value
    mask = mask >= threshold; % Apply threshold
    
    %threshold to remove pxs lower that 'x' standard deviations below mean
    %of Vascmask (e.g. the vasculature and other dark blemishes). 
    %mask = mask>=(nanmean(nanmean(mask))-(options.vasc_std*nanstd(nanstd(mask))));

    %close the image to clean up noise
    se = strel('disk',options.close_disk_size);
    mask = imclose(mask,se); 

    %Add a premade mask to outline the brain.
    if options.mask_brain_outline 
        temp = load(options.mask_brain_outline_dir); 
        %Crop the mask to the size of the image. Assumes top left alligned
        temp.brainoutline = temp.brainoutline(1:options.crop_h, 1:options.crop_w);
        mask = mask+temp.brainoutline;
        mask = mask==2;               
    end
    
    if options.verbose || options.manual_mask
        ref_masked = ref_img_aligned.*mask; 
        imshowpair(ref_img_aligned,ref_masked);
        clim([min(ref_img_aligned(:)), max(ref_img_aligned(:))/4])
        %cmin = min(min(ref_img_aligned(:)), min(ref_masked(:))); % Minimum pixel value across both images
        %cmax = max(max(ref_img_aligned(:)), max(ref_masked(:))); % Maximum pixel value across both images
        %clim([cmin, cmax])
        
        hold on; 
        scatter(options.bregma(1),options.bregma(2),'*');
        title('Cropped Img with Vascmask overlay');
    end

    %Manually Nan regions and unmask others 
    if options.manual_mask
        dlg = questdlg('Would you like to manually mask regions?', ...
        'mask ROIS','Yes','No','No');
        while 1
            switch dlg
                case 'Yes' %mask additional region
                    close all                
                    figure('name','Draw ROI around additional regions to NaN');
                    imshowpair(ref_img_aligned,ref_img_aligned.*mask);
                    %clim([min(ref_img_aligned(:)) max(min(ref_img_aligned(:)), max(ref_img_aligned(:))/4)])
                    manual_mask = impoly;
                    wait(manual_mask);               
                    manual_mask = manual_mask.createMask();
                    mask(manual_mask)=0; close; 
                    figure('name','Current masked Image')
                    imshowpair(ref_img_aligned,ref_img_aligned.*mask); 
                    %clim([min(ref_img_aligned(:)) max(min(ref_img_aligned(:)), max(ref_img_aligned(:))/4)])
               case 'No'; break %end dlg
            end
            choice = questdlg(sprintf('Would you like to mask additional regions'),...
                'mask ADDITIONAL ROIS',...
                'Yes','No','No');
            switch choice; case 'Yes';case 'No'; break; end
        end
        dlg = questdlg('Would you like to manually UNmask regions?', ...
        'mask ROIS','Yes','No','No');
        while 1 %Manuall UNmask regions 
            switch dlg
                case 'Yes' %mask additional region
                    close all
                    figure('name','Draw ROI around additional regions to NaN');
                    imshowpair(ref_img_aligned,ref_img_aligned.*mask);
                    %clim([min(ref_img_aligned(:)) max(min(ref_img_aligned(:)), max(ref_img_aligned(:))/4)])
                    manual_mask = impoly;
                    wait(manual_mask);               
                    manual_mask = manual_mask.createMask(); %his is 1 in the desired ROI          
                    mask(manual_mask)=1; close
                    figure('name','Current masked Image')
                    imshowpair(ref_img_aligned,ref_img_aligned.*mask);
                    %clim([min(ref_img_aligned(:)) max(min(ref_img_aligned(:)), max(ref_img_aligned(:))/4)])
               case 'No'; break
            end
            choice = questdlg(sprintf('Would you like to unmask additional regions'),...
                'mask ADDITIONAL ROIS',...
                'Yes','No','No');
            switch choice; case 'Yes';case 'No'; break; end
        end
    end

    if options.verbose || options.manual_mask
        close all
        imshowpair(ref_img_aligned,ref_img_aligned.*mask);
        %clim([min(ref_img_aligned(:)) max(min(ref_img_aligned(:)), max(ref_img_aligned(:))/4)])
        hold on; 
        scatter(options.bregma(1),options.bregma(2),'*');
        title('FINAL cropped Img with Vascmask overlay');
            choice = questdlg(sprintf('Would you like to flag this brain for imperfections?'),...
            'Flag brain',...
            'Yes','No','No');
       switch choice
           case 'Yes'
               options.flagged_brain = 1;
           case 'No' 
               options.flagged_brain = 0; 
       end
    elseif options.maskVasculature==0 && options.maskBrainOutline == 1 %only apply exterior mask
        mask = brainoutline;
        if options.verbose
            close all
            imshowpair(ref_img_aligned,ref_img_aligned.*mask);
            %clim([min(ref_img_aligned(:)) max(min(ref_img_aligned(:)), max(ref_img_aligned(:))/4)])
            hold on; 
            scatter(options.bregma(1),options.bregma(2),'*');
            title('FINAL cropped Img with masked Brain Overlay');
        end        
    else %apply no mask
        mask = ones(size(ref_img_aligned));
    end
    options.mask = mask; 
    options.masked_ref_img_aligned = ref_img_aligned.*mask; 

end %function end
    