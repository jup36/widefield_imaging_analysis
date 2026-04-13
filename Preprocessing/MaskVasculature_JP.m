function options = MaskVasculature_JP(ref_img_aligned, options)
% ref_img_aligned: a reference image cropped and aligned (prepro_log.cropped_alligned_img).
% options: prepro_log that contains the result from ManualAlignmentAdjust.
%
% New behavior:
%   - User can iteratively adjust options.vasc_std before manual mask/unmask.
%   - GUI shows current preview and asks whether current std is satisfactory.
%   - If "No", prompts for a new numeric std value and updates preview.
%   - Repeats until user chooses "Yes".

if options.mask_vasc % mask vasc

    % Build initial auto mask, allowing user to iteratively tune vasc_std
    [mask, options] = chooseVascStdAndBuildMask(ref_img_aligned, options);

    % Manually Nan regions and unmask others
    if options.manual_mask
        dlg = questdlg('Would you like to manually mask regions?', ...
            'mask ROIS','Yes','No','No');
        while 1
            switch dlg
                case 'Yes' % mask additional region
                    close all
                    hfig1 = figure('Name','Draw ROI around additional regions to mask');
                    [~, hBaseImage1, ~] = showImgAndMask(ref_img_aligned, ref_img_aligned.*mask, hfig1);
                    hold on;
                    scatter(options.bregma(1),options.bregma(2),'*');
                    title(sprintf('Manual MASK step (vasc\\_std = %.3g)', options.vasc_std));

                    manual_mask = impoly(gca);
                    wait(manual_mask);
                    manual_mask = manual_mask.createMask(hBaseImage1);
                    mask(manual_mask) = 0;
                    close;

                    hfig2 = figure('Name','Current masked Image');
                    [~,~,~] = showImgAndMask(ref_img_aligned, ref_img_aligned.*mask, hfig2);
                    hold on;
                    scatter(options.bregma(1),options.bregma(2),'*');
                    title('Current masked image');
                case 'No'
                    break
            end

            choice = questdlg('Would you like to mask additional regions?', ...
                'mask ADDITIONAL ROIS', ...
                'Yes','No','No');
            switch choice
                case 'Yes'
                    dlg = 'Yes';
                case 'No'
                    break
            end
        end

        dlg = questdlg('Would you like to manually UNmask regions?', ...
            'mask ROIS','Yes','No','No');
        while 1
            switch dlg
                case 'Yes' % unmask additional region
                    close all
                    hfig3 = figure('Name','Draw ROI around regions to UNmask');
                    [~, hBaseImage3, ~] = showImgAndMask(ref_img_aligned, ref_img_aligned.*mask, hfig3);
                    hold on;
                    scatter(options.bregma(1),options.bregma(2),'*');
                    title(sprintf('Manual UNMASK step (vasc\\_std = %.3g)', options.vasc_std));

                    manual_mask = impoly(gca);
                    wait(manual_mask);
                    manual_mask = manual_mask.createMask(hBaseImage3); % 1 in desired ROI
                    mask(manual_mask) = 1;
                    close

                    hfig4 = figure('Name','Current masked Image');
                    [~,~,~] = showImgAndMask(ref_img_aligned, ref_img_aligned.*mask, hfig4);
                    hold on;
                    scatter(options.bregma(1),options.bregma(2),'*');
                    title('Current masked image');
                case 'No'
                    break
            end

            choice = questdlg('Would you like to unmask additional regions?', ...
                'UNmask ADDITIONAL ROIS', ...
                'Yes','No','No');
            switch choice
                case 'Yes'
                    dlg = 'Yes';
                case 'No'
                    break
            end
        end
    end

    if options.verbose || options.manual_mask
        close all
        showImgAndMask(ref_img_aligned, ref_img_aligned.*mask);
        hold on;
        scatter(options.bregma(1),options.bregma(2),'*');
        title(sprintf('FINAL cropped Img with Vascmask overlay (vasc\\_std = %.3g)', options.vasc_std));

        choice = questdlg('Would you like to flag this brain for imperfections?', ...
            'Flag brain', ...
            'Yes','No','No');
        switch choice
            case 'Yes'
                options.flagged_brain = 1;
            case 'No'
                options.flagged_brain = 0;
        end

    elseif options.maskVasculature==0 && options.maskBrainOutline==1 % only apply exterior mask
        mask = brainoutline;
        if options.verbose
            close all
            showImgAndMask(ref_img_aligned, ref_img_aligned.*mask);
            hold on;
            scatter(options.bregma(1),options.bregma(2),'*');
            title('FINAL cropped Img with masked Brain Overlay');
        end
    else % apply no mask
        mask = ones(size(ref_img_aligned));
    end

    options.mask = mask;
    options.masked_ref_img_aligned = ref_img_aligned.*mask;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [mask, options] = chooseVascStdAndBuildMask(ref_img_aligned, options)
        % Iteratively preview auto vasculature mask and let user adjust options.vasc_std

        keepGoing = true;

        while keepGoing
            mask = computeAutoMask(ref_img_aligned, options);

            % Show preview
            hfig = figure('Name','Vasculature mask preview');
            showImgAndMask(ref_img_aligned, ref_img_aligned.*mask, hfig);
            hold on;
            scatter(options.bregma(1),options.bregma(2),'*');
            title(sprintf('Auto vasculature mask preview (current vasc_std = %.3g)', options.vasc_std));

            % Ask if satisfactory
            choice = questdlg( ...
                sprintf('Current vasc_std = %.3g.\nIs this masking satisfactory?', options.vasc_std), ...
                'Confirm vasc\_std', ...
                'Yes','No','Yes');

            switch choice
                case 'Yes'
                    keepGoing = false;
                    close(hfig);

                case 'No'
                    close(hfig);

                    validInput = false;
                    while ~validInput
                        answer = inputdlg( ...
                            sprintf('Enter a new numeric value for vasc_std.\nCurrent value: %.3g', options.vasc_std), ...
                            'Update vasc\_std', ...
                            [1 50], ...
                            {num2str(options.vasc_std)});

                        % User cancelled -> keep current value and ask again using same preview loop
                        if isempty(answer)
                            validInput = true;
                            % do nothing, keep current options.vasc_std
                            continue;
                        end

                        newVal = str2double(answer{1});
                        if isfinite(newVal) && isscalar(newVal) && newVal > 0
                            options.vasc_std = newVal;
                            validInput = true;
                        else
                            uiwait(warndlg('Please enter a positive numeric scalar for vasc_std.', ...
                                'Invalid input', 'modal'));
                        end
                    end

                otherwise
                    % User closed dialog; treat like "Yes" to avoid trapping them
                    keepGoing = false;
                    close(hfig);
            end
        end
    end

    function mask = computeAutoMask(ref_img_aligned, options)
        % Build automatic vasculature mask from current options.vasc_std

        % Identify vasculature by subtracting a large median filter
        smoothed = medfilt2(ref_img_aligned, [125 125]); % large neighborhood
        maskRaw = ref_img_aligned - smoothed;

        % Threshold using options.vasc_std SD below mean
        meanMask = nanmean(maskRaw(:));
        stdMask = nanstd(maskRaw(:));
        x = options.vasc_std;
        threshold = meanMask - (x * stdMask);
        mask = maskRaw >= threshold;

        % Close image to clean up noise
        se = strel('disk', options.close_disk_size);
        mask = imclose(mask, se);

        % Add premade brain outline mask if requested
        if options.mask_brain_outline
            temp = load(options.mask_brain_outline_dir);
            if ~isequal(size(temp.brainoutline), size(ref_img_aligned))
                temp.brainoutline = rescaleBrainOutline(temp.brainoutline, size(ref_img_aligned));
            end

            % Crop to image size
            temp.brainoutline = temp.brainoutline(1:options.crop_h, 1:options.crop_w);
            mask = mask + temp.brainoutline;
            mask = mask == 2;
        end
    end

    function [hfig, hBaseImg, hOverlay] = showImgAndMask(img, maskedImg, varargin)

        if nargin < 3 || ~isa(varargin{1}, 'matlab.ui.Figure')
            hfig = figure;
        else
            hfig = varargin{1};
        end

        % Scale main image for visibility
        minImg = min(img(:));
        maxImg = max(img(:));
        scaledImage = (img - minImg) / (maxImg - minImg + eps);

        hAxes = axes('Parent', hfig);
        hBaseImg = imshow(scaledImage, 'Parent', hAxes);
        hold on;

        % Convert maskedImg to same display scale
        minMasked = min(maskedImg(:));
        maxMasked = max(maskedImg(:)) / 3;
        scaledMaskedImg = (maskedImg - minMasked) / (maxMasked - minMasked + eps); %#ok<NASGU>

        % Blue overlay where maskedImg == 0
        blueOverlay = zeros(size(maskedImg,1), size(maskedImg,2), 3);
        blueOverlay(:,:,3) = (maskedImg == 0);

        hOverlay = imshow(blueOverlay, 'Parent', hAxes);
        set(hOverlay, 'AlphaData', 0.25);

        hold off;
        title('Image with Transparent Overlay');
    end

    function rescaledMap = rescaleBrainOutline(origMap, rescaledSize)
        rescaledMap = imresize(origMap, rescaledSize, 'nearest');
    end

end





% function options = MaskVasculature_JP(ref_img_aligned, options)
% % ref_img_aligned: a reference image cropped and aligned (prepro_log.cropped_alligned_img).
% % options: prepro_log that contains the result from ManualAlignmentAdjust.
% 
% if options.mask_vasc %mask vasc
%     %Indentify vasculature by subtracting a large median filter  (e.g. smoothed image
%     %form the base image. This excentuates dark areas that are surrounded by
%     %light areas (e.g. vasculature).
%     smoothed = medfilt2(ref_img_aligned,[125 125]); %Use large neighboorhood e.g. 1/4 of size
%     mask = ref_img_aligned-smoothed;
% 
%     % Calculate the threshold to remove pixels lower than 'x' standard deviations below the mean
%     meanMask = nanmean(mask(:)); % Calculate the mean over the entire mask
%     stdMask = nanstd(mask(:)); % Calculate the standard deviation over the entire mask
% 
%     % Apply threshold: Keep pixels 'x' standard deviations below the mean
%     x = options.vasc_std; % 'x' standard deviations specified by user
%     threshold = meanMask - (x * stdMask); % Calculate the threshold value
%     mask = mask >= threshold; % Apply threshold
% 
%     %threshold to remove pxs lower that 'x' standard deviations below mean
%     %of Vascmask (e.g. the vasculature and other dark blemishes).
%     %mask = mask>=(nanmean(nanmean(mask))-(options.vasc_std*nanstd(nanstd(mask))));
% 
%     %close the image to clean up noise
%     se = strel('disk',options.close_disk_size);
%     mask = imclose(mask,se);
% 
%     %Add a premade mask to outline the brain.
%     if options.mask_brain_outline
%         temp = load(options.mask_brain_outline_dir);
%         if ~isequal(size(temp.brainoutline), size(ref_img_aligned))
%             % rescale the brainOutline to the reference image 
%             temp.brainoutline = rescaleBrainOutline(temp.brainoutline, size(ref_img_aligned)); 
%         end
% 
%         %Crop the mask to the size of the image. Assumes top left alligned
%         temp.brainoutline = temp.brainoutline(1:options.crop_h, 1:options.crop_w);
%         mask = mask+temp.brainoutline;
%         mask = mask==2;
%     end
% 
%     if options.verbose || options.manual_mask
%         ref_masked = ref_img_aligned.*mask;
%         showImgAndMask(ref_img_aligned, ref_masked);
%         hold on;
%         scatter(options.bregma(1),options.bregma(2),'*');
%         title('Cropped Img with Vascmask overlay');
%     end
% 
%     %Manually Nan regions and unmask others
%     if options.manual_mask
%         dlg = questdlg('Would you like to manually mask regions?', ...
%             'mask ROIS','Yes','No','No');
%         while 1
%             switch dlg
%                 case 'Yes' %mask additional region
%                     close all
%                     hfig1 = figure('name','Draw ROI around additional regions to NaN');
%                     [~, hBaseImage1, ~] = showImgAndMask(ref_img_aligned, ref_img_aligned.*mask, hfig1); % started using this (2/16/24) instead of imshowpair to visualize images brighter
%                     %imshowpair(ref_img_aligned,ref_img_aligned.*mask);
%                     manual_mask = impoly(gca);
%                     wait(manual_mask);
%                     manual_mask = manual_mask.createMask(hBaseImage1);
%                     mask(manual_mask)=0; close;
%                     hfig2 = figure('name','Current masked Image');
%                     showImgAndMask(ref_img_aligned, ref_img_aligned.*mask, hfig2);
%                 case 'No'; break %end dlg
%             end
%             choice = questdlg(sprintf('Would you like to mask additional regions'),...
%                 'mask ADDITIONAL ROIS',...
%                 'Yes','No','No');
%             switch choice; case 'Yes';case 'No'; break; end
%         end
%         dlg = questdlg('Would you like to manually UNmask regions?', ...
%             'mask ROIS','Yes','No','No');
%         while 1 %Manual UNmask regions
%             switch dlg
%                 case 'Yes' %unmask additional region
%                     close all
%                     hfig3 = figure('name','Draw ROI around additional regions to NaN');
%                     [~, hBaseImage3, ~] = showImgAndMask(ref_img_aligned, ref_img_aligned.*mask, hfig3);
%                     manual_mask = impoly(gca);
%                     wait(manual_mask);
%                     manual_mask = manual_mask.createMask(hBaseImage3); %his is 1 in the desired ROI
%                     mask(manual_mask)=1; close
%                     hfig4 = figure('name','Current masked Image');
%                     showImgAndMask(ref_img_aligned, ref_img_aligned.*mask, hfig4);
%                 case 'No'; break
%             end
%             choice = questdlg(sprintf('Would you like to unmask additional regions'),...
%                 'mask ADDITIONAL ROIS',...
%                 'Yes','No','No');
%             switch choice; case 'Yes';case 'No'; break; end
%         end
%     end
% 
%     if options.verbose || options.manual_mask
%         close all
%         showImgAndMask(ref_img_aligned, ref_img_aligned.*mask);
%         %clim([min(ref_img_aligned(:)) max(min(ref_img_aligned(:)), max(ref_img_aligned(:))/4)])
%         hold on;
%         scatter(options.bregma(1),options.bregma(2),'*');
%         title('FINAL cropped Img with Vascmask overlay');
%         choice = questdlg(sprintf('Would you like to flag this brain for imperfections?'),...
%             'Flag brain',...
%             'Yes','No','No');
%         switch choice
%             case 'Yes'
%                 options.flagged_brain = 1;
%             case 'No'
%                 options.flagged_brain = 0;
%         end
%     elseif options.maskVasculature==0 && options.maskBrainOutline == 1 %only apply exterior mask
%         mask = brainoutline;
%         if options.verbose
%             close all
%             showImgAndMask(ref_img_aligned, ref_img_aligned.*mask);
%             %clim([min(ref_img_aligned(:)) max(min(ref_img_aligned(:)), max(ref_img_aligned(:))/4)])
%             hold on;
%             scatter(options.bregma(1),options.bregma(2),'*');
%             title('FINAL cropped Img with masked Brain Overlay');
%         end
%     else %apply no mask
%         mask = ones(size(ref_img_aligned));
%     end
%     options.mask = mask;
%     options.masked_ref_img_aligned = ref_img_aligned.*mask;
% end
% 
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     function [hfig, hBaseImg, hOverlay] = showImgAndMask(img, maskedImg, varargin)
% 
%         if nargin < 3 || ~isa(varargin{1}, 'matlab.ui.Figure')
%             % Prepare the figure for displaying the images
%             hfig = figure;
%         else
%             hfig = varargin{1};
%         end
% 
% 
%         % Scale the main image for better visibility
%         minImg = min(img(:));
%         maxImg = max(img(:));
%         %maxImg = max(img(:)) / 3; % Enhance contrast by adjusting the max value (BRIGHTNESS)
%         scaledImage = (img - minImg) / (maxImg - minImg);
% 
%         hAxes = axes('Parent', hfig);
% 
%         hBaseImg = imshow(scaledImage, 'Parent', hAxes); % Display the scaled base image
%         hold on;
% 
%         % Assuming maskedImg contains areas from img multiplied by a binary mask
%         % Convert maskedImg to the same scale as scaledImage for consistent display
%         minMasked = min(maskedImg(:));
%         maxMasked = max(maskedImg(:)) / 3;
%         scaledMaskedImg = (maskedImg - minMasked) / (maxMasked - minMasked);
%         % Create an RGB overlay where 0 parts of maskedImg are green
%         %greenOverlay = zeros(size(scaledMaskedImg, 1), size(scaledMaskedImg, 2), 3); % Initialize with zeros
%         blueOverlay = zeros(size(scaledMaskedImg, 1), size(scaledMaskedImg, 2), 3); % Initialize with zeros
%         % Set the green channel to 1 (or another value for a different shade of green) where maskedImg is 0
%         %greenOverlay(:,:,2) = (maskedImg == 0);
%         blueOverlay(:,:,3) = (maskedImg == 0);
%         % Overlay maskedImg with transparency
%         % Note: Ensure scaledMaskedImg is not all zeros; otherwise, it won't display
%         hOverlay = imshow(blueOverlay, 'Parent', hAxes);
%         set(hOverlay, 'AlphaData', 0.25); % Adjust alpha for desired transparency
% 
%         hold off;
%         title('Image with Transparent Overlay');
%     end
% 
%     function rescaledMap = rescaleBrainOutline(origMap, rescaledSize)
%         rescaledMap = imresize(origMap, rescaledSize, 'nearest');
%     end
% 
% 
% end