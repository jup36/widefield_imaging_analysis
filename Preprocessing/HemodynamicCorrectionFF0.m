function [dff, dff_b, dff_v] = HemodynamicCorrectionFF0(stack, opts)
%Camden MacDowell 2019
%Followed allen et al., 2017 neuron and Musall et al., 2019 Nature
%Neuro subtraction method for hemodynamic correction.

fprintf('Performing hemodynamic correction\n');
stack = SplitWavelengths(stack,opts.wavelength_pattern);

%Get hemo stack
stack_v = stack{opts.correction_wavelength};

%neural signal
stack(opts.correction_wavelength) = [];
stack_b = stack{1};

%Trim until both the same length (may be one frame discrepancy)
min_length = cellfun(@(x) size(x,3), {stack_b,stack_v},'UniformOutput',0);
stack_b = stack_b(:,:,1:min(cell2mat(min_length)));
stack_v = stack_v(:,:,1:min(cell2mat(min_length)));

%remove the masked pixels by setting to NaN
masked_pxls = (stack_b==0);
stack_b(masked_pxls)=NaN;
stack_v(masked_pxls)=NaN;

%there are large transients in both signals for the first 10 seconds of the recording, likely due to LED warm up (should add 'burn in' frames in the future);
%this can effects downstream results so replace the first 10 seconds of both with the average of the first minute; You should remove these at the end of analysis (after timelocking to behavior/ephys)
stack_b(:,:,1:3*opts.fps) = repmat(nanmean(stack_b(:,:,1:30*opts.fps),3),1,1,3*opts.fps);
stack_v(:,:,1:3*opts.fps) = repmat(nanmean(stack_v(:,:,1:30*opts.fps),3),1,1,3*opts.fps);

% Calculate the mean and standard deviation across frames (noticed an outlier frame, if any replace them!)
stack_b_frM = squeeze(nanmean(nanmean(stack_b, 1), 2));
stack_b_frM_mean = nanmean(stack_b_frM);
stack_b_frM_std = nanstd(stack_b_frM);

% Identify outlier frames (those beyond ±2 standard deviations from the mean)
outlier_frames_b = abs(stack_b_frM - stack_b_frM_mean) > 5 * stack_b_frM_std;

% Replace outlier frames with the mean of the first 30 * opts.fps frames
if any(outlier_frames_b)
    replacement_frame_b = nanmean(stack_b(:, :, 1:30 * opts.fps), 3);
    % Replicate the replacement_frame to match the number of outlier frames
    stack_b(:, :, outlier_frames_b) = repmat(replacement_frame_b, 1, 1, sum(outlier_frames_b));
end

% Calculate the mean and standard deviation across frames (noticed an outlier frame, if any replace them!)
stack_v_frM = squeeze(nanmean(nanmean(stack_v, 1), 2));
stack_v_frM_mean = nanmean(stack_v_frM);
stack_v_frM_std = nanstd(stack_v_frM);

% Identify outlier frames (those beyond ±2 standard deviations from the mean)
outlier_frames_v = abs(stack_v_frM - stack_v_frM_mean) > 5 * stack_v_frM_std;

% Replace outlier hemo frames with the mean of the first 30 * opts.fps frames
if any(outlier_frames_v)
    replacement_frame_v = nanmean(stack_v(:, :, 1:30 * opts.fps), 3);
    % Replicate the replacement_frame to match the number of outlier frames
    stack_v(:, :, outlier_frames_v) = repmat(replacement_frame_v, 1, 1, sum(outlier_frames_v));
end

%pixelwise hemodynamic correction
[nX,nY,nZ] = size(stack_b);
[stack_b, bad_col] = conditionDffMat(stack_b); % input dim: pixel x pixel x frames, output dim: frames x pixels(excluding bad ones)
[stack_v, ~] = conditionDffMat(stack_v);
stack_v_smoothed = NaN(size(stack_v));
for i = 1:size(stack_v_smoothed,2)
    temp = stack_v(:,i);
    %smooth with ~450ms gaussian
    temp = smoothdata(temp,'gaussian',floor(opts.fps/2)); %opts.fps is per wavelength (if multiplexed)
    %least square linear regression to find scaling factor and intercept
    coef = polyfit(temp,stack_b(:,i),1);
    correction = @(x) coef(1)*x+coef(2);
    stack_v_smoothed(:,i) = correction(temp);
end

%remake into full size
stack_b = conditionDffMat(stack_b,bad_col,[],[nX,nY,nZ]);
stack_v_smoothed = conditionDffMat(stack_v_smoothed,bad_col,[],[nX,nY,nZ]);

%calculate the dff for each (as in wekselblatt, 2017)
[dff_v,~] = makeDFF(stack_v_smoothed, opts, 'fractional');
[dff_b,~] = makeDFF(stack_b, opts, 'fractional');

%Division correction
dff = (dff_b./dff_v-1).*100;

end



% %     Smooth V to remove high frequency non hemodynamic activity
% %
% %
% %     least square linear regression to find scaling factor and intercept
% %     coef = polyfit(stack_v(~isnan(stack_v)),stack_b(~isnan(stack_b)),1);
% %
% %     correction equation
% %     correction = @(x) coef(1)*x+coef(2);
% %     stack_v_smoothed = reshape(correction(stack_v(:)),size(stack_v));
% %
% %     confirm that you scaled the correct direction
% %     err = nanmean(stack_v_smoothed(:)-stack_b(:));
% %     err2 = nanmean(stack_v(:)-stack_b(:));
% %
% %     calculate the dff for each
% %     [dff_v,~] = makeDFF(stack_v_smoothed, opts);
% %     [dff_b,~] = makeDFF(stack_b, opts);
% %
% %     Subtraction correction DFF = dff_b-dff_v;
% %     dff = dff_b-dff_v;
% %
