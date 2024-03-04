function [stack_b, stack_v_corrected] = getPixelwiseScaledV(stack_b, stack_v, opts)
%This function uses polyfit to 'scale' the stack_v relative to stack_b.   
% badcols and badrows are list of indices
% of cols and rows removed from time x pxls format (e.g. nan or zero from vasculature
% and mask)

% pixelwise hemodynamic correction
[nX,nY,nZ] = size(stack_b);
[stack_b, bad_col] = conditionDffMat(stack_b); % input dim: pixel x pixel x frames, output dim: frames x pixels(excluding bad ones)
[stack_v, ~] = conditionDffMat(stack_v);
stack_v_corrected = NaN(size(stack_v));
for jj = 1:size(stack_v_corrected,2)
    temp = stack_v(:,jj);

    %smooth with ~450ms gaussian
    temp = smoothdata(temp,'gaussian',floor(opts.fps/2)); %opts.fps is per wavelength (if multiplexed)
    %least square linear regression to find scaling factor and intercept
    coef = polyfit(temp, stack_b(:,jj),1);
    correction = @(x) coef(1)*x+coef(2);
    stack_v_corrected(:,jj) = correction(temp);
end

%remake into full size
stack_b = conditionDffMat(stack_b,bad_col,[],[nX,nY,nZ]);
stack_v_corrected = conditionDffMat(stack_v_corrected,bad_col,[],[nX,nY,nZ]);
end