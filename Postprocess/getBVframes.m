function [stack_b, stack_v, bFrameI, vFrameI] = getBVframes(stack)

numbFrames = size(stack, 3);

stack_b = stack(:, :, 1:2:end);
bFrameI = 1:2:numbFrames;

stack_v = stack(:, :, 2:2:end);
vFrameI = 2:2:numbFrames;

%Trim until both the same length (may be one frame discrepancy)
min_length = min(cell2mat(cellfun(@(x) size(x,3), {stack_b,stack_v},'UniformOutput',0)));
stack_b = stack_b(:,:,1:min_length);
bFrameI = bFrameI(1:min_length);

stack_v = stack_v(:,:,1:min_length);
vFrameI = vFrameI(1:min_length);

% to ensure correct assignment of blue and violet frames check their mean intensities
if nanmean(stack_b(:)) < nanmean(stack_v(:))
    stack_b_copy = stack_b;
    bFrameI_copy = bFrameI;
    stack_b = stack_v;
    bFrameI = vFrameI;
    stack_v = stack_b_copy;
    vFrameI = bFrameI_copy;
    clearvars stack_b_copy
end

%remove the masked pixels by setting to NaN
masked_pxls = (stack_b==0);
stack_b(masked_pxls)=NaN;
stack_v(masked_pxls)=NaN;

end