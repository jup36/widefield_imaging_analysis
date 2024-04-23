function [stack_b, stack_v, bFrameI, vFrameI] = getBVframesWithIndexing(stack, varargin)

numbFrames = size(stack, 3);

if nargin == 1
    bFrameI = 1:2:numbFrames;
    vFrameI = 2:2:numbFrames;
elseif nargin == 2 % if there's an optional user-defined bv index provided, use that index
    bvIndex = varargin{1}; 
    if sum(~isnan(bvIndex)) == length(bvIndex)
        bFrameI = find(bvIndex); 
        vFrameI = find(~bvIndex); 
    else
        bFrameI = 1:2:numbFrames;
        vFrameI = 2:2:numbFrames; 
    end
end

stack_b = stack(:, :, bFrameI); 
stack_v = stack(:, :, vFrameI); 

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