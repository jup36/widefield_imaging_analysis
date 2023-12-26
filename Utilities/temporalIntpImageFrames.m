function newMatrix = temporalIntpImageFrames( imgFrameMat, originalTimePoints, targetTimePoints )


% Assuming tbytDffsm is your original 3D matrix
% And tbytDat(t).frameTrel and tbytDat(t).fraceCamTrel are your time vectors

%originalTimePoints = tbytDat(t).frameTrel(valDffFrameMinI:valDffFrameMaxI);
%targetTimePoints = tbytDat(t).faceCamTrel;

% Get the size of the original matrix
[height, width, ~] = size(imgFrameMat);

% Preallocate the new 3D matrix
newMatrix = nan(height, width, length(targetTimePoints));

% Interpolate for each pixel
for i = 1:height
    for j = 1:width
        pixelTimeSeries = squeeze(imgFrameMat(i, j, :));
        if sum(~isnan(pixelTimeSeries))==length(pixelTimeSeries)
            newMatrix(i, j, :) = interp1(originalTimePoints, pixelTimeSeries, targetTimePoints, 'linear', 'extrap');
        else
            newMatrix(i, j, :) = nan(length(targetTimePoints), 1); 
        end
    end
end


end