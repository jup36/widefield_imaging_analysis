function [dffsOnTfItp1, dffsOnTfItp2, tF] = alignTwoDffTimeSeries(dff1, dff2, frameT1, frameT2, timeW)
% lickBouts = tbytDat(t).LickBoutRel;
% dffs = dffM1;
% frameT = tbytDat(t).frameTrel;
% timeW = [-0.9 5];

% get licks on the specified time frame
tF = timeW(1):0.001:timeW(end); % time frame

% get dff on the specified time frame with interpolation
ftI1 = frameT1 >= min(tF) & frameT1 <= max(tF);
dffsOnTf1 = dff1(ftI1);
frameTOnTf1 = frameT1(ftI1);

ftI2= frameT2 >= min(tF) & frameT2 <= max(tF);
dffsOnTf2 = dff2(ftI2);
frameTOnTf2 = frameT2(ftI2);

% Define the valid range based on frameTOnTf
validIdx = tF >= max(min(frameTOnTf1), min(frameTOnTf2)) & tF <= min(max(frameTOnTf1), max(frameTOnTf2));
tF = tF(validIdx); 

% Perform strict interpolation (no extrapolation)
dffsOnTfItp1 = interp1(frameTOnTf1, dffsOnTf1, tF, 'linear');
dffsOnTfItp2 = interp1(frameTOnTf2, dffsOnTf2, tF, 'linear');

end