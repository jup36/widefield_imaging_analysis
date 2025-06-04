function [hemiMaskL, hemiMaskR] = hemiMaskLR(fullMask)
%This function takes a 2D mask of mouse dorsal cortex to split it into
% left and right halves.

hemiMaskR = fullMask;
hemiMaskR(:, 1:floor(size(fullMask, 2)/2))=false;

hemiMaskL = fullMask;
hemiMaskL(:, ceil(size(fullMask, 2)/2):end)=false;
end