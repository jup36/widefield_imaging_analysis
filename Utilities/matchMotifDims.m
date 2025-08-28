function Wnew = matchMotifDims(Worg, nanpxs_org, nanpxs_new, opts) 
%This function adjusts the dimension of the vectorized spatiotemporal
% motifs to equalize the mismatch caused by different NaN pixels (nanpxs_org vs. nanpxs_new). 
fullPix = opts.originaldimensions(1)*opts.originaldimensions(2); 
Wfull = MaskTensor(Worg, nanpxs_org, [fullPix, size(Worg,2), size(Worg,3)]);
Wnew = Wfull(~ismember(1:fullPix, nanpxs_new), :, :); 
end