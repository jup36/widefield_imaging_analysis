function Wnew = matchMotifDims(Worg, nanpxs_org, nanpxs_new, opts) 
%This function adjusts the dimension of the vectorized spatiotemporal
% motifs to equalize the mismatch caused by different NaN pixels (nanpxs_org vs. nanpxs_new). 
fullPix = opts.originaldimensions(1)*opts.originaldimensions(2); 
Wfull = MaskTensor(Worg, nanpxs_org, [fullPix, size(Worg,2), size(Worg,3)]);
Wnew = Wfull(~ismember(1:fullPix, nanpxs_new), :, :); 
end

function Y = MaskTensor(X,nanpxs,dims)
%camden macdowell - timeless

if size(X,1) < dims(1)
    Y = zeros(dims); 
    temp = ones(dims(1),1);
    temp(nanpxs)=0;
    Y(temp==1,:,:) = X;
else
    Y = X;
    Y(nanpxs,:,:) = [];
end

end
