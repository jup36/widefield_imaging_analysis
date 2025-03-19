function h = montageMotif(motif_w, nanpxs)
% motif_w: pixels-by-timelag (expected)

if length(size(motif_w))==2
    if size(motif_w,1)>size(motif_w,2)
        motif_w = motif_w';
    end
end

motif_w_rs = conditionDffMat(motif_w, nanpxs); 

h = figure; 
montage(motif_w_rs, 'Size', [1,size(motif_w_rs, 3)],'DisplayRange', [0 0.5]); colormap turbo; colorbar %change 4th dimension for each motif

end