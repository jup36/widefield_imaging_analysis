function overlayMaskOutline(img, mask)
    figure; 
    
    imagesc(img);
    axis image;
    set(gca,'YDir','reverse');
    colormap turbo;
    colorbar;
    hold on;

    B = bwboundaries(mask);
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5);
    end
  

end