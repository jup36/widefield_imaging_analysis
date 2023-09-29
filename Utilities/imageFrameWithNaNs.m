function figHandle = imageFrameWithNaNs(frame, climits)
% climits = [-3 3];
% Create a new figure and store its handle
figHandle = figure;

% Create an adjusted colormap with white at the beginning
cmap = [0, 0, 0; hot(256)];

% Display the frame
imagesc(frame);

% Use the custom colormap
colormap(cmap);

% Set NaN values to a value outside the data range for visualization
clim(climits);

% Show the colorbar (optional)
colorbar;
end
